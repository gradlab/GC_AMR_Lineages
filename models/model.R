source("models/prep_data.R")

#Specify a model for the given dataset
lineage_model <- function(lin_phys, lin_data, covars, tname, N_grid_per_unit = 4, ...)
{
    t_end <- max(covars[[tname]])
    stopifnot(abs(t_end - round(t_end)) < 1e-8)
    
    n_lin <- length(lin_phys)
    stopifnot(nrow(lin_data)==n_lin)

    covars <- covars %>% 
        arrange(.data[[tname]]) %>%
        relocate(all_of(tname)) %>%
        as.matrix() %>%
        t() 

    rn_cov <-  rownames(covars[-1, ])

    grid_data <- make_grid(t_end, covars, t_end * N_grid_per_unit)
    tree_data <- make_tree_data(lin_phys, lin_data$mrst, grid_data$grid_endpts)

    structure(
        list(
            n_lin = n_lin,
            n_interactions = 0L,
            lin_interactions = list(),
            lin_data = lin_data,
            tree_data = tree_data,
            grid_data = grid_data,
            rn_cov = rn_cov
        ),
        class = "linmod"
    )
}

add_site_interaction <- function(x, site_vars, covs=NULL, sum_covs=FALSE, re=NULL, site_name=NULL)
{
    li <- lin_rt_interaction(site_vars, covs, sum_covs, re, site_name)
    x + li 
}

#Specify an interaction between a given site/locus and the associated determinants found there
#and the a subset of the treatment covariates
lin_rt_interaction <- function(site_vars, covs, sum_covs, re=NULL, site_name=NULL)
{
    structure ( 
        list(site_vars=site_vars, covs=covs, sum_covs=sum_covs, re=re, site_name=site_name, type="interaction"),
        class = "linmod_elem"
    )
}

`+.linmod` <- function(x, y)
{
    stopifnot(class(y)=="linmod_elem")
    if (y$type == "interaction")
    {
        x$n_interactions <- x$n_interactions + 1L
        if(is.null(y$site_name)) 
        {
            y$site_name <- x$n_interactions
        }

        cov_names <- y$covs
        vnames <- y$site_vars

        stopifnot(all(cov_names %in% x$rn_cov))
        stopifnot(all(vnames %in% colnames(x$lin_data)))

        dat <- as.matrix(x$lin_data[, vnames])
        
        y$design_matrix <- dat
        y$contr_names <- vnames
        y$cov_idx <- sapply(cov_names, function(z) which(x$rn_cov == z))

        y$has_re <- FALSE

        if(!is.null(y$re)) 
        {
            stopifnot(y$sums_covs == F)
            stopifnot(y$re %in% colnames(x$lin_data))
            re_factor <- as.factor(x$lin_data[[y$re]])

            y$has_re <- TRUE

            y$re_levels <- levels(re_factor)
            y$level_idx <- as.integer(re_factor)
            y$level_idx <- ifelse(is.na(y$level_idx), 0L, y$level_idx)
            y$N_levels <- length(y$re_levels)       
        }
        
        x$lin_interactions[[x$n_interactions]] <- y
    }
   
    x
}

#Generate stan files for a given model
generate_model <- function(x)
{
    foo_site_data <- function(y)
    {
        cov_names <- y$covs
        grp_sz <- length(y$contr_names)
        
        data.frame(
            n=grp_sz,
            m=length(cov_names),
            site=y$site_name
        )
    }

    foo_site_eff <- function(y)
    {
        cov_names <- y$covs
        grp_sz <- length(y$contr_names)

        data.frame(
            n=grp_sz,
            m=length(cov_names),
            site=y$site_name,
            idx_b=seq_along(y$cov_idx),
            idx_cov=y$cov_idx
        )
    }

    foo_re_data <- function(y)
    {
        cov_names <- y$covs
        grp_sz <- length(y$contr_names)
        N_lvls <- y$N_levels  
        cov_names <- y$covs
        
        data.frame(
            n=grp_sz,
            m=length(cov_names),
            N_lvls=N_lvls,
            site=y$site_name
        )
    }

    foo_re_eff <- function(y)
    {
        cov_names <- y$covs
        grp_sz <- length(y$contr_names)
        N_lvls <- y$N_levels  
        
        data.frame(
            n=grp_sz,
            N_lvls=N_lvls,
            site=y$site_name,
            idx_b=seq_along(y$cov_idx),
            idx_cov=y$cov_idx
        )
    }

    site_data <- do.call(rbind, lapply(
        x$lin_interactions, foo_site_data
    ))

    site_eff_data <- do.call(rbind, lapply(
        x$lin_interactions, foo_site_eff
    ))

    site_presence_absence <-lapply(
        x$lin_interactions, function(y) y$design_matrix
    )

    re_site_idx <- list()
    re_Z <- list()

    idx_names <- c()
    Z_names <- c()

    re_data <- data.frame()
    re_eff <- data.frame()

    for (y in x$lin_interactions){
        if (y$has_re)
        {
            re_site_idx <- c(re_site_idx, list(y$level_idx))
            re_Z <- c(re_Z, list(y$design_matrix))

            idx_names <- c(idx_names, paste0("site_", y$site_name, "_re_idx"))
            Z_names <- c(Z_names, paste0("site_", y$site_name, "_Z"))

            re_data <- rbind(re_data, foo_re_data(y))
            re_eff <- rbind(re_eff, foo_re_eff(y))
        }
    }

    names(re_site_idx) <- idx_names 
    names(re_Z) <- Z_names 


    names(site_presence_absence) <- sapply(x$lin_interactions, function(y) paste0("site_",y$site_name,"_presence_absence"))

    weighted_inter <- sapply(x$lin_interactions, function (y) y$sum_covs)
    model_txt <- make_model(site_data, 
        site_eff_data, 
        site_data[!weighted_inter, ], 
        site_data[weighted_inter, ])

    model_txt <- make_re(model_txt, re_data, re_eff)
        
    cat(model_txt,file="mod.stan")

    data_list <- c(x$tree_data, x$grid_data, site_presence_absence, re_site_idx, re_Z)
    mod <- cmdstan_model("mod.stan")
    fit <- NA

    return(list(lineage_model = x,
    data_list = data_list,
    model_txt = model_txt,
    sample = function(...)
    {
        fit <<- mod$sample(data=data_list, 
            ...
        )
    },
    fit = function()
    {
        fit
    }))
}

make_re <- function(skeleton, re_data, re_eff)
{

    re_grp_block <- make_block(re_data,
        "array[N_lin] int<lower=0> site_{site}_re_idx;",
        "matrix[N_lin, {n}] site_{site}_Z;"  
    )

    re_block <- make_block(re_data,
        "array[{n}, {N_lvls}, {m}] real site_{site}_re_b_raw;",
        "matrix[{n}, {N_lvls}] site_{site}_re_a_raw;",
        "real <lower=0> site_{site}_re_sd_b;",
        "real <lower=0> site_{site}_re_sd_a;"
    )

    re_transf_block <- make_block(re_data,
        "array[{n}] matrix[{N_lvls}, {m}] site_{site}_re_b = scale_re_b(site_{site}_re_b_raw, site_{site}_re_sd_b * 0.25);",
        "matrix[{n}, {N_lvls}] site_{site}_re_a = site_{site}_re_sd_a * 0.25 * site_{site}_re_a_raw;"
    )

    re_lin_block <- make_block(re_data,
        "matrix [N_lin, {m}] re_{site}_lin_b = lin_re_b(site_{site}_re_b, site_{site}_Z, site_{site}_re_idx);",
        "vector [N_lin] re_{site}_lin_a = lin_re_a(site_{site}_re_a, site_{site}_Z, site_{site}_re_idx);"
    )

    re_eff_block <- c(make_block(re_data,
        "lin_trends += rep_matrix(re_{site}_lin_a, N_grid);"
    ),
    make_block(re_eff,
        "lineage_coeffs[:, {idx_cov}] += re_{site}_lin_b[:, {idx_b}];"
    ))

    re_model_block <- make_block(re_data,
        "to_array_1d(site_{site}_re_a_raw) ~ normal(0, 1);",
        "to_array_1d(site_{site}_re_b_raw) ~ normal(0, 1);",
        "site_{site}_re_sd_a ~ normal(0, 4);",
        "site_{site}_re_sd_b ~ normal(0, 4);"
    )

    skeleton <- str_replace(skeleton, "###```RE_GRP_BLOCK```",glue_collapse(re_grp_block,"\n"))
    skeleton <- str_replace(skeleton, "###```RE_BLOCK```",glue_collapse(re_block,"\n"))
    skeleton <- str_replace(skeleton, "###```RE_TRANSF_BLOCK```",glue_collapse(re_transf_block,"\n"))
    skeleton <- str_replace(skeleton, "###```RE_LIN_BLOCK```",glue_collapse(re_lin_block,"\n"))
    skeleton <- str_replace(skeleton, "###```RE_EFF_BLOCK```",glue_collapse(re_eff_block,"\n"))
    skeleton <- str_replace(skeleton, "###```RE_MODEL_BLOCK```",glue_collapse(re_model_block,"\n"))

    return(skeleton)
}

make_model <- function(site_data, site_eff_data, unweighted_data, weighted_data)
{
    data_block <- c(make_block(site_data,
        "matrix [N_lin, {n}] site_{site}_presence_absence;"
    ))

    intercept_block <- c(make_block(site_data,
        "vector[{n}] site_{site}_a_ncp;"
    ))
    
    b_raw_block <- c(make_block(unweighted_data,
        "matrix [{n}, {m}] site_{site}_b_ncp;"),
    make_block(weighted_data,
        "matrix [{n}, {m}] site_{site}_b_ncp;",
        "vector [{n}] site_{site}_b_eff;",
        "real <lower=0, upper=1> site_{site}_pool;"
    ))

    b_transf_block <- c(make_block(weighted_data,
        "matrix [{n}, {m}] site_{site}_b = site_{site}_pool * rep_matrix(site_{site}_b_eff * 0.25, {m}) + (1.0 - site_{site}_pool) * site_{site}_b_ncp * 0.25;"
    ),
    make_block(unweighted_data,
        "matrix [{n}, {m}] site_{site}_b = 0.25 * site_{site}_b_ncp;"
    ))

    a_transf_block <- make_block(site_data,
        "vector [{n}] site_{site}_a = 0.25 * site_{site}_a_ncp;")

    tparams_block <- make_block(site_data,
        "vector[N_lin] site_{site}_lineage_a = site_{site}_presence_absence * site_{site}_a;",
        "matrix[N_lin, {m}] site_{site}_lineage_b = site_{site}_presence_absence * site_{site}_b;"
    )

    effects_block <- c(make_block(site_data,
        "lin_trends += rep_matrix(site_{site}_lineage_a, N_grid);"
    ),
    make_block(site_eff_data,
        "lineage_coeffs[:, {idx_cov}] += site_{site}_lineage_b[:, {idx_b}];"
    ))

    model_block <- c(make_block(site_data,
        "site_{site}_a_ncp ~ normal(0, 4);"
    ),
    make_block(unweighted_data,
        "to_array_1d(site_{site}_b_ncp) ~ normal(0, 4);"
    ),
    make_block(weighted_data,
        "site_{site}_b_eff ~ normal(0, 4);",
        "to_array_1d(site_{site}_b_ncp) ~ normal(0, 4);"
    ))

    generated_block <- c(make_block(site_data,
        "lineage_intercepts += site_{site}_lineage_a;"))


    skeleton <- read_file("models/lineage_mod_skeleton.stan")

    skeleton <- str_replace(skeleton, "###```DATA_BLOCK```",glue_collapse(data_block,"\n"))
    skeleton <- str_replace(skeleton, "###```SITE_INTERCEPT_BLOCK```",glue_collapse(intercept_block,"\n"))
    skeleton <- str_replace(skeleton, "###```B_RAW_BLOCK```",glue_collapse(b_raw_block,"\n"))
    skeleton <- str_replace(skeleton, "###```B_TRANSF_BLOCK```",glue_collapse(b_transf_block,"\n"))
    skeleton <- str_replace(skeleton, "###```A_TRANSF_BLOCK```",glue_collapse(a_transf_block,"\n"))
    skeleton <- str_replace(skeleton, "###```TPARAMS_BLOCK```",glue_collapse(tparams_block,"\n"))
    skeleton <- str_replace(skeleton, "###```EFFECTS_BLOCK```",glue_collapse(effects_block,"\n"))
    skeleton <- str_replace(skeleton, "###```MODEL_BLOCK```",glue_collapse(model_block,"\n"))
    skeleton <- str_replace(skeleton, "###```GENERATED_BLOCK```",glue_collapse(generated_block,"\n"))

    return(skeleton)
}

make_block <- function(x, ..., indent_lvl=4)
{
    stub <- list(...)
    indent <- str_dup(" ", indent_lvl)
    stub <- sapply(stub, function(s) paste0(indent, s))
    stub <- paste0(stub, collapse="\n")
    glue_collapse(glue_data(x, stub, .trim=F),"\n")
}
library(ggridges)
library(patchwork)

name_transf <- function(x) str_replace_all(x, "_", " ")

#Compute predicted growth rate effects
comp_effects <- function(fit, site, select_motifs=c(), exclude_motifs=c(), exclude_levels=c())
{
    sname <- fit$lineage_model$lin_interactions[[site]]$site_name
    contr_names <- fit$lineage_model$lin_interactions[[site]]$contr_names
    cov_names <- fit$lineage_model$lin_interactions[[site]]$covs
    cov_idx <- fit$lineage_model$lin_interactions[[site]]$cov_idx

    has_re <- fit$lineage_model$lin_interactions[[site]]$has_re
    ddf <- as_draws_df(fit$fit()$draws())

    exclude_motifs <- name_transf(exclude_motifs)
    select_motifs <- name_transf(select_motifs)

    n_cov <- length(cov_names)

    eff_df <- do.call(rbind, lapply(
        seq_along(contr_names),
        function(x) 
        {
            data.frame(
                motif=name_transf(contr_names[x]),
                level="Mean",
                cov=cov_names,
                cov_idx=cov_idx,
                vname=paste0("site_",sname,"_b[",x,",",seq_along(cov_names),"]"),
                re_name=NA,
                a_vname=paste0("site_",sname,"_a[",x,"]"),
                a_re_name=NA
            )
        }
    ))

    if (has_re) 
    {
        re_levels <- fit$lineage_model$lin_interactions[[site]]$re_levels
        re <- fit$lineage_model$lin_interactions[[site]]$re

        eff_df_re <- do.call(rbind, lapply(
            seq_along(contr_names),
            function(x) 
            {
                do.call(rbind, lapply(
                    seq_along(re_levels),
                    function(y)
                    {
                        data.frame(
                            motif=name_transf(contr_names[x]),
                            level=paste0(re, " ", re_levels[y]),
                            cov=cov_names,
                            cov_idx=cov_idx,
                            vname=paste0("site_",sname,"_b[",x,",",seq_along(cov_names),"]"),
                            re_name=paste0("site_",sname,"_re_b[",x,",",y,",",seq_along(cov_names),"]"),
                            a_vname=paste0("site_",sname,"_a[",x,"]"),
                            a_re_name=paste0("site_",sname,"_re_a[",x,",",y,"]")
                        )
                    }
                ))
            }
        ))

        eff_df <- rbind(eff_df, eff_df_re)
    }

    if (length(select_motifs) > 0)
    {
        eff_df <- eff_df %>% 
            filter(motif %in% select_motifs)
    }

    eff_df <- eff_df %>% 
        filter(!(motif %in% exclude_motifs)) %>% 
        filter(!(level %in% exclude_levels))

    grid_pts <- fit$data_list$grid_endpts
    grid_pts_to_keep <- which(abs(grid_pts %% 1) < 1e-8)
    if(grid_pts_to_keep[1] != 1)
    {
        grid_pts_to_keep <- c(1L, grid_pts_to_keep)
    }

    covs <- fit$data_list$covars[, grid_pts_to_keep]

    do.call(rbind, apply(
        eff_df, 1,
        function(x)
        {
            icept <- (unlist(ddf[, x[7]]) + if(is.na(x[8])) 0 else unlist(ddf[, x[8]])) / n_cov
            coeff <- unlist(ddf[, x[5]]) + if(is.na(x[6])) 0 else unlist(ddf[, x[6]])
            
            (matrix(coeff, nrow=length(coeff), ncol=1) %*% covs[as.integer(x[[4]]), ]) %>% 
            apply(2, function(rv) rv + icept) %>%
                as_tibble() %>%
                mutate(draw=1:nrow(ddf), 
                    motif=x[1], 
                    level=x[2], 
                    cov=x[3])# %>%
                #pivot_longer(!c(draw, motif, level, cov), 
                #    names_to="grid", 
                #    names_pattern = "(\\d+)", 
                #    names_transform = as.integer)
        }
    )) %>% 
        pivot_longer(!c(draw, motif, level, cov), 
            names_to="grid", 
            names_pattern = "(\\d+)", 
            names_transform = as.integer) %>%
        mutate(grid=grid_pts_to_keep[grid]) %>%
        group_by(draw, motif, level, grid) %>%
        summarise(y=sum(value), .groups = "keep")
}

summary_95 <- function(x)
{ 
    x <- unname(quantile(x, c(.025, .975)))
    return(data.frame(xmin=x[1], xmax=x[2]))
}

summary_50 <- function(x)
{
    x <- unname(quantile(x, c(.25, .75)))
    return(data.frame(xmin=x[1], xmax=x[2]))
}

make_dens_df <- function(data_grouped, value_col="value")
{
    dens <- data_grouped %>% 
        group_modify(
            ~ {
                d <- density(.x[[value_col]], na.rm=T)
                tibble(x=d$x, y=d$y/max(d$y))
            }
        )
    dens
}

plot_density_ridge <- function(data, ...) 
{
    dens <- make_dens_df(data)

    # Get credible interval width and median
    cred.int <- data %>% 
        summarise(n=list(summary_50(value)),
                w=list(summary_95(value)),
                m=median(value, na.rm=TRUE),
                y=0, .groups = "keep") %>% 
        unnest_wider(c(n, w),names_sep="_") %>%
        mutate(med_xlo = m-0.05*(n_xmax-n_xmin), med_xhi = m+0.05*(n_xmax-n_xmin))
    
    dens_med <- dens %>%
        left_join(cred.int %>% select(!y)) %>%
        filter(x>=med_xlo & x<=med_xhi) %>%
        select(level,x,y)

    dens %>% 
        ggplot(aes(x=x, height=y, y=level)) +
        geom_vline(xintercept=0, colour="grey50", linetype = "dashed", size=1) +
        geom_ridgeline(fill=hcl(265, 40, 60), alpha=.6, scale=0.8, ...) + 
        geom_ridgeline(data=dens_med, aes(x=x, height=y, y=level), fill="gray35", color="gray35",scale=0.8, ...) +
        geom_linerange(data=cred.int, aes(xmin=w_xmin, x=m, xmax=w_xmax, y=level), color="gray60",lwd=2) +
        geom_linerange(data=cred.int, aes(xmin=n_xmin, x=m, xmax=n_xmax, y=level), color="gray35", lwd=2) +
        labs(y="", x="Value") + 
        scale_x_continuous(n.breaks = 10) + 
        theme_bw() +
        theme(
            axis.text.x=element_text(size=rel(0.7), angle = 45, hjust=1),
            plot.margin = margin(0, 0, 0, 0, "cm"),
            legend.position="none",
            panel.grid.major.y = element_blank(), 
            panel.grid.minor.y = element_blank(),
            axis.line = element_line(size=rel(0.2), colour = "grey90"),
            plot.title = element_text(hjust = 0.5,size=rel(1.0)), ...)
}

#Comparison plots for different determinants at the same locus
plot_site_diff <- function(fit, site, offs=0, motifs_to_compare=c(), exclude_levels=c(), ci_level=c(.025, .5, .975))
{
    n_grid <- fit$data_list$N_grid
    grid_ts <- offs - fit$data_list$grid_endpts

    eff_data <- comp_effects(fit, site, c(), c(), exclude_levels) 
    
    eff_data <- eff_data %>% 
        filter(motif %in% motifs_to_compare) %>%
        mutate(y=ifelse(motif == motifs_to_compare[1], y, -y)) %>%
        group_by(level, grid, draw) %>%
        summarise(y=sum(y), motif="Difference", .groups="keep") %>% 
        ungroup() 
    #draws <- sample(max(eff_data$draw), 5)
    #
    #eff_samp <- eff_data %>% 
    #    filter(draw %in% draws) %>%
    #    mutate(x=grid_ts[as.integer(grid)])

    eff_data %>%
        group_by(motif, level, grid) %>% 
        summarise(as_tibble_row(quantile(y, ci_level),
            .name_repair = \(x) paste0('q', parse_number(x))),.groups="keep") %>% 
        mutate(x=grid_ts[as.integer(grid)]) %>%
        ggplot(aes(x=x, y=q50)) + 
        geom_line(lwd=1.5) + 
        #geom_line(data=eff_samp, aes(group=interaction(draw),x=x,y=y), alpha=.4)+
        geom_ribbon(aes(ymin = q2.5, ymax = q97.5), fill="grey35", alpha=.5) +
        labs(y= "Predicted Growth Rate Difference", x="Year") +
        geom_hline(yintercept = 0, linetype="dashed")+
        facet_grid(rows=vars(level), cols=vars(motif)) + 
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        panel.grid.minor.x=element_blank())
}

#Plot predicted effects for a given site
plot_site_effects <- function(fit, site, offs=0, relative=FALSE, select_motifs=c(), exclude_motifs=c(), exclude_levels=c(), combine_motifs=c(), ci_level=c(.025, .5, .975))
{

    n_grid <- fit$data_list$N_grid
    grid_ts <- offs - fit$data_list$grid_endpts

    eff_data <- comp_effects(fit, site, select_motifs, exclude_motifs, exclude_levels) 
    if (length(combine_motifs)>0)
    {
        mname <- paste0(combine_motifs, collapse=" + ")

        eff_data <- eff_data %>% 
            filter(motif %in% combine_motifs) %>%
            group_by(level, grid, draw) %>%
            summarise(y=sum(y), motif=mname, .groups="keep") %>% 
            ungroup() 
    }

    has_re <- fit$lineage_model$lin_interactions[[site]]$has_re

    if (relative)
    {
        eff_data <- eff_data %>% 
            group_by(motif, level, draw) %>%
            mutate(y=y-y[which.min(grid)])
    }

    plt <- eff_data %>%
        group_by(motif, level, grid) %>% 
        summarise(as_tibble_row(quantile(y, ci_level),
            .name_repair = \(x) paste0('q', parse_number(x))),.groups="keep") %>% 
        mutate(x=grid_ts[as.integer(grid)]) %>%
        ggplot(aes(x=x, y=q50)) + 
        geom_line(lwd=1.5) + 
        geom_ribbon(aes(ymin = q2.5, ymax = q97.5), fill="grey35", alpha=.5) +
        labs(y="Predicted Growth Rate Effect", x="Year") +
        geom_hline(yintercept = 0, linetype="dashed") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
            panel.grid.minor.x=element_blank())

    if(has_re) 
    {
        plt + facet_grid(cols=vars(motif), rows=vars(level))
    } else 
    {
        plt + facet_grid(cols=vars(motif), rows=NULL)
    }
}

#Plot densities for intercepts
plot_site_intercepts <- function(fit, site, exclude_motifs=c(), exclude_levels=c())
{
    sname <- fit$lineage_model$lin_interactions[[site]]$site_name
    contr_names <- fit$lineage_model$lin_interactions[[site]]$contr_names

    has_re <- fit$lineage_model$lin_interactions[[site]]$has_re
    plts <- NULL

    ddf <- as_draws_df(fit$fit()$draws())

    a_df <- do.call(rbind, lapply(
        seq_along(contr_names),
        function(x) 
        {
            data.frame(
                motif=name_transf(contr_names[x]),
                level="Mean",
                vname=paste0("site_",sname,"_a[",x,"]"),
                re_name=NA
            )
        }
    ))

    if (has_re) 
    {
        re_levels <- fit$lineage_model$lin_interactions[[site]]$re_levels
        re <- fit$lineage_model$lin_interactions[[site]]$re
        
        a_df_re <- do.call(rbind, lapply(
            seq_along(contr_names),
            function(x) 
            {
                do.call(rbind, lapply(
                    seq_along(re_levels),
                    function(y)
                    {
                        data.frame(
                        motif=name_transf(contr_names[x]),
                        level=paste0(re, " ", re_levels[y]),
                        vname=paste0("site_",sname,"_a[",x,"]"),
                        re_name=paste0("site_",sname,"_re_a[",x,",",y,"]")
                        )
                    }
                ))
            }
        ))
        a_df <- rbind(a_df, a_df_re)
    }

    a_data <- do.call(rbind, apply(
        a_df, 1,
        function(x)
        {
            data.frame(motif=x[1], 
                level=x[2], 
                value=unlist(ddf[, x[3]]) + if(is.na(x[4])) 0 else unlist(ddf[, x[4]]),
                draw=1:nrow(ddf)
                )
        }
    )) %>% 
        filter(!(motif %in% exclude_motifs)) %>% 
        filter(!(level %in% exclude_levels)) %>%
        group_by(motif, level)             

    plts <- a_data %>%
        plot_density_ridge() + 
        facet_grid(rows=vars(motif),scales="free_y")

    plts
}

#Plot densities for coefficients
plot_site_coeffs <- function(fit, site, exclude_motifs=c(), exclude_levels=c())
{
    sname <- fit$lineage_model$lin_interactions[[site]]$site_name
    contr_names <- fit$lineage_model$lin_interactions[[site]]$contr_names
    cov_names <- fit$lineage_model$lin_interactions[[site]]$covs

    has_re <- fit$lineage_model$lin_interactions[[site]]$has_re
    b_plts <- NULL

    ddf <- as_draws_df(fit$fit()$draws())

    b_df <- do.call(rbind, lapply(
        seq_along(contr_names),
        function(x) 
        {
            data.frame(
                motif=name_transf(contr_names[x]),
                level="Mean",
                cov=cov_names,
                vname=paste0("site_",sname,"_b[",x,",",seq_along(cov_names),"]"),
                re_name=NA
            )
        }
    ))

    if (has_re) 
    {
        re_levels <- fit$lineage_model$lin_interactions[[site]]$re_levels
        re <- fit$lineage_model$lin_interactions[[site]]$re
        b_df_re <- do.call(rbind, lapply(
            seq_along(contr_names),
            function(x) 
            {
                do.call(rbind, lapply(
                    seq_along(re_levels),
                    function(y)
                    {
                        data.frame(
                        motif=name_transf(contr_names[x]),
                        level=paste0(re, " ", re_levels[y]),
                        cov=cov_names,
                        vname=paste0("site_",sname,"_b[",x,",",seq_along(cov_names),"]"),
                        re_name=paste0("site_",sname,"_re_b[",x,",",y,",",seq_along(cov_names),"]")
                    )
                    }
                ))
            }
        ))
        b_df <- rbind(b_df, b_df_re)
    }

    b_data <- do.call(rbind, apply(
        b_df, 1,
        function(x)
        {
            data.frame(motif=x[1], 
                level=x[2], 
                cov=x[3], 
                value=unlist(ddf[, x[4]]) + if(is.na(x[5])) 0 else unlist(ddf[, x[5]]),
                draw=1:nrow(ddf))
        }
    )) %>% 
    filter(!(motif %in% exclude_motifs)) %>% 
    filter(!(level %in% exclude_levels)) %>%
    group_by(motif, level)

    b_plts <- b_data %>%
        split(b_data$cov) %>%
         map(~ .x %>%
            plot_density_ridge + facet_grid(rows=vars(motif)))  

    lapply(seq_along(b_plts), function(x) b_plts[[x]] + ggtitle(names(b_plts)[[x]])) %>%
        wrap_plots(nrow=1) + 
        plot_layout(guides='collect', axes='collect')
}

#Predicted lineage effects
plot_lin_effects <- function(fit)
{
    N_lin <- lnm$lineage_model$n_lin
    rn_cov <- lnm$lineage_model$rn_cov
    N_covar <- length(rn_cov)
    b_varnames <- as.vector(unlist(sapply(1:N_lin, function(x) paste0("lineage_coeffs[",x,",",1:N_covar,"]"))))
    b_new_names <- as.vector(unlist(sapply(1:N_lin, function(x) paste0("lineage_",x,"_b_",rn_cov[1:N_covar]))))
    b_cov_name <- as.vector(unlist(sapply(1:N_lin, function(x) paste0("lineage_",x,"_b_",rn_cov[1:N_covar]))))

    ddf <- as_draws_df(fit$fit()$draws())
    name_lookup <- b_varnames
    names(name_lookup) <- b_new_names
    
    ddf %>% 
        select(all_of(b_varnames)) %>%
        rename(all_of(name_lookup)) %>% 
        plot_density_ridge + 
        facet_grid(rows=vars(name), cols ,scales="free_y")
}

#Plot lineage ne(t) trajectories
plot_traj <- function(fit, with_resi=F, offs=0)
{
    n_grid <- fit$data_list$N_grid
    n_lin <- fit$data_list$N_lin
    grid_ts <- offs - fit$data_list$grid_endpts

    ddf <- as_draws_df(fit$fit()$draws())

    last_grid <- fit$data_list$grid_idx[cumsum(fit$data_list$lin_sizes)]

    traj_names <- NA
    if (with_resi)
    {
        traj_names <- lapply(1:n_lin, function(x) paste0("lin_traj[",x,",",1:n_grid,"]"))
    } else {
        traj_names <- lapply(1:n_lin, function(x) paste0("lin_traj_noresi[",x,",",1:n_grid,"]"))
    }
    traj_df <- do.call(rbind, lapply(1:n_lin, function(i) {
        df <- ddf[, traj_names[[i]]]
        colnames(df) <- 1:n_grid
        df$draw <- 1:nrow(df)
        df$lin <- i
        df %>% pivot_longer(!c(draw, lin), names_to = "x", values_to = "y")
    }))

    traj_df %>% 
        mutate(y=ifelse(as.integer(x) > last_grid[lin], NA, y)) %>%
        group_by(lin, x) %>% 
        summarise(as_tibble_row(quantile(y, c(.025, .5, .975),na.rm=T),
            .name_repair = \(x) paste0('q', parse_number(x)))) %>% 
        mutate(x=grid_ts[as.integer(x)]) %>%
        ggplot(aes(x=x, y=q50)) + 
        geom_line() + 
        geom_ribbon(aes(ymin = q2.5, ymax = q97.5), fill = "grey70",alpha=.5) +
        facet_grid(cols=vars(lin)) + 
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
}

#Plot lineage r(t) trajectories
plot_rt <- function(fit, offs=0, with_resi=F)
{
    n_grid <- fit$data_list$N_grid
    n_lin <- fit$data_list$N_lin
    grid_ts <- offs - fit$data_list$grid_endpts

    ddf <- as_draws_df(fit$fit()$draws())

    last_grid <- fit$data_list$grid_idx[cumsum(lnm$data_list$lin_sizes)]

    traj_names <- NA
    if (!with_resi){
        traj_names <- lapply(1:n_lin, function(x) paste0("lin_rt_noresi[",x,",",1:n_grid,"]"))
    } else {
        traj_names <- lapply(1:n_lin, function(x) paste0("lin_rt_resi[",x,",",1:n_grid,"]"))
    }
    
    traj_df <- do.call(rbind, lapply(1:n_lin, function(i) {
        df <- ddf[, traj_names[[i]]]
        colnames(df) <- 1:n_grid
        df$draw <- 1:nrow(df)
        df$lin <- i
        df %>% pivot_longer(!c(draw, lin), names_to = "x", values_to = "y")
    }))

    traj_df %>% 
        mutate(y=ifelse(as.integer(x) > last_grid[lin], NA, y)) %>%
        group_by(lin, x) %>% 
        summarise(as_tibble_row(quantile(y, c(.025, .5, .975), na.rm =T),
        .name_repair = \(x) paste0('q', parse_number(x)))) %>% 
        mutate(x=grid_ts[as.integer(x)]) %>%
        ggplot(aes(x=x, y=q50)) + 
        geom_line() + 
        geom_ribbon(aes(ymin = q2.5, ymax = q97.5), fill = "grey70",alpha=.5) +
        facet_grid(cols=vars(lin)) + 
        theme_bw() +
        theme(axis.text.x=element_text(angle=60, hjust=1, vjust=1))
}

#Plot global growth rate trend
plot_trend <- function(fit, offs, with_cvs = F)
{
    n_grid <- fit$data_list$N_grid
    n_lin <- fit$data_list$N_lin
    grid_ts <- offs - fit$data_list$grid_endpts

    ddf <- as_draws_df(fit$fit()$draws())
    traj_names <- paste0("global_trend[",1:n_grid,"]")
    
    df <- ddf[, traj_names]


    if(with_cvs) {
        global_eff_names <- paste0("global_effect[",1:n_grid,"]")
        df2 <- ddf[, global_eff_names]
        df <- df + df2
    }

    colnames(df) <- 1:n_grid
    df$draw <- 1:nrow(df)
    
    df %>% 
        pivot_longer(!c(draw), names_to = "x", values_to = "y") %>% 
        group_by(x) %>% 
        summarise(as_tibble_row(quantile(y, c(.025, .5, .975)),
            .name_repair = \(x) paste0('q', parse_number(x)))) %>% 
        mutate(x=grid_ts[as.integer(x)]) %>%
        ggplot(aes(x=x, y=q50)) + 
        geom_line() + 
        geom_ribbon(aes(ymin = q2.5, ymax = q97.5), fill = "grey70",alpha=.5) +
        labs(y="Mean r(t) trend",x="Year") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
            panel.grid.minor.x=element_blank())
}
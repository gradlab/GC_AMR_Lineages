
name_transf <- function(x) str_replace_all(x, "_", " ")

#Compute predicted r(t) effects for use in tables
get_effects <- function(fit, site, select_motifs=c(), exclude_motifs=c(), exclude_levels=c())
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
    covs <- fit$data_list$covars

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
                    cov=x[3])
        }
    )) %>% 
        pivot_longer(!c(draw, motif, level, cov), 
            names_to="grid", 
            names_pattern = "(\\d+)", 
            names_transform = as.integer) %>%
        group_by(draw, motif, level, grid) %>%
        summarise(y=sum(value), .groups = "keep")
}

#Summarise effect averaged across a given time windows
summarise_eff <- function(fit, site, range, relative=F, offs=0, ci_level=c(.025, .5, .975))
{
    grid_ts <- offs - fit$data_list$grid_endpts
    effs <- get_effects(fit, site) %>% 
        mutate(t=grid_ts[as.integer(grid)]) %>%
        ungroup()

    if (relative)
    {
        effs <- effs %>% 
            group_by(motif, level, draw) %>%
            mutate(y=y-y[which.min(grid)]) %>% 
            ungroup()
    }

    effs %>% 
        select(!grid) %>%
        filter((t >= range[1]) & (t < range[2])) %>%
        group_by(draw, motif, level) %>%
        summarise(y=mean(y), .groups="keep") %>%
        ungroup() %>%
        group_by(motif, level) %>%
        summarise(as_tibble_row(quantile(y, ci_level),
            .name_repair = \(x) paste0('q', parse_number(x))),.groups="keep")
}

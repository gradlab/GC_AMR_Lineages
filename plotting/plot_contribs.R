source("models/comp_effects.R")
library(ggpubr) 

#Plot contributions of loci for each lineage
plot_contribs_loci <- function(lnm, offs=2020) 
{
    df <- data.frame(lin=c(), grid=c(), Contribution=c(), y=c())
    
    site_name <- c("italic('gyrA')", "italic('penA')", "italic('ponA')", "italic('mtr')", "'rRNA 23S'")
    
    ddf <- as_draws_df(lnm$fit()$draws())

    n_grid <- lnm$data_list$N_grid
    n_lin <- lnm$data_list$N_lin

    lins <- 1:n_lin
    
    last_grid <- lnm$data_list$grid_idx[cumsum(lnm$data_list$lin_sizes)]
    first_grid <- lnm$data_list$grid_idx[c(1,1+cumsum(lnm$data_list$lin_sizes)[1:(n_lin-1)])]

    grid_ts <- offs - lnm$data_list$grid_endpts
    for (s in 1:5)
    {
        inter <- lnm$lineage_model$lin_interactions[[s]]
        effs <- get_effects(lnm, s)
        effs <- effs %>% group_by(motif, level, grid) %>% summarise(y=mean(y))

        for (lin in lins)
        {

            lin_dets <- name_transf(inter$contr_names[inter$design_matrix[lin,]])
            lin_re <- "Mean"
            if (inter$has_re) 
            {
                lin_re <- name_transf(paste0(inter$re, " ", inter$re_levels)[inter$level_idx[lin]])
            } 
            lin_contr <- effs %>% 
            ungroup() %>%
            filter(motif %in% lin_dets & level %in% lin_re) %>% 
                group_by(grid) %>%
                summarise(y=sum(y)) %>%
                mutate(Contribution = site_name[s], lin=lin)
            df <- rbind(df, lin_contr)
        }
    }
    
    rt_names <- lapply(lins, function(x) paste0("lin_rt_resi[",x,",",1:n_grid,"]"))
    resi_names <- lapply(lins, function(x) paste0("f_residual[",x,",",1:n_grid,"]"))
    resi_a_names <- paste0("lin_rt_intercepts[",lins,"]")
    sd_draws <- ddf[,"residual_sd"]

    for (lin in lins)
    {
        res_lin <- apply(unlist(sd_draws) * as.matrix(ddf[, resi_names[[lin]]]), 2, mean)
        res_df <- data.frame(lin=lin, grid=1:n_grid, Contribution="Residual", y=as.vector(res_lin))
        bgr_df <- data.frame(lin=lin, grid=1:n_grid, Contribution="'Lineage background'", y=mean(unlist(ddf[,resi_a_names[lin]])))

        df <- rbind(df, res_df, bgr_df)
    }

    leg <- c(site_name, "'Lineage background'", "Residual")

    df_effs <- df %>% 
        mutate(y=ifelse((as.integer(grid) > last_grid[lin]) | 
            (as.integer(grid) < first_grid[lin]), NA, y)) %>%
        group_by(lin, grid) %>%
        mutate(ave_rt=sum(y), 
            lin=factor(paste0("Lineage ", lin), levels=paste0("Lineage ", lins), ordered=T),
            Contribution=factor(Contribution, levels=leg, ordered=T)) %>%
        ungroup()

    has_epoch <- function(lin, ep) {
        ep_upper <- c(2007, 2010, 2012, 2019)
        ep_lower <- c(-Inf, 2007, 2010, 2012)

        return((grid_ts[last_grid[lin]] < ep_upper[ep]) & 
            (grid_ts[first_grid[lin]] >= ep_lower[ep]))
    }

    df_ci <- get_contrib_cis(lnm) %>%
        mutate(CI95 = sprintf("'%.2f (%.2f, %.2f)'",round(q50,3), round(q2.5,3), round(q97.5,3))) %>%
        mutate(CI95=ifelse(has_epoch(lin, epoch), CI95, NA)) %>%
        mutate(CI95=ifelse(!is.na(q50), CI95, NA)) 

    
    head(df_ci %>% filter(Contribution == "AMR")) %>% print

    plts <- list()

    int_ts <- abs(grid_ts-as.integer(grid_ts)) < 1e-6

    table_theme <- ttheme(
        tbody.style = tbody_style(
            parse=TRUE
        )
    )

    for (i in lins)
    {
        parse.labels <- function(x) sapply(x, function(y) parse(text = y))
        l <- paste0("Lineage ", i)
        p1 <- df_effs %>% 
            ungroup() %>%
            filter(lin == l) %>%
            mutate(grid = factor(grid), 
                grid = factor(grid, levels = rev(levels(grid)))) %>%
            ggplot(aes(x=grid, y=y, fill=Contribution)) +
                geom_col(width=1, show.legend = T) + 
                geom_line(aes(x=as.integer(grid), y=ave_rt)) +
                geom_hline(yintercept=0, linetype='dashed') +
                theme_bw() + 
                ggtitle(l) +
                labs(x="Year", y="Growth Rate", fill="Average Effect") +
                scale_x_discrete(breaks = c(1:n_grid)[int_ts],labels=grid_ts[int_ts]) +
                scale_fill_manual(limits=leg, values =brewer.pal(name ="Dark2", length(leg)),
                    breaks=leg, labels = parse.labels, drop = FALSE) +
                theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))

        p2 <- df_ci %>% 
            ungroup() %>%
            filter(lin == l) %>%
            select(epoch_s, CI95, Contribution) %>% 
            complete(epoch_s, Contribution) %>%
            mutate(CI95 = ifelse(is.na(CI95), "-", CI95)) %>%
            pivot_wider(names_from = epoch_s, values_from = CI95) %>%
            ggtexttable(rows = NULL, theme=table_theme) 

        p3 <- ggarrange(p1, p2,
            ncol = 1, nrow = 2,
            heights = c(1, 2),
            widths = c(1),
            align="v")
        plts <- c(plts, list(p3))
    }
    return(plts)
}

get_contrib_cis <- function(lnm, offs=2020)
{
    df <- data.frame(lin=c(), grid=c(), Contribution=c(), y=c())
    
    site_name <- c("italic('gyrA')", "italic('penA')", "italic('ponA')", "italic('mtr')","'rRNA 23S'")
    
    ddf <- as_draws_df(lnm$fit()$draws())

    n_grid <- lnm$data_list$N_grid
    n_lin <- lnm$data_list$N_lin

    lins <- 1:n_lin
    
    last_grid <- lnm$data_list$grid_idx[cumsum(lnm$data_list$lin_sizes)]
    grid_ts <- offs - lnm$data_list$grid_endpts

    assign_epochs <- function(x) sapply(x, function(y)
        if(y < 2007) 1L else if (y < 2010) 2L else if (y < 2012) 3L else 4L
    )

    for (s in 1:5)
    {
        inter <- lnm$lineage_model$lin_interactions[[s]]
        effs <- get_effects(lnm, s)

        for (lin in lins)
        {

            lin_dets <- name_transf(inter$contr_names[inter$design_matrix[lin,]])
            lin_re <- "Mean"
            if (inter$has_re) 
            {
                lin_re <- name_transf(paste0(inter$re, " ", inter$re_levels)[inter$level_idx[lin]])
            } 
            lin_contr <- effs %>% 
            ungroup() %>%
            filter(motif %in% lin_dets & level %in% lin_re) %>% 
                group_by(grid, draw) %>%
                summarise(y=sum(y)) %>%
                mutate(Contribution = site_name[s], lin=lin) 
            df <- rbind(df, lin_contr)
        }
    }
    
    rt_names <- lapply(lins, function(x) paste0("lin_rt_resi[",x,",",1:n_grid,"]"))
    resi_names <- lapply(lins, function(x) paste0("f_residual[",x,",",1:n_grid,"]"))
    resi_a_names <- paste0("lin_rt_intercepts[",lins,"]")
    sd_draws <- ddf[,"residual_sd"]

    for (lin in lins)
    {
        res_lin <- unlist(sd_draws) * as.matrix(ddf[, resi_names[[lin]]])
        res_df <- as.data.frame(res_lin)
        colnames(res_df) <- 1:n_grid
        res_df$draw <- 1:nrow(res_df)

        bgr_df <- as.data.frame(unlist(ddf[,resi_a_names[lin]]) %*% t(rep(1,n_grid)))
        colnames(bgr_df) <- 1:n_grid
        bgr_df$draw <- 1:nrow(bgr_df)


        res_df <- res_df %>% 
            pivot_longer(!draw, names_to = "grid", values_to = "y") %>%
            mutate(grid=as.integer(grid), Contribution="Residual", lin=lin)
        bgr_df <- bgr_df %>% 
            pivot_longer(!draw, names_to = "grid", values_to = "y") %>%
            mutate(grid=as.integer(grid), Contribution="'Lineage background'", lin=lin)

        df <- rbind(df, res_df, bgr_df)
    }


    leg <- c(site_name, "'Lineage background'", "Residual", "AMR Total", "Total")

    df %<>% mutate(epoch=assign_epochs(grid_ts[grid])) %>%
        group_by(epoch, lin, Contribution, draw) %>%
        summarise(y=mean(y))
    
    df <- df %>% group_by(epoch, lin, draw) %>%
        group_by(epoch, lin, draw) %>%
        summarise(y=sum(y), Contribution="Total") %>%
        bind_rows(df, .)

    df %>% group_by(epoch, lin, draw) %>%
        summarise(y=ifelse(sum(Contribution %in% site_name)>0,sum(ifelse(Contribution %in% site_name, y, 0)), NA), Contribution="AMR Total") %>%
        bind_rows(df, .) %>%
        group_by(epoch, lin, Contribution) %>%
        summarise(as_tibble_row(quantile(y, c(.025, .5, .975), na.rm=T),
            .name_repair = \(x) paste0('q', parse_number(x))),.groups="keep") %>%
        mutate(lin=factor(paste0("Lineage ", lin), levels=paste0("Lineage ", lins), ordered=T),
            Contribution=factor(Contribution, levels=leg, ordered=T),
            epoch_s=c("1993-2007", "2007-2010", "2010-2012", "2012-2019")[epoch])
}

comp_resi_by_yr <- function(lnm, offs=2020)
{
    df <- data.frame(lin=c(), grid=c(), Contribution=c(), y=c())

    ddf <- as_draws_df(lnm$fit()$draws())

    n_grid <- lnm$data_list$N_grid
    n_lin <- lnm$data_list$N_lin

    lins <- 1:n_lin
    
    last_grid <- lnm$data_list$grid_idx[cumsum(lnm$data_list$lin_sizes)]
    grid_ts <- offs - lnm$data_list$grid_endpts

    assign_epochs <- function(x) 
    {
        u <- 2:2020
        l <- 1:2019
        sapply(x, function(y) which((y < u) & (y >= l)))
    }
    
    rt_names <- lapply(lins, function(x) paste0("lin_rt_resi[",x,",",1:n_grid,"]"))
    resi_names <- lapply(lins, function(x) paste0("f_residual[",x,",",1:n_grid,"]"))
    resi_a_names <- paste0("lin_rt_intercepts[",lins,"]")
    sd_draws <- ddf[,"residual_sd"]

    for (lin in lins)
    {
        res_lin <- unlist(sd_draws) * as.matrix(ddf[, resi_names[[lin]]])
        res_df <- as.data.frame(res_lin)
        colnames(res_df) <- 1:n_grid
        res_df$draw <- 1:nrow(res_df)

        bgr_df <- as.data.frame(unlist(ddf[,resi_a_names[lin]]) %*% t(rep(1,n_grid)))
        colnames(bgr_df) <- 1:n_grid
        bgr_df$draw <- 1:nrow(bgr_df)


        res_df <- res_df %>% 
            pivot_longer(!draw, names_to = "grid", values_to = "y") %>%
            mutate(grid=as.integer(grid), Contribution="Residual", lin=lin)
        bgr_df <- bgr_df %>% 
            pivot_longer(!draw, names_to = "grid", values_to = "y") %>%
            mutate(grid=as.integer(grid), Contribution="'Lineage background'", lin=lin)

        df <- rbind(df, res_df, bgr_df)
    }

    df %<>% group_by(grid, lin, draw) %>%
        summarise(y=(y)) %>% 
        ungroup() %>%
        mutate(year=assign_epochs(grid_ts[grid])) %>%
        group_by(year, lin, draw) %>%
        summarise(y=mean(y)) %>% 
        ungroup() %>%
        mutate(lin=factor(paste0("Lineage ", lin), levels=paste0("Lineage ", lins), ordered=T))

    df
}
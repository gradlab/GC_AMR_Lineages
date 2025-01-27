source("models/comp_effects.R")

#Generate a dataframe for the results summary table
make_ci_dfs <- function(lnm)
{
    cip_use <- c(1993, 2007)
    cip_nouse <- c(2007, 2020)
    
    azi_use <- c(1993, 2012)
    azi_nouse <- c(2012, 2020)

    rest_use <- c(1993, 2012)
    rest_nouse <- c(2012, 2020)

    site_data <- data.frame(
        site_id=1:5,
        site_name=c("gyrA", "penA", "ponA", "mtrCDE", "rRNA 23S"),
        use_lo=c(cip_use[1], rest_use[1], rest_use[1], azi_use[1], azi_use[1]),
        use_hi=c(cip_use[2], rest_use[2], rest_use[2], azi_use[2], azi_use[2]),
        nouse_lo=c(cip_nouse[1],rest_nouse[1],rest_nouse[1], azi_nouse[1], azi_nouse[1]),
        nouse_hi=c(cip_nouse[2],rest_nouse[2],rest_nouse[2], azi_nouse[2], azi_nouse[2])
    )

    f_site_abs <- function(x)
    {
        tmp_use <- summarise_eff(lnm, x[[1]], range=c(x[[3]],x[[4]]), offs=2020, relative=F)
        tmp_use$use <- T

        tmp_nouse <- summarise_eff(lnm, x[[1]], range=c(x[[5]],x[[6]]), offs=2020, relative=F)
        tmp_nouse$use <- F

        df <- rbind(tmp_use, tmp_nouse)
        df$site <- x[[2]]
        df
    }

    f_site_rel <- function(x)
    {
        df <- summarise_eff(lnm, x[[1]], range=c(x[[3]],x[[4]]), offs=2020, relative=T)
        df$site <- x[[2]]
        df
    }

    abs_df <- do.call(rbind, 
        lapply(1:nrow(site_data),
            function(y) f_site_abs(site_data[y,])
        )
    )

    rel_df <- do.call(rbind, 
        lapply(1:nrow(site_data),
            function(y) f_site_rel(site_data[y,])
        )
    )

    list(abs_df=abs_df, rel_df=rel_df)
}

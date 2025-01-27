library(ggplot2)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
#Plot trend 
pal <- setNames(c("#E5C494", 
        "#FC8D62",
        "#E78AC3", 
        "#8DA0CB", 
        "#A6D854", 
        "#FFD92F", 
        "#66C2A5"),
    c("PEN", "FQ", "CRO_250", "CEPH_LO", "CRO_250_CO", "CEPH_LO_CO", "Other"))

scale_fill_abx <- function(...){
    ggplot2:::manual_scale(
        'fill', 
        values = pal, 
        ...
    )
}

scale_color_abx <- function(...){
    ggplot2:::manual_scale(
        'color', 
        values = pal,
        ...
    )
}

plot_usage <- function(usage_df)
{
    drugs_to_keep <- c("Cefixime 400 mg", 
        "Ceftriaxone 125 mg",
        "Ceftriaxone 250 mg",
        "Other Cephalosporins",
        "Ciprofloxacin",
        "Ofloxacin",
        "Penicillin")
    FQ <- c("Ciprofloxacin",
        "Ofloxacin")
    CEPH_LO <- c("Cefixime 400 mg",
        "Other Cephalosporins",
        "Ceftriaxone 125 mg")

    plotting_ord <- c("Other", "FQ", "CEPH_LO", "CRO_250", "CRO_250_CO", "CEPH_LO_CO", "PEN")

    lab_df <- data.frame(x=c(1989, 1993, 2007, 2010, 2012),
        y=c(100, 100, 100, 100, 100),
        lab=c("Penicillin no longer recommended",
            "Fluoroquinolones recommended",
            "Fluoroquinolones no longer recommended",
            "Azithromycin co-treatment introduced",
            "Only Ceftriaxone 250mg recommended")
    )

    usage_df %>% pivot_longer(cols = !Year, names_to = "Drug") %>% 
        mutate(Drug = factor(Drug, levels = plotting_ord, ordered=T)) %>%
        arrange((Drug)) %>%
        group_by(Year) %>%
        mutate(ymax=cumsum(value)) %>%
        mutate(ymin=ymax-value, xmin = Year, xmax=Year+1L) %>%
        ggplot(aes(fill=Drug, color=Drug)) +
            geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),linewidth=0.01) +
            geom_vline(data=lab_df, aes(xintercept=x),linetype="dashed",color="gray10", alpha=.8) + 
            geom_text(data=lab_df, aes(x=x, y=y, label = lab), size=rel(5), vjust = -1, color="gray10", hjust = 1, nudge_y=-1, angle=90,inherit.aes = FALSE) + 
            scale_fill_abx(na.value="gray65", 
                breaks = c("PEN", "FQ", "CRO_250", "CEPH_LO", "CRO_250_CO", "CEPH_LO_CO", "Other"),
                labels=c("Penicillin", "Fluoroquinolones", "Ceftriaxone 250mg", "Other Cephalosporins", "Azithromycin + Ceftriaxone 250mg", "Azithromycin + Other Cephalosporins", "Other"))+
            scale_color_abx(
                na.value="gray65",
                breaks = c("PEN", "FQ", "CRO_250", "CEPH_LO", "CRO_250_CO", "CEPH_LO_CO", "Other"),
                labels=c("Penicillin", "Fluoroquinolones", "Ceftriaxone 250mg", "Other Cephalosporins", "Azithromycin + Ceftriaxone 250mg", "Azithromycin + Other Cephalosporins", "Other"))+
            scale_x_continuous(breaks=seq(from=1988, to=2020, by=5))+
            xlim(c(1988,2020))+
            ylim(c(0,100.1))+
            labs(y="% of Primary Treatment",color="Treatment", fill="Treatment") + 
            theme_minimal() +
            theme(aspect.ratio=.6,
                axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_text(size=rel(1.8)),
                axis.title.y = element_text(size=rel(2.2)),
                legend.text = element_text(size=rel(1.6)),
                legend.title = element_text(size=rel(1.8)),
                legend.justification = "left")
}

plot_trends <- function(preva_df)
{
    plotting_ord <- c("Penicillin", "Ciprofloxacin", "Ceftriaxone", "Cefixime", "Azithromycin")
    abbrev <-  c("PEN", "FQ", "CRO_250", "CEPH_LO", "CEPH_LO_CO")
    names(abbrev) <- plotting_ord
    yr_min <- min(preva_df$Year)

    preva_df %>%
        filter(Drug %in% plotting_ord) %>%
        mutate(Drug = abbrev[Drug]) %>% 
        mutate(Drug = factor(Drug, levels = abbrev, ordered=T)) %>%
        ggplot(aes(x=Year, y=Percent, color=Drug)) +
            geom_step(lwd=1.5)+
            geom_vline(aes(xintercept=1993), linetype="dashed", color="gray25", alpha=0.8) +
            geom_vline(aes(xintercept=1989), linetype="dashed", color="gray25", alpha=0.8) +
            geom_vline(aes(xintercept=2012), linetype="dashed", color="gray25", alpha=0.8) +
            geom_vline(aes(xintercept=2010), linetype="dashed", color="gray25", alpha=0.8) +
            geom_vline(aes(xintercept=2007), linetype="dashed", color="gray25", alpha=0.8) +
            annotate("rect", xmin = 1988, xmax = yr_min, ymin = 0, ymax = 50, 
                alpha = .3, fill="gray25") +
            scale_color_abx(na.value="gray65",
            breaks = abbrev,
                labels=plotting_ord)+
            scale_x_continuous(breaks=seq(from=1988, to=2020, by=5))+
            xlim(c(1988,2020))+
            ylim(c(0,40))+
            labs(y="Isolates with Elevated MIC (%)") + 
            theme_minimal() +
            theme(aspect.ratio=.6,
                axis.text.x = element_text(size=rel(1.8)),
                axis.text.y = element_text(size=rel(1.8)),
                axis.title.y = element_text(size=rel(2.2)),
                axis.title.x = element_text(size=rel(2.2)),
                legend.text=element_text(size=rel(1.6)),
                legend.title = element_text(size=rel(1.8)),
                legend.justification = "left")
}
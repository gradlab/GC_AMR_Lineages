library(ape)
library(tidyverse)
library(cmdstanr)
library(ggplot2)
library(stringr)
library(stringi)
library(glue)
library(posterior)
library(ggplot2)
library(bayesplot)
library(gt)
library(gridExtra)

set.seed(12345)

source("models/model.R")

PTS_PER_YR = 4L 
dat_in <- read.csv("data/metadata/usage.csv")[1:3]

preva_df <- read.csv("data/metadata/NG - Prevalence of Resistance or Elevated MICs to Specific Antimicrobials by Year (GISP 2000-2022).csv")

t_begin <- 1993

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

#Data for regression model
usage_data <- dat_in %>%
    filter(Year %in% t_begin:2020) %>%
    filter(Treatment %in% drugs_to_keep) %>%
    filter(Treatment != "Penicillin") %>%
    pivot_wider(names_from=Treatment, values_from=Percent) %>%
    mutate(FQ = rowSums(across(all_of(FQ))),
        CEPH_LO = rowSums(across(all_of(CEPH_LO))),
        CRO_250 = .data[["Ceftriaxone 250 mg"]], .keep="unused") %>%
    mutate(Other = 100-rowSums(across(!all_of("Year")))) %>%
    mutate(AZI_CO = ifelse(Year < 2010, 0, CEPH_LO + CRO_250))

#Data for usage plot 
usage_data4plot <- dat_in %>%
    filter(Year %in% 1988:2020) %>%
    filter(Treatment %in% drugs_to_keep) %>%
    pivot_wider(names_from=Treatment, values_from=Percent) %>%
    mutate(FQ = rowSums(across(all_of(FQ))),
        CEPH_LO = rowSums(across(all_of(CEPH_LO))),
        CRO_250 = .data[["Ceftriaxone 250 mg"]],
        PEN = Penicillin, .keep="unused") %>%
    mutate(Other = 100-rowSums(across(!all_of("Year")))) %>%
    mutate(CEPH_LO_CO = ifelse(Year < 2010, 0, CEPH_LO),
    CRO_250_CO = ifelse(Year < 2010, 0, CRO_250),
    CEPH_LO = ifelse(Year >= 2010, 0, CEPH_LO),
    CRO_250 = ifelse(Year >= 2010, 0, CRO_250))

source("plotting/plot_usage_trend.R")
p1 <- plot_usage(usage_data4plot)
p2 <- plot_trends(preva_df)
pdf("plots/usage.pdf",13, 13)
plot(p1/p2)
dev.off()

#Remove baseline and unused covariate
usage_data <- usage_data %>%
    select(!any_of(c("CRO_250", "Other"))) %>% 
    mutate(Year = max(Year)-Year) %>%
    mutate(across(!Year, function(x) x / 100))

#Load lineages from preprocessing step
lin_data <- readRDS("data/lineages.rds")

n_lin <- length(lin_data$lin_data)
lin_trees <- lapply(lin_data$lin_data, function (x) x$phylo)
lin_alleles <- do.call(rbind,lapply(lin_data$lin_data, function (x) x$alleles))
mrsts <- sapply(lin_data$lin_data, function (x) x$mrst)

wt_types <- c("gyrA_91S/95D", 
    "pbp_421L", 
    "mtr_promoter_WT", 
    "mtr_r_GC_allele", 
    "mtr_c_GC_allele", 
    "mtr_d_GC_allele",
    "penA_A501_A",
    "penA_G543_G",
    "penA_A311_A", 
    "penA_Mosaic_No",
    "rRNA_23s_C2611T_C"
)
#Fix names
rename_lookup <- c(mtr_promoter_mosaic="mtr_promoter_divergent_promoter",
    mtrR_LOF="mtr_r_LOF",
    penA_mosaic="penA_Mosaic_Yes",
    mtr_promoter_Adel="mtr_promoter_Adel_13bp_invertedRepeat",
    ponA_421P="pbp_421P"
)
#Prepare lineage covariates
lineage_data <- as.data.frame(lin_alleles) %>% 
    select(!c(parC, penA_N513, penA_A517)) %>%
    mutate(across(everything(), ~ as.data.frame.matrix(table(row_number(), .x) * NA^(is.na(.x)) > 0))) %>% 
    unpack(where(is.data.frame), names_sep = "_") %>% 
    select(!all_of(wt_types)) %>%
    select(where(~ (sum(.x) > 0) && (sum(!.x) > 0))) %>%
    rename_with(~gsub("[A-Z](\\d+)_", "\\1", .x, perl=TRUE), matches("penA_[A-Z]\\d+.+")) %>%
    rename_with(~gsub("[A-Z](\\d+)[A-Z]_", "\\1", .x, perl=TRUE), matches("rRNA_23s_[A-Z]\\d+.+")) %>%
    rename_with(~gsub("mtr_([cd])", "mtr\\U\\1", .x, perl=TRUE), matches("mtr_[cd].*")) %>%
    rename(all_of(rename_lookup)) %>%
    as.data.frame()

lineage_data$mrst <- mrsts
#ParC baseline
parC_wt <- "86D/87S/91E"
lineage_data$parC <- factor(lin_alleles$parC, exclude = parC_wt)

gene_names <- colnames(lineage_data)
#Build Model
lnm <- lineage_model(lin_trees, lineage_data, usage_data, "Year", N_grid_per_unit=PTS_PER_YR) %>%
    add_site_interaction(gene_names[str_detect(gene_names, "(gyrA)")], "FQ", re="parC", site_name = "gyrA") %>%
    add_site_interaction(gene_names[str_detect(gene_names, "(penA)")], c("CEPH_LO"), site_name = "penA") %>%
    add_site_interaction(gene_names[str_detect(gene_names, "(ponA)")], c("CEPH_LO"), site_name = "pbp") %>%
    add_site_interaction(gene_names[str_detect(gene_names, "(mtr)")], c("FQ", "CEPH_LO", "AZI_CO"), sum_covs = T, site_name = "mtr") %>%
    add_site_interaction(gene_names[str_detect(gene_names, "(rRNA)")], c("AZI_CO"), site_name = "rRNA23S") %>%
    generate_model 
#Fit model
fitted <- lnm$sample(iter_warmup=1000,#2000, 
    iter_sampling=1000, 
    chains=4, 
    parallel_chains=4, 
    max_treedepth=14, 
    adapt_delta=0.99)
#Plot
source("models/plot_func.R")

lab_df_mtr <- data.frame(x=c(2010, 2012),
    lab=factor(c("Azithromycin Co-treatment recommended",
        "Only Ceftriaxone 250mg + Azithromycin recommended"), 
        levels=c("Azithromycin Co-treatment recommended",
        "Only Ceftriaxone 250mg + Azithromycin recommended"), ordered=T),
    col=factor(c("CEPH_LO_CO", "CRO_250_CO"), levels=c("CEPH_LO_CO", "CRO_250_CO"), ordered=T)
    )
lab_df_ceph <- data.frame(x=c(2012),
    lab=factor(c("Only Ceftriaxone 250mg + Azithromycin recommended"), 
        levels=c("Only Ceftriaxone 250mg + Azithromycin recommended"),
        ordered=T),
    col=factor(c("CRO_250_CO")))

lab_df_fq <- data.frame(x=c(2007),
    lab=factor(
        "Fluoroquinolones no longer recommended",
        levels=c(
        "Fluoroquinolones no longer recommended")),
    col=factor(c("FQ")))
xlimset <- scale_x_continuous(limits=c(1993, 2013),
        breaks=c(2013, 2008, 2003, 1998, 1993))
pdf("plots/effects_gyrA_small.pdf",7,4)
p1 <- plot_site_effects(lnm, 
    1, 
    offs=2020,
    relative=F,
    exclude_motifs=c("gyrA 91F/95N"), 
    exclude_levels=c("parC 86D/87I/91E", "parC 86D/87S/91G", "parC 86D/87N/91E", "parC 86N/87S/91E")) +
     #   "parC 86D/87S/91G")) + 
    xlimset +
    geom_vline(data=lab_df_fq, aes(xintercept=x, color=col), alpha=0.9, width=2)+
    scale_color_abx(breaks=lab_df_fq$col, labels=lab_df_fq$lab)+
    ylab(str_wrap("Predicted Growth Rate Effect Relative to Baseline Type", width=40))
plot(p1 & theme(legend.position = 'bottom', legend.title=element_blank()))
dev.off()

pdf("plots/diff_gyrA.pdf",7,7)
plot(plot_site_diff(lnm, 
    1, 
    offs=2020,
    motifs_to_compare=c("gyrA 91F/95A", "gyrA 91F/95G"), 
    exclude_levels=c("parC 86D/87I/91E", "parC 86D/87S/91G", "parC 86D/87N/91E", "parC 86N/87S/91E"))+
     #   "parC 86D/87S/91G")) + 
    xlimset +
    geom_vline(data=lab_df_fq, aes(xintercept=x, color=col), alpha=0.9, width=2)+
    scale_color_abx(breaks=lab_df_fq$col, labels=lab_df_fq$lab)+
    ggtitle("Difference Between Growth Rate Effect of gyrA 91F/95A and gyrA 91F/95G") & 
        theme(legend.position = 'bottom', legend.title=element_blank())
) 
dev.off()

pdf("plots/effects_penA_select.pdf",9,5)
ylimset <- scale_y_continuous(limits=c(-0.5, 1.25), breaks=seq(from=-.5, to=1.25, by=.25))

p_usg <- usage_data4plot %>% 
    select(c(CEPH_LO, CEPH_LO_CO, Year)) %>%
    mutate(CEPH_LO = CEPH_LO+CEPH_LO_CO, .keep="unused") %>%
    rename(value=CEPH_LO)%>%
    mutate(col=factor("CEPH_LO")) %>%
    mutate(ymax=value, ymin=0, xmin = Year, xmax=Year+1L) %>%
    ggplot(aes(fill=col)) +
    geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),linewidth=0.01) +
    scale_fill_abx(breaks="CEPH_LO", labels="Other Cephalosporins",guide="none")+
    geom_vline(data=lab_df_ceph, aes(xintercept=x, color=col), alpha=0.9, width=2)+
    scale_color_abx(breaks=lab_df_ceph$col, labels=lab_df_ceph$lab,guide="none")+
    xlimset +
    ylim(c(0,100.01))+
    labs(y=str_wrap("Other Cephalosporins (%)",15), fill="Treatment") + 
    theme_minimal() +
            theme(legend.position = "none")

p1 <- plot_site_effects(lnm, 
    2, 
    offs=2020,
    relative=F,
    select_motifs=c("penA mosaic"))+xlimset+ylimset+
    geom_vline(data=lab_df_ceph, aes(xintercept=x, color=col), alpha=0.9, width=2)+
    scale_color_abx(breaks=lab_df_ceph$col, labels=lab_df_ceph$lab)+
    ylab(str_wrap("Predicted Growth Rate Effect Relative to Baseline Type", width=30)) & theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank())
plt <-  p1 & theme(legend.position = 'bottom',legend.title=element_blank())#, widths=c(4,2), heights = c(5,5))
plot(plt / p_usg + plot_layout(heights = c(3, 1),guides = "collect")& theme(legend.position = 'bottom'))
dev.off()

## Plot combined MTR effs
ylimset <- scale_y_continuous(limits=c(-1,2.5), breaks=seq(from=-1, to=2.5, by=0.5))

pdf("plots/mtr_marginals.pdf",7,4)
p_mtr1<- plot_site_effects(lnm, 4, offs=2020, relative=F, select_motifs=c("mtrD mosaic", "mtr promoter mosaic")) +
    geom_vline(data=lab_df_mtr, aes(xintercept=x, color=col), alpha=0.9, width=2)+
    scale_color_abx(breaks=lab_df_mtr$col, labels=lab_df_mtr$lab)+
    xlimset + ylimset +
    ylab(str_wrap("Predicted Growth Rate Effect Relative to Baseline Type", width=40))
plt <- p_mtr1 & theme(legend.position = 'bottom', legend.title = element_blank())
plot(plt)
dev.off()


ylimset <- scale_y_continuous(limits=c(-1,3.5), breaks=seq(from=-0.5, to=3.5, by=0.5))

pdf("plots/mtr_combined.pdf",7,7)
p_mtr1c <- plot_site_effects(lnm, 4, offs=2020, relative=F, combine_motifs = c("mtrC mosaic", "mtrD mosaic", "mtr promoter mosaic")) +
    xlimset + ylimset +
    geom_vline(data=lab_df_mtr, aes(xintercept=x, color=col), alpha=0.9, width=2)+
    scale_color_abx(breaks=lab_df_mtr$col, labels=lab_df_mtr$lab)+
    ylab("Predicted Growth Rate Effect Relative to Baseline Type") & 
    theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_mtr2c <- plot_site_effects(lnm, 4, offs=2020, relative=F, combine_motifs = c("mtrD mosaic", "mtr promoter mosaic")) +
    geom_vline(data=lab_df_mtr, aes(xintercept=x, color=col), alpha=0.9, width=2)+
    scale_color_abx(breaks=lab_df_mtr$col, labels=lab_df_mtr$lab)+
    xlimset + ylimset +
    ylab("Predicted Growth Rate Effect Relative to Baseline Type") & 
    theme(legend.title = element_blank())

plt <- p_mtr1c + p_mtr2c +
    plot_layout(nrow=2, guides = 'collect',axis_title="collect") & theme(legend.position = 'bottom')
plot(plt)
dev.off()

###
#SUPPLEMENTARY FIGURES
###

pdf("plots/effects_gyrA.pdf",10,10)
plot(plot_site_effects(lnm, 
    1, 
    offs=2020,
    relative=F) +
    scale_x_continuous(limits=c(1993, 2013),
        breaks=c(2013, 2008, 2003, 1998, 1993)
    )  +
    ylab("Predicted Growth Rate Effect Relative to Baseline Type") +
    ggtitle("Absolute Growth Rate Effect Relative to Baseline Type")
)
dev.off()

pdf("plots/effects_gyrA_relative.pdf",10,10)
plot(plot_site_effects(lnm, 
    1, 
    offs=2020,
    relative=T) +
    scale_x_continuous(limits=c(1993, 2013),
        breaks=c(2013, 2008, 2003, 1998, 1993)

    ) +
    ylab("Predicted Change in Growth Rate Effect Relative to 2019") +
    ggtitle("Change in Growth Rate Effect Relative to 2019")
)
dev.off()

pdf("plots/effects_penA.pdf",10,10)
plot(plot_site_effects(lnm, 2, offs=2020, relative=F) +
    scale_x_continuous(limits=c(1993, 2013),
        breaks=c(2013, 2008, 2003, 1998, 1993)

    ) +
    ylab("Predicted Growth Rate Effect Relative to Baseline Type") +
    ggtitle("Absolute Growth Rate Effect Relative to Baseline Type")
)
dev.off()

pdf("plots/effects_penA_relative.pdf",10,10)
plot(plot_site_effects(lnm, 2, offs=2020, relative=T) +
    scale_x_continuous(limits=c(1993, 2013),
        breaks=c(2013, 2008, 2003, 1998, 1993)

    ) +
    ylab("Predicted Change in Growth Rate Effect Relative to 2019") +
    ggtitle("Change in Growth Rate Effect Relative to 2019")
)
dev.off()

pdf("plots/effects_ponA.pdf",7,4)
plot(plot_site_effects(lnm, 3, offs=2020, relative=F) +
    scale_x_continuous(limits=c(1993, 2013),
        breaks=c(2013, 2008, 2003, 1998, 1993)

    ) +
    xlimset +
    ylab("Predicted Growth Rate Effect Relative to Baseline Type") +
    ggtitle("Absolute Growth Rate Effect Relative to Baseline Type")
)
dev.off()

pdf("plots/effects_ponA_relative.pdf",10,10)
plot(plot_site_effects(lnm, 3, offs=2020, relative=T) +
    scale_x_continuous(limits=c(1993, 2013),
        breaks=c(2013, 2008, 2003, 1998, 1993)

    ) +
    ylab("Predicted Change in Growth Rate Effect Relative to 2019") +
    ggtitle("Change in Growth Rate Effect Relative to 2019")
)
dev.off()

pdf("plots/effects_mtr.pdf",10,10)
plot(plot_site_effects(lnm, 4, offs=2020,relative=F) +
    scale_x_continuous(limits=c(1993, 2013),
        breaks=c(2013, 2008, 2003, 1998, 1993)

    ) +
    ylab("Predicted Growth Rate Effect Relative to Baseline Type") +
    ggtitle("Absolute Growth Rate Effect Relative to Baseline Type")
)
dev.off()

pdf("plots/effects_mtr_relative.pdf",10,10)
plot(plot_site_effects(lnm, 4, offs=2020,relative=T) +
    scale_x_continuous(limits=c(1993, 2013),
        breaks=c(2013, 2008, 2003, 1998, 1993)
    ) +
    ylab("Predicted Change in Growth Rate Effect Relative to 2019") +
    ggtitle("Change in Growth Rate Effect Relative to 2019")
)
dev.off()

pdf("plots/effects_23s.pdf",10,10)
plot(plot_site_effects(lnm, 5, offs=2020,relative=F) +
    scale_x_continuous(limits=c(1993, 2013),
        breaks=c(2013, 2008, 2003, 1998, 1993)

    ) +
    ylab("Predicted Growth Rate Effect Relative to Baseline Type") +
    ggtitle("Absolute Growth Rate Effect Relative to Baseline Type")
)
dev.off()

pdf("plots/effects_23s_relative.pdf",10,10)
plot(plot_site_effects(lnm, 5, offs=2020,relative=T) +
    scale_x_continuous(limits=c(1993, 2013),
        breaks=c(2013, 2008, 2003, 1998, 1993)
    ) +
    ylab("Predicted Change in Growth Rate Effect Relative to 2019") +
    ggtitle("Change in Growth Rate Effect Relative to 2019")
)
dev.off()

pdf("plots/traj.pdf", 16,16)
plot(plot_traj(lnm,T,2020) + facet_wrap(~lin, ncol=4))
dev.off()

pdf("plots/rt.pdf", 16,16)
plot(plot_rt(lnm,2020,F) + facet_wrap(~lin, ncol=4))
dev.off()

pdf("plots/trend.pdf", 10,10)
plot(plot_trend(lnm,2020,T))
dev.off()

source("plotting/plot_contribs.R")
plts <- plot_contribs_loci(lnm)
ggsave("plots/rt_exp.pdf",  
    marrangeGrob(grobs = plts, nrow=1, ncol=1, top=NULL), width=12, height=7)

df_resi <- comp_resi_by_yr(lnm)
counts_greater <- df_resi %>% 
    group_by(lin,year) %>% 
    summarise(y=mean(y)) %>% 
    group_by(lin) %>% 
    summarise(n=sum(abs(y) > .1)) %>%
    ungroup()

source("plotting/ci_tbl.R")
ci_data <- make_ci_dfs(lnm)

abs_df <- ci_data$abs_df %>%
    mutate(use = ifelse(use, "Predicted Average Effect In Use", "Predicted Average Effect Out of Use"))

rel_df <- ci_data$rel_df %>%
    mutate(use="Average Change in Predicted Effect Relative to 2019") 

cns <- c("Average Change in Predicted Effect Relative to 2019",
    "Predicted Average Effect In Use", 
    "Predicted Average Effect Out of Use")

abs_df %>%
    rbind(rel_df) %>% 
    mutate(CI95 = sprintf("%.2f (%.2f, %.2f)",round(q50,3), round(q2.5,3), round(q97.5,3))) %>%
    mutate(d_level = (sign(q2.5) == sign(q97.5))) %>%
    mutate(d_level = d_level + d_level * (q2.5 > 0.1 | q97.5 < -0.1)) %>%
    mutate(d_level = ifelse(d_level == 0, NA , d_level)) %>%
    select(c(motif, level, CI95, use, site, d_level)) %>%
    pivot_wider(names_from = c(use), values_from = c(CI95, d_level)) %>% 
    ungroup() %>%
    gt(rowname_col = "motif",groupname_col = "site", row_group_as_column=T) %>%
    data_color(
        columns = paste0("d_level_",cns),
        target_columns = paste0("CI95_",cns),
        fn = scales::col_numeric(
            palette = "viridis",
            domain = c(1, 2)
    )
    ) %>% 
    cols_hide(paste0("d_level_",cns)) %>%
    cols_align("center") %>%
    cols_label_with(
        fn = function(x) ifelse(str_detect(x, "CI95_.+"), str_extract(x, "(?<=CI95_).+"), x)
    ) %>%
    tab_options(
        table.font.size = 7,
        quarto.use_bootstrap = TRUE,
        data_row.padding = px(1), 
        column_labels.padding = px(1), 
        heading.padding = px(1)) %>%
    gtsave("tables/res_sum.tex")


r_g_1 <- counts_greater[counts_greater$n >= 1, "lin"] %>% unlist %>% as.character
r_g_3 <- counts_greater[counts_greater$n >= 3, "lin"] %>% unlist %>% as.character

print(paste("Lineages with at least year where mean residual+background exceeded 0.1:", paste(r_g_1, collapse=", ")))
print(paste("Lineages with at 3 years where mean residual+background exceeded 0.1:", paste(r_g_3, collapse=", ")))
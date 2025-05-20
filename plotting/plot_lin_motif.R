library(RColorBrewer)
library(ggtree)
library(treeio)
library(tidyverse)
library(ggplot2)
library(stringr)
library(gt)
library(viridis)
library(gtExtras)

#Generate determinant table
make_motif_df <- function(lin_data, pena_desc)
{
    lin_alleles <- as.data.frame(do.call(rbind, lapply(lin_data, function (x) x$alleles)))
    rownames(lin_alleles) <- NULL
    colnames(lin_alleles) <- c("gyrA_allele", 
        "parC_allele", 
        "pbp_allele",
        "mtrC_allele",
        "mtrD_allele",
        "mtr_promoter",
        "mtrR_allele", 
        colnames(lin_alleles)[8:ncol(lin_alleles)])

    penA_sites <- unique(unlist(apply(penA_desc, 1, function(x) str_extract_all(x[2], "[A-Z]\\d+"))))
    penA_pos <- unlist(sapply(penA_sites, function(x) str_extract_all(x, "\\d+")))
    penA_WT <- unlist(sapply(penA_sites, function(x) str_extract_all(x, "[A-Z]")))

    parC_sites <- c("D86", "S87", "E91")
    gyrA_sites <- c("S91", "D95")
    ponA_sites <- c("L421")

    rrna_sites <- c("C2611")

    wt_for_df <- c("91S/95D", 
        "86D/87S/91E", 
        "421L", 
        "GC_allele",
        "GC_allele",
        "WT", 
        "GC_allele")
    wt_for_df <- c(wt_for_df, penA_WT)

    lin_alleles <- lin_alleles %>% 
        select(!c(penA_A517,penA_N513,penA_A311,penA_A513)) %>%
        separate_wider_delim(gyrA_allele, delim="/", names=paste0("gyrA_", gyrA_sites), cols_remove = F) %>%
        separate_wider_delim(parC_allele, delim="/", names=paste0("parC_", parC_sites), cols_remove = F) %>%
        mutate(across(c(gyrA_allele, parC_allele), ~ factor(ifelse(. %in% wt_for_df, NA, as.integer(factor(.)))))) %>%
        rename(ponA_L421=pbp_allele) %>%
        mutate(across(matches(paste0(c("gyrA", "parC", "ponA", "rRNA_23s"), "_[A-Z]\\d+")), ~ str_extract(., "[A-Z]"))) %>%
        mutate(across(starts_with("mtr"), ~ ifelse(. == "GC_allele" | . == "WT" | . =="GC", "non-m", .))) %>%
        mutate(mtr_promoter = ifelse(mtr_promoter =="divergent_promoter", "mosaic", mtr_promoter)) %>%
        mutate(mtr_promoter = ifelse(mtr_promoter =="Adel_13bp_invertedRepeat", "A-del", mtr_promoter)) %>%
        mutate(penA_Mosaic = ifelse(penA_Mosaic =="Yes", "mosaic", "non-m")) %>%
        mutate(across(matches("penA_[A-Z]\\d+"), ~ ifelse(penA_Mosaic == "mosaic", NA, .))) %>%
        rename(penA_Type = penA_Mosaic) %>%
        mutate(Lineage=factor(row_number())) %>%
        relocate(starts_with("penA"), .before= starts_with("mtr"))

        tbl <- lin_alleles %>% gt(rowname_col = "Lineage") %>%
            sub_missing()

        cn <- colnames(lin_alleles)
        
        for(locus in c("gyrA", "parC"))
        {
    
            tbl <- tbl %>% data_color(
                columns = paste0(locus,"_allele"),
                target_columns = cn[str_detect(cn, paste0(locus, "_[A-Z]\\d+"))],
                palette = "Paired"
            ) %>% 
            cols_hide(paste0(locus,"_allele")) %>%
            gt_add_divider(cn[max(which(str_detect(cn, paste0(locus, "_[A-Z]\\d+"))))],
                sides=c("right"),
                color = "gray65",
                weight = px(6))
        }

        for(locus in c("ponA", "penA", "mtr", "rRNA_23s"))
        {
            wt <- c("non-m", "GC")

            pfun <- scales::col_factor(
                palette = "Paired",
                domain = factor(1:10)
            )

            ## need to take care with mtr locus as it comes last in the table and needs renaming
            if (locus != "rRNA_23s")
            {
                tbl <- tbl %>% 
                    tab_spanner(label = md(paste0("*",locus,"*")), columns = starts_with(locus)) %>%
                    gt_add_divider(cn[max(which(str_detect(cn, paste0(locus, ".+"))))],
                        color = "gray65",
                        weight = px(6))
            } else 
            {
                tbl <- tbl %>% 
                    tab_spanner(label = md("*23S rRNA*"), columns = starts_with(locus)) %>%
                    gt_add_divider(cn[max(which(str_detect(cn, paste0(locus, ".+"))))],
                        color = "gray65",
                        weight = px(6))
            }

            for (col in cn[str_detect(cn, paste0(locus, ".+"))])
            {
                wt <- c(wt, str_extract(str_extract(col, "_[A-Z]\\d+"),"[A-Z]"))

                tbl <- tbl %>%
                    data_color(
                    columns = col,
                    #target_columns = cn[str_detect(cn, paste0(locus, "_[a-zA-Z\\d]+"))],
                    fn = function(x) pfun(factor(as.integer(factor(ifelse(x %in% wt, NA, x)))))
                    )
            }
        }

        tbl <- tbl %>%        
            tab_stubhead(label = "Lineage") %>%
            tab_spanner(label = md("*gyrA*"), columns = starts_with("gyrA")) %>%
            tab_spanner(label = md("*parC*"), columns = starts_with("parC")) %>%
            cols_label_with(
                fn = function(x) ifelse(str_detect(x, "[a-zA-Z]+_[A-Z]\\d+"), str_extract(x, "[A-Z]\\d+"), x)
            ) %>%
            cols_label("penA_Type" = "Type") %>%
            cols_label_with(fn=function(x) str_replace(str_replace(x, "_allele", ""), "_", " ")) %>%
            cols_label_with(fn=function(x) md(ifelse(str_detect(x, ".*mtr.*"), str_replace(x, "(mtr[A-Z]?)", "*\\1*"), x))) %>%
            cols_align("center") %>%
            #cols_width(everything() ~ px(30)) %>%
            tab_options(
                table.font.size = 7,
                quarto.use_bootstrap = TRUE,
                data_row.padding = px(1), 
                column_labels.padding = px(1), 
                heading.padding = px(1)) 
            
    tbl
}

plot_lin_tree <- function(tree, lin_df) {

    to <- ggtree(tree,layout="rectangular")$data %>%
        mutate(lin=lin_df[label,"lin"]) %>%
        filter(!is.na(lin_df[label,"cllab"]))%>%
        arrange(lin)

    to$yord <- rank(to$y)

    lin_df <- lin_df %>% 
        mutate(col_lab = ifelse(div, NA, as.integer(to$yord[lin]) %% 3L),
            col = ifelse(!is.na(cllab), NA, col_lab), 
            alpha = 1 - 0.2*(is.na(lin) | div | !is.na(cllab))) %>%
            mutate(col=factor(col),
                col_lab=factor(col_lab))

    p <- ggtree(tree,mrsd="2020-1-1", lwd=.8) %<+% lin_df +
        aes(color=col, alpha=alpha) + 
        geom_nodelab(mapping=aes(node=node, label=cllab, color=col_lab, x=2020+2.5, alpha=1.0),size=5.0, hjust=0)+ 
        scale_color_manual(
            breaks=factor(0L:3L),
            values=c("slategray3","steelblue4", "steelblue2", "tomato3"), 
            na.value="gray15")+
        vexpand(.02, direction=-1) +
        hexpand(.03, direction=1) + 
        scale_x_ggtree()+
        labs(x="Year")+
        scale_x_continuous(breaks=seq(from=1950, to=2020, by= 10)) +
        coord_cartesian(xlim = c(1950, 2026))+
        theme(legend.position = "none",
            axis.text.x = element_text(size=15),
            axis.title.x = element_text(size=20),
            axis.line.x = element_line()) #+
        #scale_x_ggtree()#,
            #plot.margin = unit(c(14,14,8,8), "mm"))
    p
}

plot_lin_panel <- function(lin_data) 
{
    make_motif_df(lin_data) %>% 
        pivot_longer(cols=!lin, names_to ="site", values_to="motif") %>%
        group_by(site) %>% 
        mutate(motif_id = as.integer(factor(motif))) %>%
        mutate(motif_id = factor(ifelse(motif %in% c("WT", "penA_2", "No"), NA, motif_id))) %>%
        ggplot(aes(x=lin, y=site)) +
            geom_tile(aes(fill=motif_id),color="gray25") +
            geom_text(aes(label=motif), size=2.5,angle=90) +
            scale_fill_brewer(palette ="Paired", na.value="gray95")+
            theme_minimal() + 
            labs(x="Lineage", y="Site")+
            theme(
                axis.title.x=element_text(size = 10.0),
                axis.text.y=element_text(size=8),
                legend.position = "none", 
                panel.grid = element_blank(),
                plot.margin=margin(t=0))
}


plot_lin_panel <- function(lin_data) 
{
    make_motif_df(lin_data) 
}


plot_motif_freqs <- function(lin_data) 
{
    make_motif_df(lin_data) %>% 
        pivot_longer(cols=!lin, names_to ="site", values_to="motif") %>%
        group_by(site) %>% 
        mutate(motif_id = as.integer(factor(motif))) %>%
        mutate(motif_id = factor(ifelse(motif %in% c("WT", "penA_2", "No"), NA, motif_id))) %>%
        ggplot(aes(y=site)) +
            geom_bar(aes(fill=motif_id)) +
            scale_fill_brewer(palette ="Paired", na.value="gray95")+
            theme_minimal() + 
            labs(x="Marginal Frequency")+
            theme(legend.position = "none", 
                panel.grid = element_blank(),
                axis.title.y=element_blank(),
                axis.title.x=element_text(size=10.0),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                plot.margin=margin(t=0,b=0,l=1,r=0))
}
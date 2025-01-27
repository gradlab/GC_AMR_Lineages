#Plot lineage cluster figure 
plot_lin_clust <- function(lins, anc_parents, lin_ids, nodedates, tree, allele_motif, mics)
{
    cutoffs <- c(ciprofloxacin=1, cefixime=0.125, azithromycin=4)

    tree_combined <- keep.tip(tree, unlist(sapply(lins, function(x) x$tip.label)))
    max_t <- max(nodedates[tree_combined$tip.label])
    n_tip <- length(tree_combined$tip.label)

    lin_df2 <- data.frame(node=1:(tree_combined$Nnode+length(tree_combined$tip.label)), lin=NA, div=T, cllab=NA)
    rownames(lin_df2) <- c(tree_combined$tip.label, tree_combined$node.label)

    lins_filtered <- lapply(seq_along(lins), function(x) filter_lin(lins[[x]], anc_parents[[x]]))

    lin_id <- 1
    for(li in seq_along(lins))
    {
        lin <- lins[[li]]
        non_div <- lins_filtered[[li]]

        div_tips <- lin$tip.label[!(lin$tip.label %in% non_div$tip.label)]
        non_div_nd <- drop.tip(lin, div_tips, trim.internal=T, collapse.singles=F)

        ids_non_div <- c(non_div_nd$tip.label, non_div_nd$node.label)
        ids <- c(lin$tip.label, lin$node.label)
        
        lin_root <- lin$node.label[1]

        lin_df2[ids, "lin"] <- lin_ids[li]
        lin_df2[ids_non_div, "div"] <- F
        lin_df2[lin_root, "cllab"] <- paste0("Lineage ", lin_id)
    }

    #Fill in blanks in colouring
    preord <- rev(postorder(tree_combined))
    for (i in preord)
    {
        e <- tree_combined$edge[i, ]
        cname <- rownames(lin_df2)[e[2]]
        pname <- rownames(lin_df2)[e[1]]

        if(is.na(lin_df2[cname, "lin"]) && !is.na(lin_df2[pname, "lin"]))
        {
            lin_df2[cname, "lin"] <- lin_df2[pname, "lin"]
        } 

        if(all(!is.na(lin_df2[c(cname,pname), "lin"])) && (!lin_df2[cname, "div"] && lin_df2[pname, "div"]) &&
            lin_df2[cname, "lin"] == lin_df2[pname, "lin"])
        {
            lin_df2[pname, "div"] <- lin_df2[cname, "div"]
        }
    }

    bnpr_df <- bnpr_as_df(lins_filtered, nodedates) %>% 
        as_tibble() %>%
        mutate(lineage = sapply(lineage, function (x) lin_ids[[x]])) %>%
        mutate(lineage = factor(lineage)) 

    min_t <- min(bnpr_df$t)
    #max_t <- max(bnpr_df$t)

    lin_df2 <- lin_df2 %>% 
        mutate(col = ifelse(div, NA, lin), 
            alpha = 1.0 - 0.1*(div | !is.na(cllab))) %>%
        mutate(col = ifelse(!is.na(cllab), NA, col)) %>%
        mutate(col = factor(col))

    tree_panel <- ggtree(tree_combined, mrsd=paste0(as_date(date_decimal(max_t)))) %<+% lin_df2 +
        aes(color=col, alpha=alpha) + 
        scale_color_brewer(palette ="Set2",guide="none",na.value="gray15") +
        scale_alpha_continuous(guide="none", range=c(.1,1.0)) + 
        geom_vline(aes(xintercept = 2007), linetype="dashed", color="gray25", alpha=0.8) +
        geom_vline(aes(xintercept = 2011), linetype="dashed", color="gray25", alpha=0.8) +
        scale_x_continuous(limits=c(min_t-15, max_t), expand = c(0, 0)) +
        #vexpand(.1,direction=-1) + 
        theme_tree2() +
        theme(axis.line.y = element_line(),
            axis.text.x = element_blank(),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank())

    tnames <- rev(get_taxa_name(tree_panel))

    old_names <- c("mtr_c: GC_allele",
        "mtr_c: LOF",
        "mtr_promoter: Adel_13bp_invertedRepeat",
        "mtr_promoter: WT",
        "mtr_r: GC_allele",
        "mtr_r: LOF",
        "penA: penA_2",
        "penA: penA_34",
        "penA: penA_5",
        "pbp: 421L",
        "pbp: 421P",
        "rRNA_23s_C2611T: T",
        "rRNA_23s_C2611T: C"
        )

    new_names <- c("mtrC: non-m",
        "mtrC: LOF",
        "mtr promoter: A-del",
        "mtr promoter: non-m",
        "mtrR: non-m",
        "mtrR: LOF",
        "penA: 2",
        "penA: mosaic",
        "penA: 5",
        "ponA: 421L",
        "ponA: 421P",
        "23s rRNA: 2611T",
        "23s rRNA: 2611C"
    )

    names(new_names) <- old_names
    rename_fn <- function(x) ifelse(x %in% old_names, new_names[x], x)

    motif_panel <- allele_motif %>% 
        filter(x %in% tnames) %>%
        mutate(x=factor(x, levels = tnames, ordered=TRUE)) %>%
        pivot_longer(!x, names_to = "motif", values_to = "value") %>%
        group_by(motif) %>% 
        filter(min(table(value)) > n_tip/40 && max(table(value)) < 39/40 * n_tip) %>%
        ungroup() %>%
        mutate(Presence=ifelse(value, "present", "absent")) %>%
        mutate(motif = rename_fn(motif)) %>%
        filter(!str_detect(motif, "(penA: 5)|(penA: 2)|(mtr[CD].*)"))%>%
        ggplot(aes(y=x, x=reorder(factor(motif), ifelse(str_detect(motif, ".*23S.*"), -(1:length(motif)), (1:length(motif))), order=T), fill=Presence), color=NA,width=1) +
            geom_tile() +
            scale_fill_manual(values=c("gray85", "gray25")) + 
            scale_x_discrete(position = "top") +
            theme_bw() + 
            theme(axis.title.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                axis.text.x = element_text(size=rel(1.8), angle = 90, hjust=0),
                axis.title.x = element_blank(),
                axis.line.x = element_blank(),
                axis.line.y = element_blank(),
                panel.border = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank())

    mic_panel <- mics %>% 
        filter(x %in% tnames) %>%
        mutate(x=factor(x, levels = tnames, ordered=TRUE)) %>%
        mutate(value=(value<=cutoffs[name])) %>%
        mutate(name=str_to_title(name)) %>%
        ggplot(aes(y=x, x=name, fill=value), color=NA,width=1) +
            geom_tile() +
            scale_fill_manual(breaks=c(FALSE, TRUE), 
                values=c("darkred", "steelblue4"), 
                name = "MIC <= Cutoff") + 
            scale_x_discrete(position = "top", limits=c("Ciprofloxacin", "Cefixime", "Azithromycin")) +
            theme_bw() + 
            theme(axis.title.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                axis.text.x = element_text(size=rel(1.8), angle = 90, hjust=0),
                axis.title.x = element_blank(),
                axis.line.x = element_blank(),
                axis.line.y = element_blank(),
                panel.border = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank())

    bnpr_panel <- ggplot(bnpr_df, aes(x=t, y=y, color=lineage, fill=lineage)) + 
        geom_line() + 
        geom_ribbon(aes(ymin = yhi, ymax = ylo), alpha=.5) +
        scale_color_brewer(palette ="Set2") +
        scale_fill_brewer(palette ="Set2") +
        geom_vline(aes(xintercept = 2007), linetype="dashed", color="gray25", alpha=0.8) +
        geom_vline(aes(xintercept = 2011), linetype="dashed", color="gray25", alpha=0.8) +
        scale_y_log10() +
        scale_x_continuous(limits=c(min_t-15, max_t), expand = c(0, 0)) +
        labs(x="Year", y="Effective Population Size") + 
        theme_bw() +
        theme(#axis.title.x = element_blank(),
                #axis.ticks.x = element_blank(),
                axis.text.x = element_text(size=rel(1.8)),
                axis.line.x = element_line(),
                axis.line.y = element_line(),
                axis.text.y = element_text(size=rel(1.8)),
                axis.title.y = element_text(size=rel(2.2)),
                axis.title.x = element_text(size=rel(2.2)),
                panel.border = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank())

    p <- tree_panel + ((motif_panel | mic_panel) + plot_layout(widths=c(2,1))) + bnpr_panel + guide_area()
    p + plot_layout(ncol=2, widths=c(4,2), heights = c(5,5),guides = 'collect')
}

bnpr_as_df <- function(lins_filtered, nodedates)
{
    bnpr_data <- list()

    for (i in seq_along(lins_filtered))
    {
        li <- lins_filtered[[i]]

        max_t_f <- max(nodedates[li$tip.label])
        #offs_f <- max_t_f - max_t

        bn <- BNPR(li)
        grid <- bn$grid
        xlim <- c(max(grid), min(grid))
        mask <- bn$x >= min(xlim) & bn$x <= max(xlim)
        t <- max_t_f - bn$x[mask]
        y <- bn$effpop[mask]
        yhi <- bn$effpop975[mask]
        ylo <- bn$effpop025[mask]
        
        bnpr_data[[i]] <- data.frame(t=t, y=y, yhi=yhi, ylo=ylo, lineage=i)
    }
    do.call(rbind, bnpr_data)
}

plot_count_msw <- function(lin, meta)
{
    meta[lin$tip.label, c("date","sexual_behavior")] %>%
        mutate(sexual_behavior=factor(sapply(sexual_behavior, function (x) if(x %in%  c("MSM", "MSW", "MSMW")) x else NA))) %>%
        mutate(date=as.integer(str_extract(date,"\\d\\d\\d\\d"))) %>% 
        mutate(date_d=3L*as.integer(date / 3)) %>%
        filter(!is.na(sexual_behavior)) %>%
        group_by(date_d) %>%
        count(sexual_behavior) %>%
        ungroup() %>% 
        complete(sexual_behavior, date_d, fill=list(n=0)) %>% 
        group_by(date_d) %>%
        mutate(prop=n/sum(n)) %>% 
        ggplot(aes(x=date_d, y=prop, fill=sexual_behavior)) + 
            geom_area() + 
            scale_fill_brewer(palette="Set2") +
            labs(x="Year", y="Proportion", fill="Sexual Behaviour") +
            theme_minimal()
}

plot_count_msw2 <- function(lins, meta)
{
   # all_tips <- unlist(lapply(lins, function(x) x$tip.label))

    meta[lin$tip.label, c("date","sexual_behavior")] %>%
        mutate(sexual_behavior=factor(sapply(sexual_behavior, function (x) if(x %in%  c("MSM", "MSW", "MSMW")) x else NA))) %>%
        mutate(date=as.integer(str_extract(date,"\\d\\d\\d\\d"))) %>% 
        mutate(date_d=5L*as.integer(date / 5)) %>%
        filter(!is.na(sexual_behavior)) %>%
        group_by(date_d) %>%
        mutate(tot=n()) %>%
        ungroup() %>%
        group_by(sexual_behavior, date_d) %>%
        summarise(n=ifelse(tot==0 | is.infinite(n()), 0, n()/tot), .groups = "keep") %>% 
        ggplot(aes(x=date_d, y=n, fill=sexual_behavior)) +
            geom_area() + 
            scale_fill_brewer(palette="Set2") +
            labs(x="Year", y="Count", fill="Sexual Behaviour") +
            theme_minimal()
}
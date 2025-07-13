library(ape)
library(ggplot2)
library(ggtree)
library(stringr)
library(phylodyn)
library(viridis)
library(tidyverse)
library(ggnewscale)
library(patchwork)

source("../plotting/plot_lin_motif.R")
source("extract_penA_snps.R")

size_cut <- 30 #Lineage size cutoff
time_cut <- 1980 #Lineage tMRCA cutoff

resis_alleles <- read.table("metadata/2024-02-28_gc_resistance_alleles.tsv", header=T)
mtr_alleles <- read.table("metadata/mtr_alleles.tsv", header=T)
penA_alelles <- read.table("metadata/gc_penA.txt", header=T)
pres_abs <- read.table("metadata/resistance_gene_presence_absence.tsv", header=T)
penA_desc <- read.delim("metadata/penA_descriptions.txt")
rownames(penA_desc) <- penA_desc$id

alelles_pml <-  read.csv("pastml/alleles.csv",header=T)
rownames(alelles_pml) <- alelles_pml$wgs_id

mics <- read.csv("metadata/mic_clean.csv")
rownames(mics) <- mics[,1]

azilog2<-log2(as.numeric(mics$AZI))
mics_azi <- data.frame(AZI=ifelse(is.infinite(azilog2), NA, azilog2))
rownames(mics_azi) <- rownames(mics)

meta <- data.frame(read.delim("metadata/Ng-Combined-Metadata.txt", 
    as.is=T, colClasses = "character"))
meta <- meta[!is.na(meta$wgs_id), ]
rownames(meta) <- meta$wgs_id

rrna_alleles_1 <- read.table("metadata/rRNA_allele_summary_grad.txt", header=T)
rrna_alelles_2 <- read.table("metadata/rRNA_allele_summary_rest.txt", header=T)
rrna <-rbind(rrna_alleles_1, rrna_alelles_2) 
rrna$wgs_id <- meta[match(rrna$wgs_id, meta$accession), "wgs_id"]
rrna <- rrna[!is.na(rrna$wgs_id), ]
rownames(rrna) <- rrna$wgs_id

all_ids <- resis_alleles$wgs_id
rownames(resis_alleles) <- resis_alleles$wgs_id
resis_alleles$wgs_id <- NULL

rownames(mtr_alleles) <- mtr_alleles$wgs_id
mtr_alleles$wgs_id <- NULL

rownames(penA_alelles) <- penA_alelles$wgs_id
penA_alelles$wgs_id <- NULL

anc_states_combined <- read.table("pastml/combined_ancestral_states_gyrA_filtered.tab", header=T) %>% 
    as_tibble %>% 
    left_join(read.table("pastml/combined_ancestral_states_mtr_filtered.tab", header=T) %>% as_tibble, 
        by="node") %>% 
    left_join(read.table("pastml/combined_ancestral_states_rRNA_filteredDELTRAN.tab", header=T) %>% as_tibble, 
        by="node") %>% 
    as.data.frame

anc_states_combined_with_penA <- anc_states_combined
rownames(anc_states_combined_with_penA) <- anc_states_combined_with_penA$node
anc_states_combined_with_penA$node <- NULL

penA_resis_motifs <- extract_penA_snps(anc_states_combined %>% select(penA), penA_desc)
anc_states_combined <- cbind(anc_states_combined, penA_resis_motifs)
 
rownames(anc_states_combined) <- anc_states_combined$node
anc_states_combined$node <- NULL
anc_states_combined$penA <- NULL
rownames(pres_abs) <- pres_abs$wgs_id
pres_abs$wgs_id <- NULL

mtr_names <- c("mtrPromoter", "mtrR", "mtrC", "mtrD")
gyrA_names <- c("GyrA_91", "GyrA_95")
parC_names <- c("ParC_86", "ParC_87", "ParC_91")
penA_names <- c("penA")
pbp_names <- c("PBP1_421")

timing_data <- readRDS("pastml/timing_data.rds")
tree_retained <- read.tree("pastml/tree_retained.nwk")
record <- timing_data$record
nodedates <- timing_data$nodedates

meta <- meta[tree_retained$tip.label, ]
sexbeh <- data.frame(sex=factor( sapply(meta$sexual_behavior, function (x) if(x %in%  c("MSM", "MSW", "MSMW")) x else NA)))
rownames(sexbeh) <- rownames(meta)

retained_names <- c(tree_retained$tip.label, tree_retained$node.label)
n_tip_a <- length(tree_retained$tip.label)

## Lineage stem lookup vector
stems <- rep(F, length(retained_names))
names(stems) <- retained_names

## Compute probs that a node is younger than time_cut
p_younger <- apply(record, 2, function(x) sum(x > time_cut))/nrow(record)

## Mark all oldest nodes that are younger than time_cut with probability >= .99 as potential stems
for (i in 1:nrow(tree_retained$edge))
{
    e <- tree_retained$edge[i,]
    stems[e[2]] <- ((p_younger[retained_names[e[1]]] < .99) && (p_younger[retained_names[e[2]]] >= .99))
}

## Generate subtrees younger than time cutoff and larger than sizecut
## Plot stem locations to check correctness
df_stems <- data.frame(node=1:(tree_retained$Nnode+n_tip_a), labs=retained_names, is_stem=NA)
rownames(df_stems) <- c(tree_retained$tip.label, tree_retained$node.label)
df_stems$is_stem <- sapply(1:nrow(df_stems), function(x) stems[df_stems$labs[x]])

stem_forest <- list()
for(i in 1:nrow(df_stems))
{
    if (df_stems$is_stem[i] && i > n_tip_a)
    {
        x <- extract.clade(tree_retained, df_stems$lab[i])
        if (length(x$tip.label) > size_cut)
        {
            stem_forest[[length(stem_forest) + 1]] <- x
        }
    } 
}

## For each subtree get ancestral state reconstruction from pastml result
forest_ancrec <- list()
for (i in seq_along(stem_forest))
{
    tr <- stem_forest[[i]]

    tr_nodes <- c(tr$tip.label, tr$node.label)

    forest_ancrec[[length(forest_ancrec) + 1]] <- list(
        tree=tr,
        anc_states=anc_states_combined[tr_nodes, ])
}

## Assign lineages based on subtrees and ancestral state reconstruction
assign_lins <- function(anc_tree, size_cutoff)
{
    phy <- anc_tree$tree
    anc_states <- anc_tree$anc_states

    n_tip <- length(phy$tip.label)
    n <- phy$Nnode + n_tip

    root_idx <- n_tip+1

    ord <- postorder(phy)
    n_contig <- c(rep(1L, n_tip), rep(0L, phy$Nnode))
    is_lin_anc <- rep(F, n)

    no_change <- rep(nrow(phy$edge), T)
    no_change <- no_change & apply(phy$edge, 1, function (x) all(anc_states[x[1],] == anc_states[x[2],]))
    no_change <- ifelse(is.na(no_change), FALSE, no_change)

    for (i in ord)
    {
        pa <- phy$edge[i,1]
        ch <- phy$edge[i,2]

        if (no_change[i])
        {
            n_contig[pa] <- n_contig[pa] + n_contig[ch]
        } else if (n_contig[ch] > 0) {
            is_lin_anc[ch] <- T
        }
    }
    
    if(n_contig[root_idx] > 0)
    {
        is_lin_anc[root_idx] <- T
    }

    lin_clades <- as.list(rep(list(list()),n))

    lin_MRCAs <- c()

    for(i in 1:(length(ord)+1))
    {
        node <- ifelse(i > length(ord), root_idx, phy$edge[ord[i],2])  
        if (is_lin_anc[node] && node > n_tip)
        {  
            cl <- extract.clade(phy, node)
            cl_size <- n_contig[node]

            if (cl_size > size_cutoff)
            {
                lin_MRCAs <- c(lin_MRCAs, cl$node.label[1])
            } 
        }
    }

    ### Reduce fragmentation of lineages by defining them using MRCAs
    nodes_postorder <- c(phy$edge[postorder(phy),2], n_tip + 1)
    mrcas_in_order <- nodes_postorder[sapply(nodes_postorder, 
        function (x) if (x > n_tip) phy$node.label[x-n_tip] %in% lin_MRCAs else FALSE
    )]  
    lins_refined <- lapply(mrcas_in_order, function(x) extract.clade(phy, x))
    for(i in seq_along(lins_refined))
    {
        if (i > 1)
        {
            for(j in 1:(i-1))
            {
                lins_refined[[i]] <- drop.tip(lins_refined[[i]], lins_refined[[j]]$tip.label) 
                if(is.null(lins_refined[[i]])) break
            }
        }
    }

    return(list(clades=lins_refined, n_contig=n_contig, is_lin_anc=is_lin_anc, no_change=no_change))
}

# Function to remove divergent lines of descent from lineages defined through MRCA
filter_lin <- function(x, anc_parent)
{
    iset <- c(anc_parent$tree$tip.label, anc_parent$tree$node.label)
    mrca_idx <- which(iset == x$node.label[1])
    anc_st <- anc_parent$anc_states[mrca_idx, ]
    div_desc <- c()
    for (nn in x$node.label){
        n_idx <-  which(iset == nn)
        cond <- anc_st == anc_parent$anc_states[n_idx, ]
        cond <- ifelse(is.na(cond), FALSE, cond)
        if(!all(cond))
        {
            div_desc <- union(div_desc, extract.clade(x, nn)$tip.label)
        }
    }

    div_tips <- sapply(x$tip.label,   
        function(tl) {
            ti <- which(iset == tl)
            cond <- anc_st == anc_parent$anc_states[ti, ]
            cond <- ifelse(is.na(cond), FALSE, cond)
            return(!all(cond))
        }
    )

    div_desc <- union(div_desc, x$tip.label[div_tips])
    print(length(div_desc) / length(x$tip.label))
    drop.tip(x, div_desc) 
}

lin_assignments_raw <- lapply(forest_ancrec, function(x)
    c(anc_parent = list(x), 
    lin_assignment = list(assign_lins(size_cutoff=size_cut, x))))

## Filter out subtrees with no assigned lineages
lin_assignments <- list()
for (la in lin_assignments_raw)
{
    if(length(la$lin_assignment$clades) > 0)
    {
        lin_assignments[[length(lin_assignments)+1]] <- la
    }
}

# Prepare dataframe to visualise the lineage assignment
lin_df <- data.frame(node=1:(tree_retained$Nnode+length(tree_retained$tip.label)), lin=NA, anc=FALSE, div=T, cllab=NA)
rownames(lin_df) <- c(tree_retained$tip.label, tree_retained$node.label)

lin_id <- 1
for(la in lin_assignments)
{
    subtr_ids <- c(la$anc_parent$tree$tip.label, la$anc_parent$tree$node.label)
    lin_df[subtr_ids, "anc"] <- la$lin_assignment$is_lin_anc
    for (lin in la$lin_assignment$clades)
    {
        non_div <- filter_lin(lin, la$anc_parent)
        
        div_tips <- lin$tip.label[!(lin$tip.label %in% non_div$tip.label)]
        non_div_nd <- drop.tip(lin, div_tips, trim.internal=T, collapse.singles=F)
        ids_non_div <- c(non_div_nd$tip.label, non_div_nd$node.label)

        ids <- c(lin$tip.label, lin$node.label)
        lin_root <- lin$node.label[1]
        lin_df[ids, "lin"] <- lin_id
        lin_df[ids_non_div, "div"] <- F
        lin_df[lin_root, "cllab"] <- paste0("Lineage ", lin_id)
        lin_id <- lin_id + 1
    }
}

#Fill in blanks in colouring
preord <- rev(postorder(tree_retained))
for (i in preord)
{
    e <- tree_retained$edge[i, ]
    cname <- rownames(lin_df)[e[2]]
    pname <- rownames(lin_df)[e[1]]

    if(is.na(lin_df[cname, "lin"]) && !is.na(lin_df[pname, "lin"]))
    {
        lin_df[cname, "lin"] <- lin_df[pname, "lin"]
    }

    if((!lin_df[cname, "div"] && lin_df[pname, "div"]) &&
    all(!is.na(lin_df[c(cname, pname), "lin"] )) &&
    lin_df[cname, "lin"] == lin_df[pname, "lin"])
    {
        lin_df[pname, "div"] <- lin_df[cname, "div"]
    }
}

lin_df$lin <- factor(lin_df$lin)

# Lineage preview plotting function
plot_func <- function(x, anc_parent)
{
    major_pena_names <- union(names(sort(table(penA_alelles), T))[1:10], "unknown")
    major_penA <- sapply(tree_retained$tip.label, function(y) 
    {
        if(penA_alelles[y, 1] %in% major_pena_names) penA_alelles[y, 1] else "Other"
    })

    max_t <- max(nodedates[tree_retained$tip.label])
    max_t_lin <- max(nodedates[x$tip.label])
    offs <- max_t_lin - max_t
    
    lin_filtered <- filter_lin(x, anc_parent)
    max_t_f <- max(nodedates[lin_filtered$tip.label])
    offs_f <- max_t_f - max_t

    alleles_merged <- cbind(resis_alleles[tree_retained$tip.label, c(gyrA_names,parC_names,pbp_names)], 
        mtr_alleles[tree_retained$tip.label, mtr_names], penA=major_penA, rRNA_23s_C2611T=rrna[tree_retained$tip.label, 'X23S_C2611T'])
    alleles_merged2 <- alleles_merged %>%
        mutate(across(everything(), ~ as.data.frame.matrix(table(row_number(), .x) *   
            NA^(is.na(.x)) > 0))) %>% 
        unpack(where(is.data.frame), names_sep = ": ") %>% as.data.frame()
    rownames(alleles_merged2) <- rownames(alleles_merged)

    p1 <- ggtree(x,ladderize=T, mrsd=paste0(as_date(date_decimal(max_t_lin)))) + 
        theme_tree2() + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1.0, size=6.0)) +
        ggtitle(names(sort(table(penA_alelles[x$tip.label,1]), decreasing = T))[1]) + 
        scale_x_ggtree(breaks=2020:time_cut)
    p2 <- gheatmap(p1, alleles_merged2, offset=rel(8), width=1.0, colnames=F, legend_title="Alleles",color=NA) + 
        scale_fill_viridis_d(option="C",limits=c(TRUE,FALSE)) + scale_x_ggtree()

    p2 <- p2 + new_scale_fill()
    p2 <- gheatmap(p2, mics_azi, offset=rel(2), width=0.1, colnames=F, legend_title="AZI_MIC",color=NA) + 
        scale_fill_viridis(option="D", limits=c(-7,8))

    plot(p2)

    plot_BNPR(BNPR(lin_filtered), xlim=c(2020-time_cut,offs_f), xaxt='n')
    axis(1, at = c((2020-time_cut):0), labels=c(time_cut:2020))
}

# Plot all lineage previews 
pdf("trees.pdf",10,10)
lin_id <-1
for (la in lin_assignments) 
{
    for (lin in la$lin_assignment$clades)
    {
        print(paste("Lineage", lin_id))
        plot_func(lin, la$anc_parent)
        lin_id<-lin_id+1
    }
}
dev.off()

# Function to extract lineage data
get_lin_data <- function(x, anc_parent)
{
    lin_filtered <- filter_lin(x, anc_parent)

    max_t_lin <- max(nodedates[lin_filtered$tip.label])
    mrst <- 2020 - max_t_lin

    lin_alleles <- anc_states_combined[lin_filtered$node.label[1], ]

    lin_mics <- meta[lin_filtered$tip.label, ] %>% 
        select(c(ceftriaxone, ciprofloxacin, cefixime, azithromycin)) %>% 
        pivot_longer(everything()) %>% 
        mutate(value = as.numeric(value))

    sex_beh <- sexbeh[lin_filtered$tip.label, ]

    list(phylo=lin_filtered, 
        mrst=mrst, 
        alleles=lin_alleles, 
        lin_mics=lin_mics, 
        sex_beh=sex_beh, 
        dates=as.numeric(str_extract(meta[lin_filtered$tip.label, "date"], "^\\d{4}"))
        )
}

# Extract data for all lineage
lin_data <- list()
for (la in lin_assignments) 
{
    for (lin in la$lin_assignment$clades)
    {
        lin_data[[length(lin_data)+1]] <- get_lin_data(lin, la$anc_parent)
    }
}

# Plot lineage summary
source("../plotting/plot_lin_motif.R")

mic_for_clust <- meta[tree_retained$tip.label, ] %>% 
        select(c(ciprofloxacin, cefixime, azithromycin)) %>% 
        mutate(x=tree_retained$tip.label) %>%
        pivot_longer(!c(x)) %>% 
        mutate(value = as.numeric(value))

ltrees <- plot_lin_tree(tree_retained, lin_df, mic_for_clust)
pdf("../plots/lin_tree.pdf", 12, 16)
plot(ltrees$p1)
dev.off()

pdf("../plots/lin_tree_mic.pdf", 12, 16)
plot(ltrees$p2)
dev.off()

make_motif_df(lin_data) %>% 
    as_latex() %>% 
    writeLines("../tables/determinant_table.tex")


lins_unfiltered <- list(lins=list(), anc_parent=list(), lin_id=list())
lin_id <-1
for (la in lin_assignments) 
{
    for (lin in la$lin_assignment$clades)
    {
        lins_unfiltered$lins[[lin_id]] <- lin
        lins_unfiltered$anc_parent[[lin_id]] <- la$anc_parent
        lins_unfiltered$lin_id[[lin_id]] <- lin_id
        lin_id<-lin_id+1
    }
}

allele_motif <- anc_states_combined_with_penA[tree_retained$tip.label, ] %>%
        mutate(across(everything(), ~ as.data.frame.matrix(table(row_number(), .x) *   
            NA^(is.na(.x)) > 0))) %>% 
        unpack(where(is.data.frame), names_sep = ": ") %>% as.data.frame()
allele_motif$x <- tree_retained$tip.label

source("../plotting/plot_lin_clust.R")
lins_to_extract <- c(20:22)

#Lineage cluster plot
p <- plot_lin_clust(lapply(lins_to_extract, function(x) lins_unfiltered$lins[[x]]), 
    lapply(lins_to_extract, function(x) lins_unfiltered$anc_parent[[x]]),
    lapply(lins_to_extract, function(x) lins_unfiltered$lin_id[[x]]), nodedates, tree_retained, allele_motif, mic_for_clust)

pdf("../plots/lin_clust.pdf", 13, 13)
plot(p)
dev.off()

#Lineage 20 sexual beahviours
pdf("../plots/sex_beh.pdf", 13, 13)
plot(plot_count_msw(lins_unfiltered$lins[[20]], meta))
dev.off()

#Run linear model to verify that MICs for 91F/95A are indistinguishable from those of 91F/95G
dat4lm <- alelles_pml %>% filter(gyrA %in% c("91F/95A", "91F/95G"), 
        parC == "86D/87R/91E") %>%
        select(wgs_id, gyrA, mtr_r, mtr_c, mtr_d, mtr_promoter) %>%
        mutate(mic = meta[wgs_id,"ciprofloxacin"]) %>%
        filter(!is.na(mic) & mic > 0) %>%
        mutate(mic=log2(as.numeric(mic))) 

dat4lm$mtr_r <- relevel(factor(dat4lm$mtr_r), "GC_allele")
dat4lm$mtr_c <- relevel(factor(dat4lm$mtr_c), "GC_allele")
dat4lm$mtr_d <- relevel(factor(dat4lm$mtr_d), "GC_allele")
dat4lm$mtr_promoter <- relevel(factor(dat4lm$mtr_promoter), "WT")

lm(data=dat4lm, formula=mic~gyrA+ mtr_r+mtr_c+mtr_d+mtr_promoter) %>% 
    summary() %>% 
    print()


# Write lineage data
saveRDS(list(lin_data=lin_data), "lineages.rds")

rrna_df <- data.frame(node=1:(tree_retained$Nnode+length(tree_retained$tip.label)))
rrna_df$rrna <- anc_states_combined[c(tree_retained$tip.label,tree_retained$node.label), "rRNA_23s_C2611T"]

# Plot distribution of 23s mutations
pdf("../plots/rrna_tree.pdf",10,10)
p <- ggtree(tree_retained,mrsd="2020-1-1", lwd=.8) %<+% (rrna_df %>% rename(rRNA_23S_C2611T=rrna)) +
        aes(color=rRNA_23S_C2611T) 
p <- gheatmap(p, rrna %>% select(X23S_C2611T), offset=rel(8), width=0.1, colnames=F, legend_title="rrna copies",color=NA) + 
        scale_fill_viridis(option="C") + labs(fill = "rRNA_23s_C2611T Copy Count") + scale_x_ggtree()
plot(p)
dev.off()
library(ape)
library(stringr)
library(tidyverse)
library(lubridate)

#Read inputs 
resis_alleles <- read.table("../metadata/2024-02-28_gc_resistance_alleles.tsv", header=T)
mtr_alleles <- read.table("../metadata/mtr_alleles.tsv", header=T)
penA_alelles <- read.table("../metadata/gc_penA.txt", header=T)
timing_res <- readRDS("../big_timed_filtered.rds")

rrna_alleles_1 <- read.table("../metadata/rRNA_allele_summary_grad.txt", header=T)
rrna_alelles_2 <- read.table("../metadata/rRNA_allele_summary_rest.txt", header=T)
rrna <- rbind(rrna_alleles_1, rrna_alelles_2) 
meta <- data.frame(read.delim("../metadata/Ng-Combined-Metadata.txt"))
rrna$wgs_id <- meta[match(rrna$wgs_id, meta$accession), "wgs_id"]

rrna <- rrna[!is.na(rrna$wgs_id), ]
rownames(rrna) <- rrna$wgs_id
rrna[, 2:4] <- rrna[, 2:4] > 0
rrna[, "X23S_C2611T"] <- ifelse(rrna[, "X23S_C2611T"], "T", "C")

all_ids <- resis_alleles$wgs_id
rownames(resis_alleles) <- resis_alleles$wgs_id
resis_alleles$wgs_id <- NULL

rownames(mtr_alleles) <- mtr_alleles$wgs_id
mtr_alleles$wgs_id <- NULL

rownames(penA_alelles) <- penA_alelles$wgs_id
penA_alelles$wgs_id <- NULL

mtr_names <- c("mtrPromoter", "mtrR", "mtrC", "mtrD")
gyrA_names <- c("GyrA_91", "GyrA_95")
parC_names <- c("ParC_86", "ParC_87", "ParC_91")
penA_names <- c("penA")
pbp_names <- c("PBP1_421")

alleles <- data.frame(wgs_id = all_ids)
is_missing <- function(x) sapply(x, function (y) is.na(y) | any(str_detect(y, c("(?i)unknown","(?i)duplicate"))))
make_prefix <- function(x)
{
    y <- str_extract(str_extract(x, "_\\d+"),"\\d+")
    ifelse(is.na(y), "", y)
}
splat_allele <- function(names, df)
{
    function(x) ifelse(all(!is_missing(df[x, names])), paste(paste(make_prefix(names), df[x, names], sep=""), collapse="/"), NA)
}

alleles <- within(alleles, 
    gyrA <- sapply(wgs_id, splat_allele(gyrA_names, resis_alleles))
)

alleles <- within(alleles, 
    parC <- sapply(wgs_id, splat_allele(parC_names, resis_alleles))
)

alleles <- within(alleles, 
    pbp <- sapply(wgs_id, splat_allele(pbp_names, resis_alleles))
)

alleles <- within(alleles, 
    penA <- sapply(wgs_id, splat_allele(penA_names, penA_alelles))
)

alleles <- within(alleles, 
    mtr_promoter <- sapply(wgs_id, splat_allele(mtr_names[1], mtr_alleles))
)
alleles <- within(alleles, 
    mtr_r <- sapply(wgs_id, splat_allele(mtr_names[2], mtr_alleles))
)
alleles <- within(alleles, 
    mtr_c <- sapply(wgs_id, splat_allele(mtr_names[3], mtr_alleles))
)
alleles <- within(alleles, 
    mtr_d <- sapply(wgs_id, splat_allele(mtr_names[4], mtr_alleles))
)

alleles <- within(alleles,
    rRNA_23s_C2611T <- sapply(wgs_id, splat_allele("X23S_C2611T", rrna))
)

rownames(alleles) <- alleles$wgs_id
attributes(timing_res$tree)$order <- NULL

record <- timing_res$record
tree_dated <- timing_res$tree
tree_retained <- makeNodeLabel(tree_dated)

# Ensure that tree is consistent & filter out 
## Find best root rows as returned by bactdating
n1 <- tree_retained$Nnode + length(tree_retained$tip.label)
bestroot <- as.numeric(names(sort(table(record[floor(nrow(record) / 2):nrow(record), 'root']),decreasing=T)[1]))
bestrows <- intersect(floor(nrow(record) / 2):nrow(record), which(record[, 'root']==bestroot))
record <- record[bestrows, 1:n1]
colnames(record) <- c(tree_retained$tip.label, tree_retained$node.label)

write.tree(tree_retained, "tree_retained.nwk")
t2 <- tree_retained
tree_retained <- read.tree("tree_retained.nwk")

## Check RF distance to ensure that the tree is the one we expect
stopifnot(dist.topo(unroot(t2),unroot(tree_retained)) < 1e-6)
stopifnot(dist.topo(unroot(t2),unroot(tree_retained), 'score') < 1e-6)
stopifnot(all(abs(sort(node.depth.edgelength(tree_retained)) - sort(node.depth.edgelength(t2))) < 1e-6))

## Remove tips with missing typing data
nodedates <- timing_res$rootdate + node.depth.edgelength(tree_retained)
names(nodedates) <- c(tree_retained$tip.label, tree_retained$node.label)
alleles <- alleles[tree_retained$tip.label, ]

missing_names <- names(which(apply(apply(alleles, 1, is.na),2,any)))
tree_retained <- drop.tip(tree_retained, missing_names)

nodedates <- nodedates[c(tree_retained$tip.label, tree_retained$node.label)]
record <- record[,c(tree_retained$tip.label, tree_retained$node.label)]
alleles <- alleles[tree_retained$tip.label, ]

write.csv(alleles, "alleles.csv")
saveRDS(list(record=record, nodedates=nodedates), "timing_data.rds")
write.tree(tree_retained, "tree_retained.nwk")
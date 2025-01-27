library(stringr)
library(tidyverse)

extract_penA_snps <- function(penA_alleles, penA_desc)
{
    penA_SNPS <- unique(unlist(apply(penA_desc, 1, function(x) str_extract_all(x[2], "[A-Z]\\d+[A-Z]"))))
    penA_sites <- unique(unlist(sapply(penA_desc, function(x) str_extract_all(x, "[A-Z]\\d+"))))

    penA_motif <- t(apply(penA_alleles, 1,
    function(x) {
        ii <- str_extract(x, "\\d+")
        c(str_detect(penA_desc[ii, 2], penA_SNPS), ifelse(penA_desc[ii, 3] == "yes", TRUE, FALSE))
    }
    ))

    colnames(penA_motif) <- c(penA_SNPS, "Mosaic")
    penA_patterns <- matrix("", nrow=nrow(penA_motif), ncol=length(penA_sites)+1L)
    colnames(penA_patterns) <- c(penA_sites, "Mosaic")

    for (x in colnames(penA_motif)) {
        WT_type <- ""
        mut_type <- ""
        pat_name <- ""
        if (x!= "Mosaic")
        {
            pat_name <- str_extract(x, "[A-Z]\\d+")
            WT_type <- str_extract_all(x, "[A-Z]")[[1]][1]
            mut_type <- str_extract_all(x, "[A-Z]")[[1]][2]
        } else {
            WT_type <- "No"
            mut_type <- "Yes"
            pat_name <- "Mosaic"
        }
        penA_patterns[, pat_name] <- ifelse(penA_patterns[, pat_name]=="", WT_type, penA_patterns[, pat_name])
        penA_patterns[, pat_name] <- ifelse(penA_motif[, x], mut_type, penA_patterns[, pat_name])
    }
    rownames(penA_patterns) <- rownames(penA_motif)
    colnames(penA_patterns) <- paste0("penA_", colnames(penA_patterns))
    return(penA_patterns)
}

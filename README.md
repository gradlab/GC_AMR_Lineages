# Code necessary to reproduce statistical analyses featured in _Quantifying the impact of genetic determinants of antibiotic resistance on bacterial lineage dynamics_
## Layout
This repository contains all the code and data that is needed to reproduce the results featured in _Quantifying the impact of genetic determinants of antibiotic resistance on bacterial lineage dynamics_. The input data provided consists of all sample metadata, phylogenetic trees, and usage data that was used as an input for the various analyses featured in the accompanying manuscript.
* The top level directory contains several subdirectories as well as the file ```run_model.R```. This file contains the main analyses and uses inputs generated during lineage calling steps.
* The ```snakemake_pipeline``` directory contains the genome assembly and variant calling pipeline. You probably don't want to run this as it requires access to HPC. All the required outputs are already included within the data directory.
* The directories ```plots/``` and ```tables/``` contain plots and tables generated during the various steps that make up the analysis.
* The ```plotting/``` directory contains scripts used to generate specialised figures used within the manuscript.
* The ```models/``` directory contains files that make up the phylodynamic model, as well as generic plotting functions for visualising the model's outputs, along with preprocessing functions.
* The ```data/``` directory contains data used by the phylodynamic model, as well as lineage calling scripts.
  + ```data/metadata``` contains metadata and covariate files.
  + ```data/pastml``` contains raw ```pastml``` output files as well as a script for pre-processing pastml inputs.
  + ```data/split_lineages.R``` is the main script used to assign lineages and generate input files for the phylodynamic model.

## Usage
The execution order is as follows:
* Generate lineage data by navigating to ```data/``` and running ```split_lineages.R```.
* Run the phylodynamic model by running ```run_model.R``` from the top-level directory.

## Sequence Acession Numbers and Metadata
The accession numbers and metadata can be found under ```data/metadata/Ng-Combined-Metadata.txt```

## Dependencies
The following packages are required to execute the software in this repository and available on CRAN:
* ```tidyverse```
* ```cmdstanr```
* ```ape```
* ```stringi```
* ```glue```
* ```gt```
* ```ggplot2```
* ```gridExtra```
* ```viridis```
* ```ggnewscale```
* ```patchwork```
* ```ggridges```
* ```RColorBrewer```
* ```gtExtras```
* ```posterior```

The following packages are required and have to be installed from their respective ```github```repositories:
* ```phylodyn``` [https://github.com/mdkarcher/phylodyn](https://github.com/mdkarcher/phylodyn)

The following packages are required and available on Bioconductor:
* ```treeio```
* ```ggtree```

Here's the detailed session info including all package versions:
```
R version 4.4.0 (2024-04-24)
Platform: aarch64-apple-darwin20
Running under: macOS 15.6

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggtree_3.11.2      treeio_1.27.1      phylodyn_0.9.02    posterior_1.5.0   
 [5] gtExtras_0.5.0     RColorBrewer_1.1-3 ggridges_0.5.6     patchwork_1.2.0   
 [9] ggnewscale_0.4.10  viridis_0.6.5      viridisLite_0.4.2  gridExtra_2.3     
[13] gt_0.11.0          glue_1.7.0         stringi_1.8.3      ape_5.8           
[17] cmdstanr_0.7.1     lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1     
[21] dplyr_1.1.4        purrr_1.0.2        readr_2.1.5        tidyr_1.3.1       
[25] tibble_3.2.1       ggplot2_3.5.1      tidyverse_2.0.0   

loaded via a namespace (and not attached):
 [1] gtable_0.3.5         tensorA_0.36.2.1     xfun_0.43           
 [4] processx_3.8.4       lattice_0.22-6       paletteer_1.6.0     
 [7] tzdb_0.4.0           yulab.utils_0.1.4    vctrs_0.6.5         
[10] tools_4.4.0          ps_1.7.6             generics_0.1.3      
[13] parallel_4.4.0       fansi_1.0.6          pkgconfig_2.0.3     
[16] ggplotify_0.1.2      checkmate_2.3.1      distributional_0.4.0
[19] lifecycle_1.0.4      compiler_4.4.0       munsell_0.5.1       
[22] fontawesome_0.5.2    ggfun_0.1.4          htmltools_0.5.8.1   
[25] lazyeval_0.2.2       pillar_1.9.0         cachem_1.0.8        
[28] abind_1.4-5          nlme_3.1-164         aplot_0.2.2         
[31] tidyselect_1.2.1     digest_0.6.35        rematch2_2.1.2      
[34] fastmap_1.1.1        grid_4.4.0           colorspace_2.1-0    
[37] cli_3.6.2            magrittr_2.0.3       utf8_1.2.4          
[40] withr_3.0.0          scales_1.3.0         backports_1.4.1     
[43] timechange_0.3.0     hms_1.1.3            memoise_2.0.1       
[46] knitr_1.46           gridGraphics_0.5-1   rlang_1.1.3         
[49] Rcpp_1.0.12          tidytree_0.4.6       xml2_1.3.6          
[52] jsonlite_2.0.0       R6_2.5.1             fs_1.6.4            
```

# Code necessary to reproduce statistical analyses featured in _Quantifying the impact of genetic determinants of antibiotic resistance on bacterial lineage dynamics_
## Layout
This repository contains all the code and data that is needed to reproduce the results featured in _Quantifying the impact of genetic determinants of antibiotic resistance on bacterial lineage dynamics_. The input data provided consists of all sample metadata, phylogenetic trees, and usage data that was used as an input for the various analyses featured in the accompanying manuscript.
* The top level directory contains several subdirectories as well as the file ```run_model.R```. This file contains the main analyses and uses inputs generated during lineage calling steps.
* The ```snakemake_pipeline``` directory contains the genome assembly and variant calling pipeline
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

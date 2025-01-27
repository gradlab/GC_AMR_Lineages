#!/bin/bash
#SBATCH -p hsph
#SBATCH -c 20
#SBATCH -t 3-00:00
#SBATCH --mem-per-cpu=8G
#SBATCH -o treetime_pseudogenomes.out
#SBATCH -e treetime_pseudogenomes.err
#SBATCH --mail-type=END
#SBATCH --mail-user=dhelekal@hsph.harvard.edu 

treetime --aln combined.filtered_polymorphic_sites.fasta --tree tree_per_site.tre --dates dates_treetime.csv  --coalescent skyline --outdir treetime_out/
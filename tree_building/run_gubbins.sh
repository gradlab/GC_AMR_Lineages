#!/bin/bash
#SBATCH -p intermediate
#SBATCH -c 20
#SBATCH -t 14-00:00
#SBATCH --mem-per-cpu=8G
#SBATCH -o gubbins_pseudogenomes.out
#SBATCH -e gubbins_pseudogenomes.err
#SBATCH --mail-type=END
#SBATCH --mail-user=dhelekal@hsph.harvard.edu 

run_gubbins.py --first-tree-builder rapidnj --tree-builder iqtree --first-model JC --model GTR --threads 20 --iterations 10 combined.fasta

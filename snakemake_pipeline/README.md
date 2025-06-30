# *N. gonorrhoeae* assembly pipeline implemented in snakemake

## Pipeline Description

This pipeline performs QC, de novo assembly, mapping to a reference genome, annotation, and resistance-associated allele calling using *N. gonorrhoeae* genomic data
sequenced on the Illumina platform as input.

![Flow chart of pipeline](PipelineOverview.png)

## Installing snakemake
1. [Install conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
2. [Install snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

## Installing biopython
Follow instructions here to install biopython: https://biopython.org/wiki/Download

## Running the pipeline

1. Make a directory for your current project
2. Copy the contents of this directory (`Snakefile`, `start_snakemake.sh`, `conda_envs/`, `reference_sequences/`, `slurm_genomics_pipeline/`, `penA_alleles/`)
3. Edit the snakemake config (`slurm_genomics_pipeline/config.yaml`) with your conda installation path, specifically the line that starts with conda-prefix:
4. Edit the submit script (`start_snakemake.sh`) with your email address (you can also update the time/queue/etc. depending on how many isolates you are assembling)
5. Add executable permissions to `slurm-status.py` if needed (`chmod +x slurm_genomics_pipeline/slurm-status.py`)
6. Make sure your channel priorities are set to `flexible` rather than strict (`channel_priority: flexible` in your `.condarc` file). Note that if you don't have a .condarc file, you should first run `conda config` to create one.
7. Make a directory called `fastqs` with all your fastqs. Paired end files should be named using [sample_name]\_1.fastq.gz and [sample_name]\_2.fastq.gz naming convention.
8. Submit job using `sbatch start_snakemake.sh`

## Downloading fastqs

1. Search for the Bioproject number in the "Enter accession" text box in the upper right hand corner on the ENA website (ENA and NCBI are mirrored): https://www.ebi.ac.uk/ena/browser/home.
2. Scroll down to where it says "Download report: JSON TSV" and download the report in TSV format.
3. Edit the file to just be a list of URLs (replacing the ; between the paired URLs with a new line character). Usually you have to add ftp:// to the beginning of each line as well.
4. Download the URLs with `wget -i text_file_with_urls.txt`. May need to submit this as a batch job for large datasets.

Note that you can also use prefetch + fasterq-dump (default split-3) from sra-tools (https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump), but would additionally need to gzip the fastq files after downloading as that is what the pipeline expects.

## Troubleshooting tips
If the pipeline stops due to an error, you can try the following troubleshooting tips: 
1. Restart the pipeline.
2. Run `snakemake -n` to see which jobs still need to be run and check to see if error files exist for these jobs in `logs/slurm/*.err`. Also look at the progression in the flow diagram to determine which jobs may be the source of the error. 
3. Most common error is not enough time or memory for the variant rule in the Snakefile. Can increase the time or memory in the Snakefile by editing the multiplier (i.e. in line 148 increase the time multiplier from 60 to 240). The Snakefile tries 3 attempts for the time and memory allocation before throwing an error. 
4. Sometimes if the job needs to stop and restart, it may leave temporary files that cause it throw an error when it restarts. Delete the `mapping/*tmp*`, `variants(_16S,_23S)/*.fasta` files before restarting

## Saving relevant outputs
After finishing assemblies, uploading relevant outputs to `/n/grad_lab2/Lab/gonococcus/datasets` under a new folder with the name format `firstauthorlastname_yearofpublication_description` where the description is typically the country or smaller scale location, but can also be another short (1-2 word) description (i.e. historic, dgi, etc).

* Copy the folders `itol/`,`pseudogenomes/`,`qc/`,`resistance/`. 
* Copy the folder `mapping/` and rename it as `bams/`. Copy the folder `variants/` and rename is as `vcfs/`. 
* Copy `assemblies/filtered/*_contigs_filtered.fa` to a folder named `assemblies/`.
* Copy `annotations/*/*.gff` to a folder named `annotations/`.
* Move 23S and 16S VCFs to `all_23S_variant_calling` and `all_16S_variant_calling` respectively

Updated poppunk database (`gc_clusters/`) should be copied to `/n/grad_lab2/Lab/gonococcus/analyses/gc_clusters/` unless it's a special analysis.

## Updating metadata file
After finishing assemblies and uploading relevant outputs, update the metadata file using instructions [here](https://github.com/gradlab/lab-handbook/wiki/Cluster-Manual#updating-the-metadata-table) and [here](https://github.com/gradlab/gc_genomics/tree/master/metadata). Commit and push the changes, and in the commit message close the relevant issue. 

Another thing that's good to do is to spot check that the AMR phenotypes in the metadata are what we would expect from the mutations that were called by the pipeline. This is to ensure that there weren't sample mix-ups along the way! 

## Pipeline output

```
annotations/
    A directory for each sample with output from [prokka][https://github.com/tseemann/prokka#output-files]
blastdb/
    A blast database of contigs from assemblies
blast_results/
    Blast results from resistance associated genes
itol/
    ITOL annotation files from resistance associated SNPs (those called from pseudogenomes)
logs/
    Pipeline log files
mapping/
    BAM files for each sample mapped to NCCP11945
mapping_16S/
    BAM files for each sample mapped to the 16S rRNA sequence
mapping_23S/
    BAM files for each sample mapped to a single copy of 23S rRNA
pseudogenomes/
    NCCP11945 reference genome replaced with high quality variant calls or missing data for each sample
qc/
    bamqc/
        bamqc output for each sample
    fastqc/
        fastqc output for each sample
    multiqc_bamqc_data/
        summarized bamqc output
    multiqc_bamqc.html
        HTML representation of summarized bamqc output
    multiqc_bamqc_data/
        summarized fastqc output
    multiqc_bamqc.html
        HTML representation of summarized fastqc output
    qc_summary.txt
        Tab delimited file with the following columns:
            wgs_id: sample name
            assembly_length: total length of assembled contigs after they have been filtered removing short/low coverage contigs
            assembly_coverage: average coverage of filtered contigs as reported in contig headers
            contigs: total number of contigs
            genes: total number of annotated genes
            reference_coverage: average coverage of reads mapped to NCCP11945 reference
            reference_percentage_mapped: the percentage of total reads that mapped to NCCP11945 reference
            percent_missing: the percentage of sites in the reference genome that could not be confidently called in pseudogenome
     
resistance/
    TSVs describing the presence of previously described AMR associated alleles
variants/
    VCFs for each sample (NCCP11945)
variants_16S/
    VCFs for 16S rRNA locus
variants_23S/
    VCFs for 23S rRNA locus
```

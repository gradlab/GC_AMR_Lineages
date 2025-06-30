#!/usr/bin/env python

import sys
import argparse
import os
from datetime import datetime
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SearchIO
from io import StringIO
from collections import Counter

def get_args():
    parser = argparse.ArgumentParser(description='ID penA alleles in assemblies')
    parser.add_argument("assemblydir", help="Assembly directory")
    parser.add_argument("name", help="Data set name")
    parser.add_argument("blastdb", help="Path to blastdb")
    parser.add_argument("penA_descriptions", help="File describing penA alleles (i.e. mosaicism)")
    return parser.parse_args()

def get_sample_name(file_name):
    if file_name.endswith("_contigs_filtered.fa"):
        sample_name = file_name.replace("_contigs_filtered.fa", "")
    elif file_name.endswith("_with_plasmids.fna"):
        sample_name = file_name.replace("_with_plasmids.fna", "")
    elif file_name.endswith(".fna"):
        sample_name = file_name.replace(".fna", "")
    elif file_name.endswith(".fasta"):
        sample_name = file_name.replace(".fasta", "")
    elif file_name.endswith(".fa"):
        sample_name = file_name.replace(".fa", "")
    else:
        print(f"Warning: could not parse sample name from {file_name}")
        sample_name = file_name
    return(sample_name)

def run_blast(assembly, blastdb):
    blastn_cline = NcbiblastnCommandline(query = assembly,
                                         db = blastdb,
                                         evalue = 0.001,
                                         outfmt = 5)
    blast_out, blast_err = blastn_cline()
    blast_result = SearchIO.parse(StringIO(blast_out), "blast-xml")
    contigs_with_hits = [br for br in blast_result if len(br) > 0]
    if len(contigs_with_hits) > 1:
        print(f"Warning: multiple contigs encoding penA in {assembly}.")
        return("unknown")
    elif len(contigs_with_hits) == 0:
        return("unknown")
    top_hit = contigs_with_hits[0][0]
    penA_id = top_hit.id
    penA_length = top_hit.seq_len
    top_hsp = top_hit[0]
    ident = top_hsp.ident_num
    if penA_length == ident:
        return(penA_id.split(".")[0])
    else:
        return("unknown")

def summarize_allele_counts(penA_alleles, penA_descriptions):
    penA_counts = Counter(penA_alleles.values())
    mosaic = []
    with open(penA_descriptions, "r") as infile:
        for i,line in enumerate(infile):
            if i > 0:
                line = line.strip().split("\t")
                if line[2] == "yes":
                    mosaic.append("penA_" + line[0])
    penA_counts_sorted = sorted(penA_counts, key=penA_counts.get, reverse=True)
    top_mosaics = []
    for allele in penA_counts_sorted:
        if len(top_mosaics) == 10:
            break
        if allele in mosaic:
            top_mosaics.append(allele)
        else:
            continue
    return top_mosaics, mosaic

def write_itol(allele_dict, colors, top_mosaics, mosaics, name):
    with open(f"itol_{name}_penA.txt", "w") as itolfile:
        itolfile.write(f"DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tpenA\n")
        itolfile.write(f"COLOR\t{colors[1]}\nLEGEND_TITLE\tpenA\n")
        colors = colors[:len(top_mosaics)] + ["#2F4F4F", "#778899", "#ffffff"]
        color_string = "\t".join(colors)
        shape_string = "\t".join(["1"]*len(colors))
        allele_string = "\t".join(top_mosaics + ["Other Mosaic penA", "Nonmosaic penA", "Unknown"])
        itolfile.write(f"LEGEND_SHAPES\t{shape_string}\n")
        itolfile.write(f"LEGEND_COLORS\t{color_string}\n")
        itolfile.write(f"LEGEND_LABELS\t{allele_string}\n")
        itolfile.write("BORDER_WIDTH\t0.25\nBORDER_COLOR\t#CCCCCC\nDATA\n")
        for strain,allele in allele_dict.items():
            if allele in top_mosaics:
                color = colors[top_mosaics.index(allele)]
            elif allele in mosaics:
                color = "#2F4F4F"
            elif allele == "unknown":
                color = "#ffffff"
            else:
                color = "#778899"
            itolfile.write(f"{strain}\t{color}\t{allele}\n")

def write_summary_table(penA_alleles, name):
    with open(f"{name}_penA.txt", "w") as outfile:
        outfile.write("wgs_id\tpenA\n")
        for s in penA_alleles:
            outfile.write(f"{s}\t{penA_alleles[s]}\n")

def main():
    args = get_args()
    assemblies = os.listdir(args.assemblydir)
    penA_alleles = {}
    for assembly in assemblies:
        sample_name = get_sample_name(assembly)
        penA_alleles[sample_name] = run_blast(args.assemblydir + assembly, args.blastdb)

    top_mosaics, mosaics = summarize_allele_counts(penA_alleles, args.penA_descriptions)
    write_itol(penA_alleles, ["#beab3e","#b94973","#64aa53","#5c3788","#45c097","#6c81d9","#8e863a","#c26abb","#c36c34","#b74943"],top_mosaics, mosaics, args.name)
    write_summary_table(penA_alleles, args.name)

if __name__ == "__main__":
    main()

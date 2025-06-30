#!/usr/bin/env python

import sys
import argparse
import re
from Bio import SearchIO
from collections import defaultdict
from datetime import datetime

def get_args():
    parser = argparse.ArgumentParser(description="Update missing data in resistance allele table from de novo assemblies. Must be run after snakemake pipeline.")
    parser.add_argument("resistance_tsv", help="output of resistance_alleles_pseudogenomes.py")
    return parser.parse_args()

def read_alleles(res_alleles_file):
    res_loci_dict = defaultdict(dict)
    with open(res_alleles_file, "r") as infile:
        for i,line in enumerate(infile):
            if i == 0:
                line = line.strip().split('\t')
                loci = line[1:]
            else:
                line = line.strip().split('\t')
                for j,allele in enumerate(line[1:]):
                    res_loci_dict[line[0]][loci[j]] = allele
    return res_loci_dict,loci

def penA_alleles(res_loci_dict):
    missing_ids = []
    for wgs_id, allele_dict in res_loci_dict.items():
        penA_alleles = [allele_dict[locus] for locus in ["PBP2_501", "PBP2_512", "PBP2_516", "PBP2_542", "PBP2_551"]]
        if "NA" in penA_alleles:
            missing_ids.append(wgs_id)
    blast_output = SearchIO.read("blast_results/penA_gc_blast.xml", "blast-xml")
    for hit in blast_output:
        wgs_id = re.split('_\d+_', hit.id)[0]
        if wgs_id in missing_ids:
            penA_peptide = hit[0].hit.seq.ungap('-').translate()
            # check for insertion
            if "DDTHV" in str(penA_peptide):
                res_loci_dict[wgs_id]["PBP2_501"] = penA_peptide[501]
                res_loci_dict[wgs_id]["PBP2_512"] = penA_peptide[512]
                res_loci_dict[wgs_id]["PBP2_516"] = penA_peptide[516]
                res_loci_dict[wgs_id]["PBP2_542"] = penA_peptide[542]
                res_loci_dict[wgs_id]["PBP2_551"] = penA_peptide[551]
            else:
                res_loci_dict[wgs_id]["PBP2_501"] = penA_peptide[500]
                res_loci_dict[wgs_id]["PBP2_512"] = penA_peptide[511]
                res_loci_dict[wgs_id]["PBP2_516"] = penA_peptide[515]
                res_loci_dict[wgs_id]["PBP2_542"] = penA_peptide[541]
                res_loci_dict[wgs_id]["PBP2_551"] = penA_peptide[550]
    return res_loci_dict


def mtrD_alleles(res_loci_dict):
    missing_ids = []
    for wgs_id, allele_dict in res_loci_dict.items():
        mtrD_alleles = [allele_dict[locus] for locus in ["MtrD_823"]]
        if "NA" in mtrD_alleles:
            missing_ids.append(wgs_id)
    blast_output = SearchIO.read("blast_results/mtrD_gc_blast.xml", "blast-xml")
    for hit in blast_output:
        wgs_id = re.split('_\d+_', hit.id)[0]
        if wgs_id in missing_ids:
            mtrD_peptide = hit[0].hit.seq.ungap('-').translate()
            res_loci_dict[wgs_id]["MtrD_823"] = mtrD_peptide[822]
    return res_loci_dict


args = get_args()
res_loci_dict,loci = read_alleles(args.resistance_tsv)
res_loci_dict = penA_alleles(res_loci_dict)
res_loci_dict = mtrD_alleles(res_loci_dict)

with open('{0}_gc_resistance_alleles.tsv'.format(datetime.strftime(datetime.now(), '%Y-%m-%d')), "w") as outfile:
    outfile.write("wgs_id\t" + "\t".join(loci) + "\n")
    for s in res_loci_dict:
        a = [res_loci_dict[s][l] for l in loci]
        outfile.write("{0}\t{1}\n".format(s, "\t".join(a)))

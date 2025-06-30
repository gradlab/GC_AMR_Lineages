#!/usr/bin/env python

import sys
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="Filter sites that have an alteranative genotype call (1/1) but no alternate allele")
    parser.add_argument("vcf", help="VCF to filter")
    parser.add_argument("--replace-name", help="Replace SAMPLE in header with sample name",
                        action='store_true')
    return parser.parse_args()

def filter_vcf(vcf_file, replace):
    filter_vcf_file = vcf_file.split(".")[0] + "_filter.vcf"
    sample = vcf_file.split(".")[0].strip("_pilon").split("/")[-1]
    with open(vcf_file, "r") as infile:
        with open(filter_vcf_file, "w") as outfile:
            for line in infile:
                if line[0:2] == "##":
                    outfile.write(line)
                elif line[0] == "#":
                    if replace:
                        line = line.replace("SAMPLE", sample)
                    outfile.write(line)
                else:
                    record = line.strip().split()
                    ALT = record[4]
                    GENO = record[9]
                    if ALT != "." or GENO == "0/0":
                        outfile.write(line)

args = get_args()
filter_vcf(args.vcf, args.replace_name)

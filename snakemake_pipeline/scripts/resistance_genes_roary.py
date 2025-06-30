#!/usr/bin/env python

import sys
import argparse
from datetime import datetime

def get_args():
    parser = argparse.ArgumentParser(description="Create resistance gene matrix and itol files")
    parser.add_argument("csv", help="Gene presence/absence CSV from roary")
    return parser.parse_args()

color_palettes = {"bugn_2":["#b2e2e2", "#238b45"],
                  "bupu_2":["#b3cde3", "#88419d"],
                  "gnbu_2":["#bae4bc", "#2b8cbe"],
                  "orrd_2":["#fdcc8a", "#d7301f"],
                  "pubu_2":["#bdc9e1", "#0570b0"],
                  "pubugn_2":["#bdc9e1", "#02818a"],
                  "purd_2":["#d7b5d8", "#ce1256"],
                  "rdpu_2":["#fbb4b9", "#ae017e"],
                  "ylgn_2":["#c2e699", "#238443"],
                  "ylgnbu_2":["#a1dab4", "#225ea8"],
                  "ylorbr_2":["#fed98e", "#cc4c02"],
                  "ylorrd_2":["#fecc5c", "#e31a1c"],
                  "bugn_3":["#ccece6", "#66c2a4", "#006d2c"],
                  "bupu_3":["#bfd3e6", "#8c96c6", "#810f7c"],
                  "gnbu_3":["#ccebc5", "#7bccc4", "#0868ac"],
                  "orrd_3":["#fdd49e", "#fc8d59", "#b30000"],
                  "pubu_3":["#d0d1e6", "#74a9cf", "#045a8d"],
                  "purd_3":["#d4b9da", "#df65b0", "#980043"],
                  "rdpu_3":["#fcc5c0", "#f768a1", "#7a0177"],
                  "ylgn_3":["#d9f0a3", "#78c679", "#006837"],
                  "ylgnbu_6":["#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84"],
                  "pubugn_6":["#d0d1e6", "#a6bddb", "#67a9cf", "#3690c0", "#02818a", "#016450"]}

def read_csv(roary_csv):
    res_genes = ["bla", "tetA_1", "tetA_2", "tet(M)", "ermC\'", "traC"]
    res_alleles = {}
    with open(roary_csv, "r") as infile:
        for i,line in enumerate(infile):
            if i == 0:
                line = line.strip().split(",")
                strains = [x.strip('\"') for x in line[14:]]
            else:
                line = line.strip().split(",")
                gene = line[0].strip('\"')
                if gene in res_genes:
                    alleles = [x.strip('\"') for x in line[14:]]
                    res_alleles[gene] = alleles
    return(strains, res_alleles)

def write_itol(locus_name, color):
    presence = []
    for a in res_alleles[locus_name]:
        if a == "":
            presence.append("0")
        else:
            presence.append("1")
    locus_dict = dict(zip(strains,presence))
    with open("itol_{0}.txt".format(locus_name), "w") as itol_file:
        itol_file.write("DATASET_COLORSTRIP\n\n")
        itol_file.write("SEPARATOR TAB\n\n")
        itol_file.write("DATASET_LABEL\t{0}\nCOLOR\t{1}\n\n".format(locus_name, color[-1]))
        itol_file.write("LEGEND_TITLE\t{0}\nLEGEND_SHAPES\t1\t1\nLEGEND_COLORS\t{1}\nLEGEND_LABELS\tAbsent\tPresent\n\n".format(locus_name, "\t".join(color)))
        itol_file.write("BORDER_WIDTH\t0.25\nBORDER_COLOR\t#CCCCCC\n\n")
        itol_file.write("DATA\n")
        for strain,allele in locus_dict.items():
            if allele == "0":
                c = color[0]
            else:
                c = color[1]
            itol_file.write("{0}\t{1}\t{2}\n".format(strain, c, allele))

    return(locus_dict)


args = get_args()


strains, res_alleles = read_csv(args.csv)

res_loci_dict = {}
res_loci_dict["bla"] = write_itol("bla", color_palettes["bugn_2"])
res_loci_dict["tetA_1"] = write_itol("tetA_1", color_palettes["gnbu_2"])
res_loci_dict["tetA_2"] = write_itol("tetA_2", color_palettes["pubu_2"])
res_loci_dict["tetM"] = write_itol("tet(M)", color_palettes["purd_2"])
res_loci_dict["ermC\'"] = write_itol("ermC\'", color_palettes["ylgn_2"])
res_loci_dict["GGI"] = write_itol("traC", color_palettes["rdpu_2"])

with open('{0}_gc_resistance_genes.tsv'.format(datetime.strftime(datetime.now(), '%Y-%m-%d')), "w") as outfile:
    loci = list(res_loci_dict.keys())
    outfile.write("Accession\t" + "\t".join(loci) + "\n")
    for s in strains:
        a = [res_loci_dict[l][s] for l in loci]
        outfile.write("{0}\t{1}\n".format(s, "\t".join(a)))

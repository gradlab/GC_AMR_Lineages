#!/usr/bin/env python

from Bio import SeqIO
from collections import defaultdict
import re

def call_LOF(seq_record, peptide_length):
    peptide = str(seq_record.seq.ungap('-').translate()).split('*')[0]
    if len(peptide) < 0.8*peptide_length:
        return "LOF"
    else:
        return "full_length"

mtr_mosaic_dictionary = defaultdict(lambda: defaultdict(lambda: "NA"))

genes = ["mtrR", "mtrC", "mtrD", "mtrE"]
peptide_lengths = {"mtrR":211, "mtrC":413, "mtrD":1068, "mtrE":468}

for g in genes:
    for rec in SeqIO.parse(f"blast_results/{g}_gc.fa", "fasta"):
        rec.id = re.split('_\d+_length', rec.id)[0]
        lof = call_LOF(rec, peptide_lengths[g])
        if lof == "LOF":
            mtr_mosaic_dictionary[rec.id][g] = "LOF"
        else:
            mtr_mosaic_dictionary[rec.id][g] = "GC_allele"
    for rec in SeqIO.parse(f"blast_results/{g}_gc_mosaics.fa", "fasta"):
        rec.id = re.split('_\d+_length', rec.id)[0]
        lof = call_LOF(rec, peptide_lengths[g])
        if lof == "LOF":
            mtr_mosaic_dictionary[rec.id][g] = "LOF"
        else:
            mtr_mosaic_dictionary[rec.id][g] = "mosaic"
    for rec in SeqIO.parse(f"blast_results/{g}_gc_lengthOutliers.fa", "fasta"):
        rec.id = re.split('_\d+_length', rec.id)[0]
        lof = call_LOF(rec, peptide_lengths[g])
        if lof == "LOF":
            mtr_mosaic_dictionary[rec.id][g] = "LOF"
        else:
            mtr_mosaic_dictionary[rec.id][g] = "lengthOutlier"

mtr_promoter_dictionary = {}

for i,rec in enumerate(SeqIO.parse("blast_results/mtrPromoter_gc.fa", "fasta")):
    rec.seq = rec.seq.ungap('-')
    rec.id = re.split('_\d+_length', rec.id)[0]
    if "TAAAAAAG" in rec.seq and "CTTTTTA" in rec.seq and "AACCT" in rec.seq:
        mtr_promoter_dictionary[rec.id] = "WT"
    elif "AATCT" in rec.seq:
        mtr_promoter_dictionary[rec.id] = "Tpos120"
    elif "CTTTTA" in rec.seq or "CTTTTC" in rec.seq:
        mtr_promoter_dictionary[rec.id] = "Adel_13bp_invertedRepeat"
    elif "TAAAAAAAAG" in rec.seq:
        mtr_promoter_dictionary[rec.id] = "TTins_13bp_invertedRepeat"
    elif "TAAAAAAAG" in rec.seq:
        mtr_promoter_dictionary[rec.id] = "Tins_13bp_invertedRepeat"
    elif "TAAAAAG" in rec.seq:
        mtr_promoter_dictionary[rec.id] = "Tdel_13bp_invertedRepeat"
    elif "CTTTGTA" in rec.seq:
        mtr_promoter_dictionary[rec.id] = "AtoC_13bp_invertedRepeat"
    else:
        print(rec.id)

with open("resistance/mtr_alleles.tsv", "w") as outfile:
    outfile.write("wgs_id\tmtrPromoter\tmtrR\tmtrC\tmtrD\tmtrE\n")
    for sample in mtr_mosaic_dictionary:
        mtr_alleles = mtr_mosaic_dictionary[sample]
        if "mtrR" not in mtr_alleles:
            mtr_alleles["mtrR"] = "blast_hit_missing_possible_mosaic"
        alleles = [mtr_alleles[mtr_region] for mtr_region in genes]
        allele_string = "\t".join(alleles)
        if sample not in mtr_promoter_dictionary:
            mtr_promoter_dictionary[sample] = "divergent_promoter"
        outfile.write(f"{sample}\t{mtr_promoter_dictionary[sample]}\t{allele_string}\n")

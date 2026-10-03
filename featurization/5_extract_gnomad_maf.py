#!/usr/bin/env python
"""
gnomAD frequencies of the prepared variants (the gnomad_MAF feature), defined as in training: INFO/AF of
the gnomAD v3.1 genomes record with the same two alleles, in either orientation. This is the frequency of
the record's ALT allele, not a minor-allele frequency; a record without AF counts as 0. Variants not in
gnomAD get no value, and the scoring step sets them to the median over the scored variants, as training did.
The training values are those of gnomAD v3.1 genomes (later releases differ for some variants).

usage:   python 5_extract_gnomad_maf.py OUT_DIR --vcf 'DIR/gnomad.genomes.v3.1.sites.{chrom}.vcf.bgz'
         --vcf: the per-chromosome gnomAD v3.1 genomes sites VCFs, with their .tbi indexes (README); a URL
         works too, as tabix reads only the variants' positions (the indexes are saved in OUT_DIR/gnomad_index)
input:   OUT_DIR/variants.tsv (1_prepare_variants.py)
output:  OUT_DIR/gnomad_MAF.tsv   variant_id, gnomad_MAF
"""
import argparse
import os
import sys

import pandas as pd
import pysam

ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
ap.add_argument("out_dir")
ap.add_argument("--vcf", required=True, help="per-chromosome gnomAD v3.1 genomes VCFs, with {chrom} for the chromosome (e.g. chr7)")
a = ap.parse_args()
if "{chrom}" not in a.vcf:
    sys.exit("--vcf needs {chrom}, e.g. 'gnomad/gnomad.genomes.v3.1.sites.{chrom}.vcf.bgz'")


def af_text(rec):
    """INFO/AF as written in the VCF ('.' or absent -> 0), parsed from the text as training did."""
    info = str(rec).split("\t", 8)[7]
    value = next((kv[3:] for kv in info.split(";") if kv.startswith("AF=")), ".")
    return 0.0 if value == "." else float(value.split(",")[0])


out_dir = os.path.abspath(a.out_dir)
variants = pd.read_csv(f"{out_dir}/variants.tsv", sep="\t")
if "://" in a.vcf:                       # htslib saves the index of a remote VCF in the working directory
    os.makedirs(f"{out_dir}/gnomad_index", exist_ok=True)
    os.chdir(f"{out_dir}/gnomad_index")
rows = []
for chrom, v in variants.groupby("CHR"):
    vcf = pysam.VariantFile(a.vcf.format(chrom=chrom))
    for r in v.itertuples():
        same, swapped = None, None
        for rec in vcf.fetch(chrom, r.BP - 1, r.BP):
            if rec.pos != r.BP or not rec.alts:
                continue
            af = af_text(rec)
            if same is None and (rec.ref, rec.alts[0]) == (r.REF, r.ALT):
                same = af
            elif swapped is None and (rec.ref, rec.alts[0]) == (r.ALT, r.REF):
                swapped = af
        maf = same if same is not None else swapped
        if maf is not None:
            rows.append((r.variant_id, maf))
out = pd.DataFrame(rows, columns=["variant_id", "gnomad_MAF"])
out.to_csv(f"{out_dir}/gnomad_MAF.tsv", sep="\t", index=False)
print(f"wrote {out_dir}/gnomad_MAF.tsv: {len(out)} of {len(variants)} variants in gnomAD")

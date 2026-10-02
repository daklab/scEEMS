#!/usr/bin/env python
"""
MAGMA gene analysis for one chromosome, against the AD GWAS of one population with the ADSP reference
panel of that population for LD. Three kinds of gene annotation are tested:

  emagma      the scEEMS annotations of make_magma_files.py, for one cell type: "prediction" (variants
              the model predicts as eQTLs, each assigned to one gene) and "pip" (fine-mapped variants,
              PIP > 0.10); European GWAS (Bellenguez) or ADGC AFR, AMR or EAS
  knn         the TSS-distance control of make_magma_files_knn.py, for one cell type; European GWAS
  positional  the gene-body annotation of make_positional_annotation.py; European GWAS

Gene p-values use MAGMA's default SNP-wise mean model; --gene-settings adap-permp=10000 adds adaptive
permutation p-values (column PERMP, which varies between runs).

usage:   python run_magma.py CHR emagma COHORT [POP]     POP: EUR (default), AFR, AMR or EAS
         python run_magma.py CHR knn COHORT
         python run_magma.py CHR positional
outputs: {magma_dir}/{cohort}/MAGMA_output{,_ADGC_POP}/{cohort}_chr{CHR}_MAGMA_{prediction,pip}.genes.out
         {magma_dir}/knn/{cohort}/MAGMA_output/{cohort}_chr{CHR}_MAGMA_{prediction,pip}.genes.out
         {magma_dir}/positional/MAGMA_output/positional_chr{CHR}_MAGMA.genes.out
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import cohort_name, path

chrom = int(sys.argv[1])
kind = sys.argv[2]
assert kind in ("emagma", "knn", "positional"), "annotation must be emagma, knn or positional"
pop = sys.argv[4] if kind == "emagma" and len(sys.argv) > 4 else "EUR"
assert pop in ("EUR", "AFR", "AMR", "EAS"), "POP must be EUR, AFR, AMR or EAS"

MAGMA_DIR = path("magma_dir")
bfile = f"{path('plink_ref_dir')}/{pop}/plink/ADSP_{pop}_chr{chrom}"
sumstats = f"{MAGMA_DIR}/sumstats/" + ("bellenguez_MAGMA_sumstats.txt" if pop == "EUR"
                                      else f"ADGC_{pop}_MAGMA_sumstats.txt")


def magma(annot, out):
    cmd = [path("magma_binary"), "--bfile", bfile, "--gene-annot", annot, "--pval", sumstats, "ncol=N",
           "--gene-settings", "adap-permp=10000", "--out", out]
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


if kind == "positional":
    out_dir = f"{MAGMA_DIR}/positional/MAGMA_output"
    os.makedirs(out_dir, exist_ok=True)
    magma(f"{MAGMA_DIR}/positional/bellenguez_MAGMA.genes.annot", f"{out_dir}/positional_chr{chrom}_MAGMA")
else:
    cohort = cohort_name(sys.argv[3])
    root = f"{MAGMA_DIR}/{cohort}" if kind == "emagma" else f"{MAGMA_DIR}/knn/{cohort}"
    out_dir = f"{root}/MAGMA_output" + ("" if pop == "EUR" else f"_ADGC_{pop}")
    os.makedirs(out_dir, exist_ok=True)
    for arm in ("prediction", "pip"):
        magma(f"{root}/{arm}/{cohort}_chr{chrom}_MAGMA.genes.annot", f"{out_dir}/{cohort}_chr{chrom}_MAGMA_{arm}")

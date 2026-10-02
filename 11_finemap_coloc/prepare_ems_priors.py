#!/usr/bin/env python
"""
Per-gene EMS scores for the EMS fine-mapping prior (finemap_and_coloc.R, prior EMS), from the
expression modifier scores of Wang et al. (2021) for one GTEx tissue: whole blood for microglia and
brain frontal cortex (BA9) for the other cell types by default.

For each gene of the cell type's region list that has step 6 predictions, the EMS scores of the gene's
cis variants (the variants of its prediction table) are written; variants without a score get p_random
in finemap_and_coloc.R, where p_random = mean(ems / ems_normalized) over the tissue's table, as in the
EMS paper.

usage:   python prepare_ems_priors.py CELL [--tissue TISSUE]
input:   {ems_dir}/ems_top_{TISSUE}.tsv.bgz (columns v, g, ems, ems_normalized; from the EMS release)
output:  {finemap_dir}/{CELL}/ems_vector/{gene_id}.tsv (variant_id, ems), p_random.txt, tissue.txt,
         coverage.tsv
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import cohort_name, path

TISSUE = {"Mic": "Whole_Blood", "Ast": "Brain_Frontal_Cortex_BA9", "Exc": "Brain_Frontal_Cortex_BA9",
          "Inh": "Brain_Frontal_Cortex_BA9", "Oli": "Brain_Frontal_Cortex_BA9", "OPC": "Brain_Frontal_Cortex_BA9"}

ap = argparse.ArgumentParser()
ap.add_argument("cell", choices=list(TISSUE))
ap.add_argument("--tissue", default=None, help="GTEx tissue of the EMS table (default: see TISSUE)")
args = ap.parse_args()
tissue = args.tissue or TISSUE[args.cell]
out_dir = os.path.join(path("finemap_dir"), args.cell, "ems_vector")
os.makedirs(out_dir, exist_ok=True)


def translate_variant_id(v):
    """chr10_100006504_T_C_b38 -> chr10:100006504:T:C"""
    v = v.replace("_b38", "")
    parts = v.split("_")
    return ":".join(parts) if len(parts) == 4 else v


ems = pd.read_csv(os.path.join(path("ems_dir"), f"ems_top_{tissue}.tsv.bgz"), sep="\t", compression="gzip")
ems["variant_id"] = ems["v"].map(translate_variant_id)
ems["gene_id"] = ems["g"].str.split(".", n=1).str[0]                # drop the Ensembl version
ems = ems[["variant_id", "gene_id", "ems", "ems_normalized"]].drop_duplicates(subset=["variant_id", "gene_id"])
p_random = float((ems["ems"] / ems["ems_normalized"]).mean())
with open(os.path.join(out_dir, "p_random.txt"), "w") as fh:
    fh.write(f"{p_random:.6e}\n")
with open(os.path.join(out_dir, "tissue.txt"), "w") as fh:
    fh.write(f"{tissue}\n")
print(f"[{args.cell}] {tissue}: {len(ems):,} EMS rows, p_random = {p_random:.6e}", flush=True)

regions = pd.read_csv(os.path.join(path("eqtl_data_dir"), args.cell, "phenotype",
                                   f"snuc_pseudo_bulk.{args.cell}.mega.normalized.log2cpm.region_list.txt"),
                      sep="\t")
pred_dir = os.path.join(path("predictions_dir", cohort=cohort_name(args.cell)), "weighted_full")
ems_by_gene = dict(tuple(ems.groupby("gene_id")))

coverage, n_skipped = [], 0
for gene_id in regions["ID"].unique():
    pred = os.path.join(pred_dir, f"{gene_id}_predictions.tsv")
    if not os.path.exists(pred):
        n_skipped += 1
        continue
    window = pd.read_csv(pred, sep="\t", usecols=["variant_id"]).drop_duplicates()
    gene_ems = ems_by_gene.get(gene_id)
    if gene_ems is None:
        merged = window.assign(ems=float("nan"))
    else:
        merged = window.merge(gene_ems[["variant_id", "ems"]], on="variant_id", how="left")
    n_with = int(merged["ems"].notna().sum())
    merged.loc[merged["ems"].notna(), ["variant_id", "ems"]].to_csv(
        os.path.join(out_dir, f"{gene_id}.tsv"), sep="\t", index=False)
    coverage.append({"gene_id": gene_id, "n_in_window": len(window), "n_with_ems": n_with,
                     "coverage": n_with / len(window) if len(window) else 0.0})

cov = pd.DataFrame(coverage)
cov.to_csv(os.path.join(out_dir, "coverage.tsv"), sep="\t", index=False)
print(f"[{args.cell}] {len(cov):,} genes written ({n_skipped:,} without predictions), "
      f"mean coverage {cov['coverage'].mean():.3f}")

#!/usr/bin/env python
"""
Feature table of the prepared variant-gene pairs, computed as steps 2 (annotate_variants.py) and 3
(create_gene_datasets.py) computed the training data, from the scores of steps 2-3 of this folder and the
files of the data release (featurization/):
  Enformer      diff_32 of the 4,675 non-CAGE tracks (2_score_enformer.py)
  cell type     whether the variant position lies in each of the 40 cell-type ATAC/promoter/enhancer
                annotations (celltype_annotations.bed.gz), and per cell type the product of its
                promoter/enhancer/ATAC union (500 bp) annotation with its TF score (largest Enformer ChIP
                diff_32 of the cell type's TFs, truncated to an integer)
  TF            largest, smallest and largest absolute Enformer ChIP diff_32 of each cell type's TFs
  ABC           per cell type, the largest ABC score of the gene's enhancer links within 1,024 bp of the
                variant (abc_scores.tsv.gz); 0 if there is none
  ChromBPNet    per cell type, score and assay, the largest (>= 0), smallest (<= 0) and largest absolute
                score over the variant's peaks (3_score_chrombpnet.py); 0 for a variant near no peak
  distance      gene_TSS - position, its absolute value and its log (genes.tsv)
  baseline      0, their value throughout the training data
Variant type, gene constraint, gnomAD frequency and GPN-STAR are added at scoring (7_score_variants.py), as
in training.

usage:   python 6_build_features.py OUT_DIR
input:   OUT_DIR/variants.tsv, pairs.tsv, enformer.parquet, chrombpnet.tsv
output:  OUT_DIR/features.parquet   one row per variant-gene pair
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
import featurize as F
from config import path

CELLS = ["microglia", "astrocyte", "oligodendrocyte", "neuron"]
ASSAYS = ["ATAC", "H3K27ac", "H3K4me3"]
SCORES = ["log_counts_diff_chrombpnet", "log_probs_diff_abs_sum_chrombpnet", "probs_jsd_diff_chrombpnet"]
ABC_WINDOW = 1024

# TFs of each cell type: for the cell-type products (step 2) ...
TF_STEP2 = {
    "microglia": ["ETS2", "FLI1", "SPI1", "IRF8", "PRDM1", "CEBPA", "CEBPB", "CEBPD", "CEBPE", "KLF11", "KLF2",
                  "TFEC", "BACH1", "BATF", "BATF3", "RUNX1", "RUNX2", "RUNX3", "LYL1", "TAL1", "ELK3", "MEF2C",
                  "CREB3L2", "ATF4", "MAFB", "SALL1", "eGFP-SALL1", "SMAD5", "MEF2A", "MEF2B", "SMAD2", "USF1",
                  "STAT3", "eGFP-MAFG", "NFYB", "NRF1", "CREB1", "IRF1", "SOX9"],
    "neuron": ["ASCL1", "ASCL5", "BHLHA15", "BHLHE22", "NEUROD1", "NEUROD2", "NEUROD6", "TWIST2", "EGR4", "SP8",
               "SP9", "MEIS3", "PKNOX2", "MSC", "KLF5", "KLF8", "TBR1", "HLF"],
    "astrocyte": ["MEIS2", "NFIB", "TGIF1", "EMX2", "LHX2", "RFX2", "RFX4", "RORA", "RORB", "SOX1", "SOX2",
                  "SOX21", "SOX9", "NFATC4", "SOX5", "POU3F2", "POU3F3", "POU3F4", "FOXG1", "FOXO1", "SP5"],
    "oligodendrocyte": ["SOX10", "SOX13", "SOX21", "SOX3", "SOX6", "SOX8", "NHLH2", "NFIX", "SP7", "NFE2",
                        "E2F1", "CREB5", "POU3F3", "MYCN"],
}
# ... and for the TF scores (step 3), which add further TFs
TF_STEP3 = {
    "microglia": TF_STEP2["microglia"] + ["SP1", "MAX", "ETS1", "HINFP", "ZNF416", "MEF2D", "NFIA"],
    "neuron": TF_STEP2["neuron"] + ["MEF2B", "SP1", "NRF1", "GATA3", "MEF2C", "EGR2", "GFI1B", "TEAD1"],
    "astrocyte": TF_STEP2["astrocyte"] + ["RFX1", "MYB", "NFIC", "FOXK2", "ATF3", "SP1", "ATF1", "ZNF416", "MEF2D"],
    "oligodendrocyte": TF_STEP2["oligodendrocyte"] + ["CTCF", "SP1", "RFX1", "FOSL2", "MYB", "E2F7", "ATF7",
                                                      "RUNX1", "USF1", "YY2", "GFI1B", "FOXJ2", "ETV1"],
}

ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
ap.add_argument("out_dir")
a = ap.parse_args()
out_dir = a.out_dir

# ---- Enformer: variant-level, without the CAGE tracks
targets = pd.read_csv(path("targets_file"), sep="\t")
cage = {f"diff_32_{i}" for i in targets.loc[targets["description"].str.contains("CAGE"), "identifier"]}
chip = targets[targets["sum_stat"].str.contains("mean") & targets["description"].str.contains("CHIP:")].copy()
chip["TF"] = chip["description"].str.split(":", expand=True)[1]


def tf_columns(tfs):
    return [f"diff_32_{i}" for i in chip.loc[chip["TF"].isin(tfs), "identifier"]]


variants = pd.read_csv(f"{out_dir}/variants.tsv", sep="\t")
enformer = pd.read_parquet(f"{out_dir}/enformer.parquet").drop(columns=["SNP"])
enformer = enformer[[c for c in enformer.columns if c not in cage]]
df = variants.merge(enformer, on=["CHR", "BP", "REF", "ALT"], how="inner")
if len(df) < len(variants):
    print(f"warning: {len(variants) - len(df)} variants have no Enformer scores and are left out", flush=True)

# ---- cell-type annotations: the variant position inside an interval ([BP-1, BP) overlaps [start, end))
new = {}
annotations = pd.read_csv(path("celltype_annotations_file"), sep="\t")
for name, ann in annotations.groupby("annotation", sort=False):
    hit = np.zeros(len(df), dtype=int)
    for chrom, iv in ann.groupby("chrom"):
        iv = iv.sort_values("start")
        rows = np.flatnonzero(df["CHR"].to_numpy() == chrom)
        bp = df["BP"].to_numpy()[rows]
        i = np.searchsorted(iv["start"].to_numpy(), bp - 1, side="right") - 1
        hit[rows] = ((i >= 0) & (iv["end"].to_numpy()[np.maximum(i, 0)] >= bp)).astype(int)
    new[name] = hit
for c in F.load_columns_dict()["baseline"]:
    new[c] = np.zeros(len(df), dtype=int)
for cell in CELLS:
    new[f"none-none-{cell}"] = df[tf_columns(TF_STEP2[cell])].max(axis=1).astype("int").to_numpy()
    new[f"{cell}-enhancer_promoter_union_atac_500-{cell}"] = \
        (new[f"{cell}_enhancer_promoter_union_atac_500"] * new[f"none-none-{cell}"]).astype(int)
    cols = tf_columns(TF_STEP3[cell])
    new[f"max_{cell}_TF_score"] = df[cols].max(axis=1).to_numpy()
    new[f"min_{cell}_TF_score"] = df[cols].min(axis=1).to_numpy()
    new[f"abs_max_{cell}_TF_score"] = df[cols].abs().max(axis=1).clip(lower=0).to_numpy()

# ---- ChromBPNet: per variant, summaries over its peaks
cbp = pd.read_csv(f"{out_dir}/chrombpnet.tsv", sep="\t")
for cell in CELLS:
    for score in SCORES:
        for assay in ASSAYS:
            s = cbp.loc[(cbp["cell_type"] == cell) & (cbp["assay"] == assay)].dropna(subset=[score])
            g = s.groupby("variant_id")[score]
            base = f"{cell}-{score}"
            new[f"{base}_max_{assay}"] = df["variant_id"].map(g.max().clip(lower=0)).fillna(0).to_numpy()
            new[f"{base}_min_{assay}"] = df["variant_id"].map(g.min().clip(upper=0)).fillna(0).to_numpy()
            new[f"{base}_abs_max_{assay}"] = df["variant_id"].map(s[score].abs().groupby(s["variant_id"]).max()).fillna(0).to_numpy()
df = pd.concat([df, pd.DataFrame(new, index=df.index)], axis=1)

# ---- variant-gene pairs: distance to the TSS and ABC
pairs = pd.read_csv(f"{out_dir}/pairs.tsv", sep="\t").merge(
    pd.read_csv(path("genes_file"), sep="\t")[["gene_id", "gene_TSS"]], on="gene_id", how="left")
df = pairs.merge(df, on="variant_id", how="inner")
df["distance_TSS"] = df["gene_TSS"].astype(float) - df["BP"]          # float, as in step 3 (mean TSS)
df["abs_distance_TSS"] = df["distance_TSS"].abs()
with np.errstate(divide="ignore"):
    df["abs_distance_TSS_log"] = np.log(df["abs_distance_TSS"])

abc = pd.read_csv(path("abc_scores_file"), sep="\t")
abc = abc[abc["gene_id"].isin(set(df["gene_id"]))]
groups = {k: g for k, g in abc.groupby(["cell_type", "gene_id", "chrom"])}
for cell in CELLS:
    best = []
    for r in df[["gene_id", "CHR", "BP"]].itertuples(index=False):
        g = groups.get((cell, r.gene_id, r.CHR))
        if g is None:
            best.append(0.0)
            continue
        lo, hi = r.BP - ABC_WINDOW, r.BP + ABC_WINDOW
        s, e = g["start"], g["end"]
        near = ((s >= lo) & (s <= hi)) | ((e >= lo) & (e <= hi)) | ((s <= lo) & (e >= hi))
        best.append(g.loc[near, "abc_score"].max() if near.any() else 0.0)
    df[f"ABC_score_{cell}"] = best

df = df.drop(columns=["gene_TSS"])
df.to_parquet(f"{out_dir}/features.parquet", index=False)
print(f"wrote {out_dir}/features.parquet: {len(df)} variant-gene pairs, {df.shape[1]} columns")

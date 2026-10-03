#!/usr/bin/env python
"""
Score the featurized variant-gene pairs with the published scEEMS models of the data release
(model_training/models/{cohort}/). For each cell type, a pair is scored by the model that was trained
without the variant's chromosome, as the released predictions were. The feature matrix is built by
shared/featurize.py, as in training and inference: the features of 6_build_features.py, plus the variant
type, the gene constraint of the data release, gnomad_MAF (5_extract_gnomad_maf.py) and |gpn_star_llr|
(4_score_gpn_star.py).

usage:   python 7_score_variants.py OUT_DIR [--cell-types Mic Ast ...] [--models weighted_full ...]
         models: weighted_full (Weighted (Full), scEEMS; the default), unweighted_full, weighted_restricted
input:   OUT_DIR/features.parquet, gnomad_MAF.tsv, gpn_star.tsv
outputs: OUT_DIR/predictions.tsv          variant_id, SNP, gene_id, gene_name, pred_<cell type>_<model>
         OUT_DIR/feature_matrix.parquet   variant_id, gene_id and the 4,840 model features (the same
                                          features, in the same order, for every cell type's models)
"""
import argparse
import os
import pickle
import sys

import joblib
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
import featurize as F
from config import COHORTS, cohort_name, path

MODELS = ["weighted_full", "unweighted_full", "weighted_restricted"]

ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
ap.add_argument("out_dir")
ap.add_argument("--cell-types", nargs="+", default=[c.split("_")[0] for c in COHORTS])
ap.add_argument("--models", nargs="+", choices=MODELS, default=["weighted_full"])
a = ap.parse_args()

feat = pd.read_parquet(f"{a.out_dir}/features.parquet")
maf = pd.read_csv(f"{a.out_dir}/gnomad_MAF.tsv", sep="\t")
gpn = pd.read_csv(f"{a.out_dir}/gpn_star.tsv", sep="\t").dropna(subset=["gpn_star_llr"])
gpn = dict(zip(gpn["variant_id"], gpn["gpn_star_llr"]))
column_dict = F.load_columns_dict()
df = F.attach_aux(feat, dict(glof=F.load_gene_lof(), maf=maf))

out = df[["variant_id", "SNP", "gene_id", "gene_name"]].copy()
for cell in a.cell_types:
    cohort = cohort_name(cell)
    model_dir = path("published_models_dir", cohort=cohort)
    with open(f"{model_dir}/feature_cols.pkl", "rb") as fh:
        fc = pickle.load(fh)
    missing = [c for c in fc["cols"] if c not in df.columns]
    if missing:
        sys.exit(f"features.parquet lacks {len(missing)} model features, e.g. {missing[:5]}")
    X, FEATS, _, _ = F.build_X(df, column_dict, gpn, cols=fc["cols"], abscols=fc["abscols"])
    if cell == a.cell_types[0]:
        pd.concat([df[["variant_id", "gene_id"]], X], axis=1).to_parquet(f"{a.out_dir}/feature_matrix.parquet", index=False)
    for m in a.models:
        pred = pd.Series(index=df.index, dtype=float)
        for chrom, rows in df.groupby("_chrom").groups.items():
            pred[rows] = joblib.load(f"{model_dir}/{m}_{chrom}.joblib").predict_proba(X.loc[rows, FEATS])[:, 1]
        out[f"pred_{cell}_{m}"] = pred
out.to_csv(f"{a.out_dir}/predictions.tsv", sep="\t", index=False)
print(f"wrote {a.out_dir}/predictions.tsv ({len(out)} variant-gene pairs) and feature_matrix.parquet")

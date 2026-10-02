#!/usr/bin/env python
"""
TreeSHAP attribution of the scEEMS (weighted_full) model for one cell type and one chromosome, summed
over feature categories.

Variants explained: the model's predicted eQTLs on the chromosome, i.e. variant-gene pairs with pred_prob
above the cell type's tau* (the threshold selected by the S-LDSC analysis of step 9, read from
tau_star.json). Each chromosome is explained by the LOCO model that held it out, so the attributions
describe out-of-fold predictions, as in training and evaluation.

Features are rebuilt exactly as the model was trained (shared/featurize.py with the columns of
feature_cols.pkl). The |SHAP| of each variant is summed within 11 categories that partition the 4,840
features: ABC, CRE, chromBPNet, Enformer, TF, abs_gpn (GPN-STAR), gene_lof, variant_type, baseline,
distance, gnomad_MAF. Variants are split into promoter-like (|distance to TSS| <= 10 kb) and
enhancer-like (> 10 kb).

usage:   python shap_analysis.py COHORT CHR            (CHR: 1-22)
output:  {shap_dir}/shap_category_chr{CHR}.parquet: variant_id, gene_id, distance_TSS, region, tau_star,
         and one summed |SHAP| column per category
"""
import json
import os
import pickle
import sys

import joblib
import numpy as np
import pandas as pd
import pyarrow.dataset as ds
from catboost import Pool

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
import featurize as F
from config import cohort_name, path

cohort = cohort_name(sys.argv[1])
c = f"chr{int(sys.argv[2])}"
ENH_THRESH = 10000                      # |distance_TSS| > 10 kb = enhancer-like
CAT_ORDER = ["ABC", "CRE", "chromBPNet", "Enformer", "TF", "abs_gpn",
             "gene_lof", "variant_type", "baseline", "distance", "gnomad_MAF"]

with open(f"{path('aggregate_dir')}/tau_star.json") as fh:
    PRED_THRESH = json.load(fh)[cohort]
MDIR = path("model_dir", cohort=cohort)
PRED = f"{path('predictions_parquet_dir', cohort=cohort)}/weighted_full/predictions.parquet"
VAR = path("all_variants_dir", cohort=cohort)
MODEL = f"{MDIR}/weighted_full_{c}.joblib"
OUT_DIR = path("shap_dir", cohort=cohort)
for p in (PRED, MODEL, f"{MDIR}/feature_cols.pkl"):
    if not os.path.exists(p):
        sys.exit(f"missing input: {p}")
os.makedirs(OUT_DIR, exist_ok=True)

# predicted eQTLs on this chromosome (chr is the partition key of the prediction dataset)
preds = ds.dataset(PRED, partitioning="hive").to_table(
    columns=["variant_id", "gene_id", "pred_prob"],
    filter=(ds.field("chr") == c) & (ds.field("pred_prob") > PRED_THRESH)).to_pandas()
if preds.empty:
    sys.exit(f"[{cohort} {c}] no variants with pred_prob > {PRED_THRESH}")
print(f"[{cohort} {c}] predicted eQTLs (pred_prob > {PRED_THRESH}): {len(preds):,} over "
      f"{preds.gene_id.nunique():,} genes", flush=True)

# the trained feature matrix for exactly these variant-gene pairs, one gene at a time
with open(f"{MDIR}/feature_cols.pkl", "rb") as fh:
    fc = pickle.load(fh)
cols, abscols, FEATS = fc["cols"], fc["abscols"], fc["FEATS"]
aux = F.load_aux(maf_chrom=c)
gpn = F.load_gpn_map(chrom=c)
Xs, metas = [], []
for gene_id, grp in preds.groupby("gene_id", sort=False):
    gpath = f"{VAR}/{gene_id}"
    if not os.path.exists(gpath):
        print(f"  skip {gene_id}: no all_variants table", flush=True)
        continue
    try:
        df = F.load_gene_variants(gpath, gene_id, aux)
    except Exception as e:              # unreadable (e.g. truncated) table
        print(f"  skip {gene_id}: unreadable ({type(e).__name__})", flush=True)
        continue
    df = df.merge(grp[["variant_id"]], on="variant_id", how="inner")
    if df.empty:
        continue
    X, _, _, _ = F.build_X(df, aux["column_dict"], gpn, cols=cols, abscols=abscols)
    Xs.append(X)
    metas.append(df[["variant_id", "gene_id", "distance_TSS"]])
if not Xs:
    sys.exit(f"[{cohort} {c}] no variants could be featurized")
X = pd.concat(Xs, ignore_index=True)[FEATS]
meta = pd.concat(metas, ignore_index=True).reset_index(drop=True)
del Xs, metas
print(f"[{cohort} {c}] feature matrix {X.shape[0]:,} x {X.shape[1]:,}", flush=True)

# TreeSHAP with the LOCO model that held out this chromosome
clf = joblib.load(MODEL)
assert list(clf.feature_names_) == list(FEATS), "feature order differs from the model's"
sv = clf.get_feature_importance(Pool(X), type="ShapValues")[:, :-1]     # drop the expected-value column
np.abs(sv, out=sv)

# the 8 weighted categories plus baseline, distance and gnomad_MAF partition the features exactly
f2c = F.category_map(FEATS, aux["column_dict"])
base = set(aux["column_dict"].get("baseline", []))
dist = set(aux["column_dict"].get("distance", []))
for f in FEATS:
    if f not in f2c:
        f2c[f] = "gnomad_MAF" if f == "gnomad_MAF" else "baseline" if f in base else "distance" if f in dist else "other"
unassigned = [f for f in FEATS if f2c[f] == "other"]
assert not unassigned, f"{len(unassigned)} features outside every category: {unassigned[:10]}"
idx = {cat: [i for i, f in enumerate(FEATS) if f2c[f] == cat] for cat in CAT_ORDER}
assert sum(len(v) for v in idx.values()) == len(FEATS), "categories do not partition the features"

out = meta.copy()
out["region"] = np.where(np.abs(out["distance_TSS"]) > ENH_THRESH, "enhancer", "promoter")
out["tau_star"] = PRED_THRESH
for cat in CAT_ORDER:
    out[cat] = sv[:, idx[cat]].sum(axis=1)
out.to_parquet(f"{OUT_DIR}/shap_category_{c}.parquet", index=False)
n_enh = int((out.region == "enhancer").sum())
print(f"[{cohort} {c}] enhancer-like {n_enh:,} | promoter-like {len(out) - n_enh:,} -> "
      f"{OUT_DIR}/shap_category_{c}.parquet", flush=True)

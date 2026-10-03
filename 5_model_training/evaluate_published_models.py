#!/usr/bin/env python
"""
Evaluate the published scEEMS models of one cell type on the held-out test sets of the data release.

The data release contains the models of the manuscript (model_training/models/{cohort}/): for every
held-out chromosome the three models weighted_full (Weighted (Full), scEEMS), unweighted_full and
weighted_restricted, and feature_cols.pkl, the feature columns they were trained on. Each chromosome's
test set is scored by the models that held that chromosome out, with features built exactly as in training
(shared/featurize.py), and the pooled AUPRC over all chromosomes is reported for each model: the held-out
AUPRCs of the manuscript. No training is needed; a cell type takes a few minutes.

usage:   python evaluate_published_models.py COHORT            COHORT: Mic or Mic_mega_eQTL
outputs: {published_models_eval_dir}/test_pred_chr{N}.parquet    same columns as train_loco.py
         {published_models_eval_dir}/auprc.json                  pooled AUPRC, all variants and SNVs only
"""
import json
import os
import pickle
import sys

import joblib
import pandas as pd
from sklearn.metrics import average_precision_score

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
import featurize as F
from config import cohort_name, path

MODELS = ["weighted_full", "unweighted_full", "weighted_restricted"]

cohort = cohort_name(sys.argv[1])
model_dir = path("published_models_dir", cohort=cohort)
out_dir = path("published_models_eval_dir", cohort=cohort)
os.makedirs(out_dir, exist_ok=True)
with open(f"{model_dir}/feature_cols.pkl", "rb") as fh:
    fc = pickle.load(fh)

parts = []
for c in F.ALL_CHR:
    aux = F.load_aux(maf_chrom=c)
    dft = F.load_training_tables(cohort, "test", [c], aux)
    if not len(dft):                        # microglia chr21 has no test set
        print(f"[{cohort} {c}] no test set", flush=True)
        continue
    X, _, _, _ = F.build_X(dft, aux["column_dict"], F.load_gpn_map(chrom=c), cols=fc["cols"], abscols=fc["abscols"])
    part = pd.DataFrame({"chrom": c,
                         "variant_id": dft["variant_id"].to_numpy(),
                         "gene_id": dft["gene_id"].to_numpy(),
                         "y_true": dft["label"].to_numpy(),
                         "is_snv": dft["is_SNP"].to_numpy()})
    for m in MODELS:
        part[f"pred_{m}"] = joblib.load(f"{model_dir}/{m}_{c}.joblib").predict_proba(X)[:, 1]
    part.to_parquet(f"{out_dir}/test_pred_{c}.parquet", index=False)
    parts.append(part)
    print(f"[{cohort} {c}] scored {len(part)} test variant-gene pairs", flush=True)

df = pd.concat(parts, ignore_index=True)
y = df["y_true"].to_numpy().astype(int)
snv = df["is_snv"].to_numpy().astype(bool)
res = {"cohort": cohort, "n_chromosomes": len(parts), "n_test": int(len(df)), "n_pos": int(y.sum()),
       "pooled_auprc": {m: float(average_precision_score(y, df[f"pred_{m}"])) for m in MODELS},
       "pooled_auprc_snv": {m: float(average_precision_score(y[snv], df.loc[snv, f"pred_{m}"])) for m in MODELS}}
with open(f"{out_dir}/auprc.json", "w") as fh:
    json.dump(res, fh, indent=2)
print(json.dumps(res, indent=2))

#!/usr/bin/env python
"""
Leave-one-chromosome-out (LOCO) training for one cell type and one held-out chromosome.

Trains the three models of the manuscript on every chromosome except the held-out one and saves them for
inference (step 6):
    weighted_full        Weighted (Full)        train/,            selected category weights (scEEMS)
    unweighted_full      Unweighted (Full)      train/,            every category weight 1
    weighted_restricted  Weighted (Restricted)  train_restricted/, selected category weights
The category weights are read from best_configs_{cohort}.json: by default the file of the data release
(config setting feature_weights_file), or the file given as the third argument, e.g. the output of
select_feature_weights.py. A held-out odd chromosome uses the weights selected on the even chromosomes and
vice versa, so the weights never saw the chromosome they are evaluated on.

If the held-out chromosome has a test set (PIP > 0.9 positives), the three models score it for the
held-out AUPRC (evaluate_auprc.py). Microglia chr21 has no test set.

usage:   python train_loco.py COHORT CHR [WEIGHTS_JSON]       COHORT: Mic or Mic_mega_eQTL; CHR: 1-22
outputs: {model_dir}/{model}_chr{CHR}.joblib, {model_dir}/feature_cols.pkl
         {test_predictions_dir}/test_pred_chr{CHR}.parquet
"""
import json
import os
import pickle
import sys
import tempfile
import time

import joblib
import numpy as np
import pandas as pd
from catboost import CatBoostClassifier

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
import featurize as F
from config import cohort_name, path

cohort = cohort_name(sys.argv[1])
c = f"chr{int(sys.argv[2])}"                                   # held-out chromosome
train_chrs = [x for x in F.ALL_CHR if x != c]

MODELS = path("model_dir", cohort=cohort)
TESTPRED = path("test_predictions_dir", cohort=cohort)
for d in (MODELS, TESTPRED):
    os.makedirs(d, exist_ok=True)

weights_file = sys.argv[3] if len(sys.argv) > 3 else path("feature_weights_file", cohort=cohort)
with open(weights_file) as fh:
    cfg = json.load(fh)
weights = cfg["even"]["best"]["weights"] if int(c[3:]) % 2 == 1 else cfg["odd"]["best"]["weights"]

# ------------------------------------------------------------------ training data (all but c)
t0 = time.time()
aux = F.load_aux()
gpn = F.load_gpn_map()
df = F.load_training_tables(cohort, "train", train_chrs, aux)
cols, abscols = F.compute_cols(df, aux["column_dict"])
FEATS = list(cols) + [F.GPN_COL]
fd, tmp = tempfile.mkstemp(dir=MODELS, suffix=".pkl")         # atomic: 22 jobs write the same file
os.close(fd)
with open(tmp, "wb") as fh:
    pickle.dump({"cols": cols, "abscols": abscols, "FEATS": FEATS}, fh)
os.replace(tmp, f"{MODELS}/feature_cols.pkl")                  # inference reindexes to this exact order

X, _, _, _ = F.build_X(df, aux["column_dict"], gpn, cols=cols, abscols=abscols)
y, w = df["label"].to_numpy(), F.sample_weights(df)
print(f"[{cohort} {c}] train n={len(df):,} features={len(FEATS)} ({time.time() - t0:.0f}s)", flush=True)


def fit(category_weights, Xtr, ytr, wtr):
    fw = F.feature_weights(category_weights, FEATS, aux["column_dict"])
    return CatBoostClassifier(**F.CONS, feature_weights=fw).fit(Xtr, ytr, sample_weight=wtr)


# ------------------------------------------------------------------ the three models
models = {"unweighted_full": fit(F.NULL_WEIGHTS, X, y, w),
          "weighted_full": fit(weights, X, y, w)}
dfr = F.load_training_tables(cohort, "train_restricted", train_chrs, aux)
if len(dfr) and dfr["label"].sum():
    Xr, _, _, _ = F.build_X(dfr, aux["column_dict"], gpn, cols=cols, abscols=abscols)
    models["weighted_restricted"] = fit(weights, Xr, dfr["label"].to_numpy(), F.sample_weights(dfr))
for name, clf in models.items():
    joblib.dump(clf, f"{MODELS}/{name}_{c}.joblib")
print(f"[{cohort} {c}] saved {sorted(models)} ({time.time() - t0:.0f}s)", flush=True)

# ------------------------------------------------------------------ held-out test set
dft = F.load_training_tables(cohort, "test", [c], aux)
if not len(dft):
    print(f"[{cohort} {c}] no test set ({time.time() - t0:.0f}s)", flush=True)
    sys.exit(0)

Xte, _, _, _ = F.build_X(dft, aux["column_dict"], gpn, cols=cols, abscols=abscols)
out = pd.DataFrame({"chrom": c,
                    "variant_id": dft["variant_id"].to_numpy(),
                    "gene_id": dft["gene_id"].to_numpy(),
                    "y_true": dft["label"].to_numpy(),
                    "is_snv": dft["is_SNP"].to_numpy()})
for name in ("weighted_full", "unweighted_full", "weighted_restricted"):
    out[f"pred_{name}"] = models[name].predict_proba(Xte)[:, 1] if name in models else np.nan
out.to_parquet(f"{TESTPRED}/test_pred_{c}.parquet", index=False)
print(f"[{cohort} {c}] scored {len(out)} test variant-gene pairs ({time.time() - t0:.0f}s)", flush=True)

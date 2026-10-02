#!/usr/bin/env python
"""
Score the variant-gene pairs of the non-European MAGMA panel with the scEEMS (weighted_full) model, for
one cell type and one chromosome.

The non-European panel is the set of variants in the ADSP African, admixed American and East Asian
reference panels, paired with genes and featurized like the main variant set (multi_ancestry_predictions_dir:
MAGMA_features_{cohort}_chr{N}.parquet with the features, MAGMA_predictions_{cohort}_chr{N}.tsv.gz with the
row-aligned variant and gene columns). The features do not depend on the model, so only the scoring is done
here: the GPN-STAR feature is joined from the GPN-STAR scores of the data release and the scores of the
non-European SNVs (noneur_gpn_star/), the columns are put in the model's order, and the pairs are scored by
the model that held out this chromosome. make_magma_files.py adds the pairs that pass tau* to the
prediction annotation.

usage:   python score_noneur.py COHORT CHR
output:  {magma_dir}/noneur_predictions/{cohort}/MAGMA_predictions_{cohort}_chr{CHR}.tsv.gz
"""
import os
import sys

import joblib
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
import featurize as F
from config import cohort_name, path

cohort = cohort_name(sys.argv[1])
c = f"chr{int(sys.argv[2])}"

SRC = f"{path('multi_ancestry_predictions_dir')}/{cohort}"
FEAT = f"{SRC}/MAGMA_features_{cohort}_{c}.parquet"           # features, row-aligned with META
META = f"{SRC}/MAGMA_predictions_{cohort}_{c}.tsv.gz"         # variant and gene columns
NONEUR_GPN = path("noneur_gpn_star_file")
MODEL = f"{path('model_dir', cohort=cohort)}/weighted_full_{c}.joblib"
OUT_DIR = f"{path('magma_dir')}/noneur_predictions/{cohort}"
OUT = f"{OUT_DIR}/MAGMA_predictions_{cohort}_{c}.tsv.gz"
META_OUT = ["CHR", "BP", "SNP", "REF", "ALT", "gene_id", "TargetGene", "distance_TSS", "TargetGeneTSS", "ABC.Score"]

for p in (FEAT, META, MODEL, NONEUR_GPN):
    if not os.path.exists(p):
        sys.exit(f"missing input: {p}")
os.makedirs(OUT_DIR, exist_ok=True)

Xbase = pd.read_parquet(FEAT)
meta = pd.read_csv(META, sep="\t")
assert len(Xbase) == len(meta), f"row mismatch: features {len(Xbase)} vs variants {len(meta)}"

# GPN-STAR: scores of the data release for this chromosome, plus the non-European SNV scores
gpn_map = F.load_gpn_map(chrom=c)
nz = pd.read_parquet(NONEUR_GPN, columns=["variant_id", F.GPN_COL])
gpn_map.update(dict(zip(nz["variant_id"].to_numpy(), nz[F.GPN_COL].to_numpy())))
variant_id = meta["CHR"].astype(str) + ":" + meta["BP"].astype(str) + ":" + meta["REF"] + ":" + meta["ALT"]
gpn = np.abs(variant_id.map(gpn_map).to_numpy())                 # NaN for variants without a score

# the model's feature order, gpn_star_llr last
clf = joblib.load(MODEL)
FEATS = list(clf.feature_names_)
X = Xbase.reindex(columns=[f for f in FEATS if f != F.GPN_COL], fill_value=0).copy()
X[F.GPN_COL] = gpn
meta["pred_prob"] = clf.predict_proba(X[FEATS])[:, 1]

meta[META_OUT + ["pred_prob"]].to_csv(OUT, sep="\t", index=False, compression="gzip")
print(f"[{cohort} {c}] {len(meta):,} variant-gene pairs | GPN-STAR score for {np.isfinite(gpn).mean():.1%} "
      f"-> {OUT}", flush=True)

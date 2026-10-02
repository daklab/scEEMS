#!/usr/bin/env python
"""
Score every cis-variant of one gene with the three models of step 5.

For the gene in row IDX of the cell type's gene list (make_gene_lists.py), the gene's all_variants table
(step 3) is featurized exactly as in training (shared/featurize.py, with the columns and order saved in
feature_cols.pkl) and scored by the models that held out the gene's chromosome, so no gene is scored by a
model that saw its chromosome in training.

usage:   python model_inference.py COHORT IDX          (IDX = 1..n_genes, a row of the gene list)
outputs: {predictions_dir}/{model}/{gene_id}_predictions.tsv for each model, with columns
         variant_id, chr, pos, ref, alt, pip, gene_id, pred_prob
"""
import os
import pickle
import sys

import joblib
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
import featurize as F
from config import cohort_name, path

MODELS = ["weighted_full", "unweighted_full", "weighted_restricted"]

cohort = cohort_name(sys.argv[1])
idx = int(sys.argv[2])
genes = pd.read_csv(path("gene_list_file", cohort=cohort), sep="\t")
row = genes.iloc[idx - 1]                                       # 1-based array index -> row
gene_id, gene_chr = row["gene_id"], row["chr"]
model_dir = path("model_dir", cohort=cohort)
pred_dir = path("predictions_dir", cohort=cohort)

# featurize the gene's variants with the training-time feature columns (chromosome-specific MAF + GPN-STAR)
with open(f"{model_dir}/feature_cols.pkl", "rb") as fh:
    fc = pickle.load(fh)
aux = F.load_aux(maf_chrom=gene_chr)
gpn = F.load_gpn_map(chrom=gene_chr)
df = F.load_gene_variants(f"{path('all_variants_dir', cohort=cohort)}/{gene_id}", gene_id, aux)
X, _, _, _ = F.build_X(df, aux["column_dict"], gpn, cols=fc["cols"], abscols=fc["abscols"])

info = df[["variant_id", "chr", "pos", "ref", "alt", "pip", "gene_id"]].reset_index(drop=True)
scored = []
for model in MODELS:
    mf = f"{model_dir}/{model}_{gene_chr}.joblib"
    if not os.path.exists(mf):
        print(f"[{cohort} {gene_id} {gene_chr}] no {mf} -- {model} skipped", flush=True)
        continue
    out = info.copy()
    out["pred_prob"] = joblib.load(mf).predict_proba(X)[:, 1]
    os.makedirs(f"{pred_dir}/{model}", exist_ok=True)
    out.to_csv(f"{pred_dir}/{model}/{gene_id}_predictions.tsv", sep="\t", index=False)
    scored.append(model)
print(f"[{cohort} {gene_id} {gene_chr}] scored {len(X)} variants with {scored}", flush=True)

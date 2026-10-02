#!/usr/bin/env python
"""
Size-matched S-LDSC, step 1: for one cell type, pick the top N baseline variants under each ranking, so
that every annotation of the comparison contains exactly N variants and has the same Prop._SNPs.

Why: enrichment grows as an annotation gets smaller, so comparing annotations of different sizes (e.g. a
credible-set annotation with tens of thousands of variants against a predicted-eQTL annotation with a few
thousand) partly measures size rather than how well each one localizes heritability. N = 5,000 is smaller
than every credible-set annotation of every cell type, so no ranking runs out of variants.

Ranking happens among the baseline variants (after the join), so every annotation ends up with exactly N
variants. Per chromosome the top KEEP_PER_CHR variants by each score are kept, then the genome-wide top N
is taken from the pooled set, which is exact because KEEP_PER_CHR >= N.

pred_prob is not on the same scale across models (the weighted_restricted model predicts a much rarer
event), so each model is ranked on its own scores; rankings are never pooled across models.

Rankings (maximum over genes per variant in every case):
  pred_weighted_full, pred_unweighted_full, pred_weighted_restricted   the three models, eQTL + other genes
  pip                   the FunGen-xQTL fine-mapping PIP, eQTL + other genes
  pip_{prior}           the credible-set members of each fine-mapping prior (step 11), eQTL genes that
                        were fine-mapped under all five priors

usage:   python select_topn.py COHORT [N=5000]
output:  {sldsc_dir}/top{N}_selection/{cohort}_top{N}.tsv.gz (CHR, BP, A1, A2, score, ranking)
"""
import os
import sys

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import cohort_name, path

cohort = cohort_name(sys.argv[1])
cell = cohort.replace("_mega_eQTL", "")
N = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
KEEP_PER_CHR = max(20000, 4 * N)

BASELINE = path("baseline_annot_dir")
BANK = path("predictions_parquet_dir", cohort=cohort)
FM = f"{path('aggregate_dir')}/finemapping_comparison"
OUT_DIR = f"{path('sldsc_dir')}/top{N}_selection"
os.makedirs(OUT_DIR, exist_ok=True)

MODELS = {"pred_weighted_full": "weighted_full",
          "pred_unweighted_full": "unweighted_full",
          "pred_weighted_restricted": "weighted_restricted"}
PRIORS = ["uniform", "EMS", "scEEMS_Weighted_Full", "scEEMS_Weighted_Restricted", "scEEMS_Unweighted_Full"]

# baseline keys per chromosome, so ranking happens only among variants the annotation can contain
base = {}
for ch in range(1, 23):
    b = pd.read_csv(f"{BASELINE}/baseline_chr{ch}.annot.gz", sep="\t", usecols=["BP", "A1", "A2"])
    base[ch] = b.drop_duplicates()


def top_from_bank(model, score_col):
    """Max score per variant, restricted to the baseline, top KEEP_PER_CHR per chromosome."""
    parts = []
    for ch in range(1, 23):
        t = ds.dataset(f"{BANK}/{model}/predictions.parquet", format="parquet", partitioning="hive").to_table(
            columns=["pos", "ref", "alt", score_col], filter=ds.field("chr") == f"chr{ch}").to_pandas()
        t = t.rename(columns={"pos": "BP", "alt": "A1", "ref": "A2", score_col: "score"})
        t = t.groupby(["BP", "A1", "A2"], as_index=False)["score"].max()
        t = base[ch].merge(t, on=["BP", "A1", "A2"])
        parts.append(t.nlargest(min(KEEP_PER_CHR, len(t)), "score").assign(CHR=ch))
    return pd.concat(parts, ignore_index=True)


def top_from_cs(prior):
    """Max PIP per variant among the prior's credible-set members, restricted to the baseline."""
    cons = pd.read_csv(f"{FM}/gene_consensus.tsv", sep="\t")
    cons = set(cons.loc[cons["in_all_priors"] & (cons["cell_type"] == cell), "gene_id"])
    d = pd.read_csv(f"{FM}/{cell}_{prior}_cs_variants.tsv", sep="\t")
    d = d[d["gene_id"].isin(cons)]
    p = d["variant_id"].str.split(":", expand=True)
    d = d.assign(CHR=p[0].str.replace("chr", "", regex=False), BP=p[1].astype(np.int64),
                 A2=p[2], A1=p[3], score=d["pip"])
    d = d[d["CHR"].str.isdigit()].astype({"CHR": int})
    d = d.groupby(["CHR", "BP", "A1", "A2"], as_index=False)["score"].max()
    return pd.concat([base[ch].merge(d[d.CHR == ch], on=["BP", "A1", "A2"]) for ch in range(1, 23)],
                     ignore_index=True)


rankings = {k: top_from_bank(v, "pred_prob") for k, v in MODELS.items()}
rankings["pip"] = top_from_bank("weighted_full", "pip")     # the PIP is the same in every model's dataset
for pr in PRIORS:
    rankings[f"pip_{pr}"] = top_from_cs(pr)

out = []
for name, df in rankings.items():
    sel = df.nlargest(N, "score")
    if len(sel) < N:
        print(f"  WARNING {name}: only {len(sel):,} candidates, fewer than N={N}", flush=True)
    out.append(sel.assign(ranking=name)[["CHR", "BP", "A1", "A2", "score", "ranking"]])
    print(f"  {name:<34} selected {len(sel):>6,}  score cutoff {sel['score'].min():.6g}", flush=True)

res = pd.concat(out, ignore_index=True)
dest = f"{OUT_DIR}/{cohort}_top{N}.tsv.gz"
res.to_csv(dest, sep="\t", index=False, compression="gzip")
print(f"wrote {dest}  ({len(res):,} rows, {res.ranking.nunique()} rankings)")

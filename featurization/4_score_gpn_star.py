#!/usr/bin/env python
"""
GPN-STAR scores of the prepared SNVs. The models use |gpn_star_llr|, the log-likelihood ratio of ALT versus
REF under the primate GPN-STAR model (gpn-star-hg38-p243-200m), which scored the training variants with
GPN-STAR's own variant scoring (vep, window 256). Each SNV's score is
  - its training score, if it is a training variant (data release: gpn_star_scores_all.parquet);
  - otherwise the score from GPN-STAR's published genome-wide scores of the same model. These are
    calibrated for mutation rate, and the calibration is undone with the model's calibration table:
        gpn_star_llr = llr_calibrated + llr_neutral_mean[5-bp reference context centred on the SNV, ALT]
    On chromosome 22 this reproduces the training scores to within 0.01 (the published scores are rounded).
Insertions, deletions, multi-base substitutions and SNVs without a published score get no score (NaN),
which the models treat as missing, as in training.

usage:   python 4_score_gpn_star.py OUT_DIR --scores 'DIR/llr_{chrom}.parquet' --calibration llr.parquet --fasta GRCh38.fa
         --scores: the per-chromosome LLR files of gpn-star-hg38-p243-200m (songlab/gpn-star-scores);
         --calibration: calibration_table/llr.parquet of songlab/gpn-star-hg38-p243-200m (README)
input:   OUT_DIR/variants.tsv (1_prepare_variants.py)
output:  OUT_DIR/gpn_star.tsv   variant_id, gpn_star_llr, source (training, gpn-star-scores or none)
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
import featurize as F

ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
ap.add_argument("out_dir")
ap.add_argument("--scores", required=True, help="per-chromosome LLR files, with {chrom} for the chromosome (e.g. chr7)")
ap.add_argument("--calibration", required=True, help="calibration_table/llr.parquet of the model")
ap.add_argument("--fasta", required=True, help="GRCh38 FASTA (featurization/README.md)")
a = ap.parse_args()
if "{chrom}" not in a.scores:
    sys.exit("--scores needs {chrom}, e.g. 'gpn-star-scores/data/gpn-star-hg38-p243-200m/llr/llr_{chrom}.parquet'")

variants = pd.read_csv(f"{a.out_dir}/variants.tsv", sep="\t")
neutral = pd.read_parquet(a.calibration).set_index("pentanuc_mut")["llr_neutral_mean"]
fa = pysam.FastaFile(a.fasta)
snv = (variants["REF"].str.len() == 1) & (variants["ALT"].str.len() == 1)
out = []
for chrom, v in variants[snv].groupby("CHR"):
    v = v.copy()
    training = F.load_gpn_map(chrom=chrom)
    swapped = v["CHR"] + ":" + v["BP"].astype(str) + ":" + v["ALT"] + ":" + v["REF"]
    v["gpn_star_llr"] = v["variant_id"].map(training)
    v.loc[v["gpn_star_llr"].isna(), "gpn_star_llr"] = -swapped.map(training)     # the training ID had REF/ALT swapped
    v["source"] = np.where(v["gpn_star_llr"].notna(), "training", "none")
    todo = v["gpn_star_llr"].isna()
    if todo.any():
        pos = sorted(v.loc[todo, "BP"].unique().tolist())
        s = pq.read_table(a.scores.format(chrom=chrom), columns=["pos", "ref", "alt", "llr_calibrated"],
                          filters=[("pos", "in", pos)]).to_pandas().drop_duplicates(["pos", "ref", "alt"])
        w = v[todo].merge(s, left_on=["BP", "REF", "ALT"], right_on=["pos", "ref", "alt"], how="left")
        context = [fa.fetch(chrom, p - 3, p + 2).upper() + "_" + alt for p, alt in zip(w["BP"], w["ALT"])]
        llr = (w["llr_calibrated"].astype(float) + pd.Series(context).map(neutral).astype(float)).to_numpy()
        v.loc[todo, "gpn_star_llr"] = llr
        v.loc[todo, "source"] = np.where(np.isnan(llr), "none", "gpn-star-scores")
    out.append(v[["variant_id", "gpn_star_llr", "source"]])
if (~snv).any():
    out.append(variants.loc[~snv, ["variant_id"]].assign(gpn_star_llr=np.nan, source="none"))
out = pd.concat(out) if out else pd.DataFrame(columns=["variant_id", "gpn_star_llr", "source"])
out.to_csv(f"{a.out_dir}/gpn_star.tsv", sep="\t", index=False)
print(f"wrote {a.out_dir}/gpn_star.tsv: " + ", ".join(f"{n} {k}" for k, n in out["source"].value_counts().items()))

#!/usr/bin/env python
"""
ChromBPNet variant effects of the prepared variants, computed as for the training data, in four brain cell
types (microglia, astrocyte, oligodendrocyte, neuron) for chromatin accessibility (ATAC: ChromBPNet) and
two histone marks (H3K27ac, H3K4me3: BPNet):
  - a variant is paired with every peak of the cell type and assay whose centre is within 1,024 bp;
  - the 2,114 bp sequence centred on the peak centre is predicted with the reference and with the
    alternate allele by each of the 5 fold models;
  - per model, log_counts_diff = alternate - reference predicted log counts,
    log_probs_diff_abs_sum = sum over positions of |log alternate - log reference profile probability|
    and probs_jsd_diff = Jensen-Shannon distance of the two profiles, both signed like log_counts_diff;
    each score is averaged over the 5 models.
6_build_features.py summarizes a variant's peaks per cell type, assay and score (largest absolute value).

usage:   python 3_score_chrombpnet.py OUT_DIR --fasta GRCh38.fa [--device cuda|cpu]
input:   OUT_DIR/variants.tsv; chrombpnet_models/ and chrombpnet_peaks.tsv.gz of the data release
output:  OUT_DIR/chrombpnet.tsv   one row per variant, cell type, assay and peak
Run in the scEEMS_chrombpnet environment (environment_chrombpnet.yml). A GPU is used if available.
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import torch
from bpnetlite import BPNet
from bpnetlite.bpnet import ControlWrapper
from scipy.spatial.distance import jensenshannon
from tangermeme.utils import one_hot_encode

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
import variant_sequence as VS
from config import path

CELLS = ["microglia", "astrocyte", "oligodendrocyte", "neuron"]
ASSAYS = ["ATAC", "H3K27ac", "H3K4me3"]
SEQUENCE_LENGTH = 2114
MAX_DISTANCE = 1024                                  # variant to peak centre
IGNORE = list("QWERYUIOPSDFHJKLZXVBNM")              # one-hot encoded as 0 (e.g. N)
SCORES = ["log_counts_diff_chrombpnet", "log_probs_diff_abs_sum_chrombpnet", "probs_jsd_diff_chrombpnet"]

ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
ap.add_argument("out_dir")
ap.add_argument("--fasta", required=True, help="GRCh38 FASTA (featurization/README.md)")
ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
a = ap.parse_args()
torch.set_float32_matmul_precision("high")
torch.backends.cuda.matmul.allow_tf32 = True

variants = pd.read_csv(f"{a.out_dir}/variants.tsv", sep="\t")
genome = VS.Genome(a.fasta)
peaks = pd.read_csv(path("chrombpnet_peaks_file"), sep="\t")
peaks["peak_center"] = peaks["start"] + (peaks["end"] - peaks["start"]) // 2
model_dir = path("chrombpnet_models_dir")


def load_models(cell, assay):
    """The 5 fold models: ChromBPNet without its bias model (ATAC) or BPNet (histone marks)."""
    if assay == "ATAC":
        return [BPNet.from_chrombpnet(f"{model_dir}/ATAC_{cell}_fold_{f}.chrombpnet_nobias.h5").to(a.device)
                for f in range(5)]
    return [ControlWrapper(torch.load(f"{model_dir}/{assay}_{cell}_fold_{f}.final.torch", map_location=a.device,
                                      weights_only=False)).to(a.device) for f in range(5)]


def encode(seq):
    return one_hot_encode(seq, ignore=IGNORE).to(torch.float32).unsqueeze(0).to(a.device)


def predict(model, x):
    """Profile probabilities (1, length) and predicted log counts (1,)."""
    profile, counts = model(x)
    return (torch.softmax(profile.squeeze(1), dim=1).detach().cpu().numpy(),
            counts.detach().cpu().numpy().reshape(1, ))


def variant_effects(models, x_ref, x_alt):
    per_model = []
    with torch.no_grad():
        for m in models:
            ref_prob, ref_count = predict(m, x_ref)
            alt_prob, alt_count = predict(m, x_alt)
            counts = alt_count - ref_count
            probs = np.sum(np.abs(np.log(alt_prob) - np.log(ref_prob)), axis=1) * np.sign(counts)
            jsd = np.array([jensenshannon(x, y) for x, y in zip(alt_prob, ref_prob)]) * np.sign(counts)
            per_model.append((counts, probs, jsd))
    return [np.mean([p[k] for p in per_model], axis=0).item() for k in range(3)]


# ---- variant-peak pairs: peak centre within 1,024 bp of the variant
pairs = []
for (cell, assay, chrom), p in peaks.groupby(["cell_type", "assay", "chrom"]):
    v = variants[variants["CHR"] == chrom]
    if v.empty:
        continue
    p = p.sort_values("peak_center")
    centers = p["peak_center"].to_numpy()
    for r in v.itertuples():
        lo = np.searchsorted(centers, r.BP - MAX_DISTANCE, "left")
        hi = np.searchsorted(centers, r.BP + MAX_DISTANCE, "right")
        for pk in p.iloc[lo:hi].itertuples():
            pairs.append((cell, assay, r.variant_id, r.CHR, r.BP, r.REF, r.ALT, pk.start, pk.end, pk.peak_center))
pairs = pd.DataFrame(pairs, columns=["cell_type", "assay", "variant_id", "CHR", "BP", "REF", "ALT",
                                     "peak_start", "peak_end", "peak_center"])
pairs["distance"] = pairs["BP"] - pairs["peak_center"]
print(f"{len(variants)} variants, {len(pairs)} variant-peak pairs | device {a.device}", flush=True)

# ---- scores, one cell type and assay (5 models) at a time
rows, t0 = [], time.time()
for (cell, assay), g in pairs.groupby(["cell_type", "assay"], sort=False):
    models = load_models(cell, assay)
    for r in g.itertuples(index=False):
        start, end = VS.window(r.peak_center - 1, r.peak_center, SEQUENCE_LENGTH)
        reference, alternate = VS.sequences(genome, r.CHR, start, end, r.BP, r.REF, r.ALT)
        rows.append(list(r) + variant_effects(models, encode(reference), encode(alternate)))
    print(f"  {cell} {assay}: {len(g)} pairs ({time.time() - t0:.0f} s)", flush=True)
    del models
pd.DataFrame(rows, columns=list(pairs.columns) + SCORES).to_csv(f"{a.out_dir}/chrombpnet.tsv", sep="\t", index=False)
print(f"wrote {a.out_dir}/chrombpnet.tsv")

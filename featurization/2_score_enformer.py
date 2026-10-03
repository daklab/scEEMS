#!/usr/bin/env python
"""
Enformer variant effects of the prepared variants, computed as for the training data: the published
Enformer model (TensorFlow Hub deepmind/enformer/1) predicts the 393,216 bp reference and alternate
sequences centred on the variant, and for each of the 5,313 human output tracks
    diff_32 = sum over the central 32 output bins (4,096 bp) of the reference prediction
              minus the same sum for the alternate prediction.
The sequences are built as in the official Enformer usage notebook (variant_sequence.py).

usage:   python 2_score_enformer.py OUT_DIR --fasta GRCh38.fa [--model HANDLE_OR_DIR]
input:   OUT_DIR/variants.tsv (1_prepare_variants.py); the track list targets_human.txt (data release)
output:  OUT_DIR/enformer.parquet   CHR, SNP, BP, REF, ALT and diff_32_<track identifier> for all tracks
Run in the scEEMS_enformer environment (environment_enformer.yml). A GPU is used when TensorFlow finds
one (a few seconds per variant); on CPU allow about a minute per variant.
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import tensorflow as tf
import tensorflow_hub as hub

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
import variant_sequence as VS
from config import path

SEQUENCE_LENGTH = 393_216
N_BINS = 896
CENTRAL_BINS = 32

ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
ap.add_argument("out_dir")
ap.add_argument("--fasta", required=True, help="GRCh38 FASTA (featurization/README.md)")
ap.add_argument("--model", default="https://tfhub.dev/deepmind/enformer/1",
                help="TensorFlow Hub handle of Enformer, or the directory of a downloaded copy")
a = ap.parse_args()

targets = pd.read_csv(path("targets_file"), sep="\t")
columns = [f"diff_32_{i}" for i in targets["identifier"]]
variants = pd.read_csv(f"{a.out_dir}/variants.tsv", sep="\t")
genome = VS.Genome(a.fasta)
model = hub.load(a.model).model
print(f"{len(variants)} variants | GPUs: {len(tf.config.list_physical_devices('GPU'))}", flush=True)

lo, hi = N_BINS // 2 - CENTRAL_BINS // 2, N_BINS // 2 + CENTRAL_BINS // 2
diffs, t0 = [], time.time()
for i, v in enumerate(variants.itertuples(), 1):
    start, end = VS.window(v.BP - 1, v.BP - 1, SEQUENCE_LENGTH)
    reference, alternate = VS.sequences(genome, v.CHR, start, end, v.BP, v.REF, v.ALT)
    ref = model.predict_on_batch(VS.one_hot(reference)[np.newaxis])["human"].numpy()     # (1, 896, 5313)
    alt = model.predict_on_batch(VS.one_hot(alternate)[np.newaxis])["human"].numpy()
    diffs.append((ref[:, lo:hi, :].sum(axis=1) - alt[:, lo:hi, :].sum(axis=1))[0])
    if i % 10 == 0 or i == len(variants):
        print(f"  {i}/{len(variants)} variants ({(time.time() - t0) / i:.1f} s each)", flush=True)

out = pd.concat([variants[["CHR", "SNP", "BP", "REF", "ALT"]].reset_index(drop=True),
                 pd.DataFrame(np.array(diffs, dtype=np.float32).reshape(len(variants), len(columns)), columns=columns)],
                axis=1)
out.to_parquet(f"{a.out_dir}/enformer.parquet", index=False)
print(f"wrote {a.out_dir}/enformer.parquet ({len(out)} variants x {len(columns)} tracks)")

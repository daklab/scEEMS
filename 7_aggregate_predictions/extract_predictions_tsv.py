#!/usr/bin/env python
"""
Export the scEEMS (weighted_full) predictions of one cell type and one chromosome as a bgzipped,
tabix-indexed TSV, the format of the predictions/ folder of the data release:
    #CHROM  POS  ID  REF  ALT  GENE_ID  PIP  PRED_PROBABILITY
CHROM has no "chr" prefix; rows are sorted by position, gene and variant; PIP is the FunGen-xQTL
fine-mapping PIP carried through from step 1. Query a region with, e.g.,
    tabix predictions_Mic_mega_eQTL_11.tsv.gz 11:86000000-86200000

usage:   python extract_predictions_tsv.py CHR COHORT          (CHR: 1-22)
output:  {release_export_dir}/predictions/{cohort}/predictions_{cohort}_{CHR}.tsv.gz (+ .tbi)
"""
import os
import sys

import pyarrow.dataset as ds
import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import cohort_name, path

chrom = sys.argv[1]
cohort = cohort_name(sys.argv[2])
bank_dir = f"{path('predictions_parquet_dir', cohort=cohort)}/weighted_full/predictions.parquet/chr=chr{chrom}"
write_dir = f"{path('release_export_dir')}/predictions/{cohort}"
os.makedirs(write_dir, exist_ok=True)

output_file = f"{write_dir}/predictions_{cohort}_{chrom}.tsv"
output_bgz = f"{output_file}.gz"
output_tbi = f"{output_bgz}.tbi"

# skip finished chromosomes, so a resubmitted array only redoes failed tasks; clear partial output
if os.path.exists(output_tbi):
    print(f"{output_tbi} already exists; skipping chromosome {chrom}")
    sys.exit(0)
for f in (output_file, output_bgz):
    if os.path.exists(f):
        os.remove(f)

df = ds.dataset(bank_dir, format="parquet").to_table(
    columns=["variant_id", "pos", "ref", "alt", "gene_id", "pip", "pred_prob"]).to_pandas()
assert len(df) > 0, f"no predictions in {bank_dir}"
assert df["variant_id"].str.startswith(f"chr{chrom}:").all(), "variant from another chromosome"

df = df.sort_values(["pos", "gene_id", "variant_id"], kind="mergesort").reset_index(drop=True)
out = df.rename(columns={"pos": "POS", "variant_id": "ID", "ref": "REF", "alt": "ALT",
                         "gene_id": "GENE_ID", "pip": "PIP", "pred_prob": "PRED_PROBABILITY"})
out.insert(0, "CHROM", chrom)
out = out[["CHROM", "POS", "ID", "REF", "ALT", "GENE_ID", "PIP", "PRED_PROBABILITY"]]
with open(output_file, "w") as f:
    f.write("#" + "\t".join(out.columns) + "\n")
    out.to_csv(f, sep="\t", index=False, header=False)
n_rows = len(out)
del df, out

# bgzip + tabix index: sequence column 1, position columns 2-2, '#' header line skipped
pysam.tabix_index(output_file, seq_col=0, start_col=1, end_col=1, meta_char="#", force=True)
n_back = sum(1 for _ in pysam.TabixFile(output_bgz).fetch(chrom))
assert n_back == n_rows, f"read back {n_back} rows, wrote {n_rows}"
print(f"{output_bgz}: {n_rows:,} rows")

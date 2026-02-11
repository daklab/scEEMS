# Step 7: Aggregate Predictions

Combine per-gene predictions into chromosome-level Parquet files.

## Overview

Merges individual gene prediction TSV files from Step 6 into a single Parquet dataset partitioned by chromosome, enabling efficient access for downstream analyses.

## Scripts

| Script | Description |
|--------|-------------|
| `aggregate_predictions.py` | Merge predictions into chromosome-level Parquet |
| `run_pipeline.sh` | Run aggregation for all cohorts |

## Inputs

- Per-gene prediction TSV files from Step 6

## Outputs

- `{data_dir}/training_data/{cohort}/predictions_parquet_catboost/predictions.parquet/` (partitioned by chromosome)

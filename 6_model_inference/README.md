# Step 6: Model Inference

Score all variants genome-wide using trained CatBoost models.

## Overview

For each gene, loads the LOCO-matched trained model (from Step 5) and predicts the probability that each variant in the gene's cis window is a causal eQTL variant.

## Scripts

| Script | Description |
|--------|-------------|
| `model_inference.py` | Load trained model and predict on all variant-gene pairs |
| `run_pipeline.sh` | Run inference for all genes |

## Inputs

- Trained models from Step 5
- Per-gene featurized variant data from Step 3

## Outputs

- `{data_dir}/training_data/{cohort}/predictions_catboost/{gene_id}_predictions.tsv`

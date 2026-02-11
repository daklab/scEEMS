# Step 4: Create Training Data

Sample positive and negative training examples from featurized gene datasets.

## Overview

Creates balanced training datasets by:
1. Identifying **positive** variants with high PIPs (above threshold) in credible sets
2. Sampling **negative** variants with low PIPs, stratified by variant type (SNP, insertion, deletion) to maintain the same ratio as positives
3. Merging variant features from Step 3

## Sampling Strategy

- **Train set**: PIP > 0.1 for positives, PIP < 0.01 for negatives
- **Test set**: PIP > 0.9 for positives, PIP < 0.01 for negatives
- **Negative-to-positive ratio (NPR)**: 10:1 (configurable)
- **Stratification**: Maintains SNP/indel ratio between positive and negative sets

## Scripts

| Script | Description |
|--------|-------------|
| `create_training_datasets.py` | Stratified sampling of variant-gene pairs |
| `data_params.yaml` | PIP thresholds for train/test splits |
| `run_pipeline.sh` | Run sampling for all chromosomes |

## Inputs

- Featurized gene datasets from Step 3
- PIP data from Step 1

## Outputs

- `{data_dir}/training_data/{cohort}/training_data/{split}_NPR_{N}_PIP_{pos}_{neg}/annotated_data_{cohort}_chr{N}.parquet`

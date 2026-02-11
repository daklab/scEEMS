# Step 5: Model Training

Train CatBoost classifiers for eQTL variant prediction using Leave-One-Chromosome-Out (LOCO) cross-validation.

## Overview

Trains the **conservative weighted** CatBoost model variant, which is used for final genome-wide inference. Key design choices:

- **LOCO CV**: One model per held-out chromosome (22 total per cell type)
- **Conservative regularization**: Depth 5, L2 regularization 5.0, bagging temperature 1.0
- **PIP-weighted samples**: Positive examples weighted by their PIP values
- **Feature weighting**: ChromBPNet, TF binding, and Enformer diff features weighted 10x

## Model Architecture

- **Algorithm**: CatBoost (gradient boosted decision trees)
- **Loss**: Logloss (binary classification)
- **Features**: ~200 features from variant annotations, cell-type epigenomics, TF binding, and genomic context

## Scripts

| Script | Description |
|--------|-------------|
| `train_model.py` | Train single CatBoost model for one held-out chromosome |
| `data_params.yaml` | PIP thresholds for train/test splits |
| `model_config.yml` | Model hyperparameters and feature weight configuration |
| `run_pipeline.sh` | Train models for all 22 chromosomes |

## Inputs

- Training data from Step 4
- Gene loss-of-function scores (Excel file)
- gnomAD MAF data
- Feature column dictionary (pickle)

## Outputs

- `{data_dir}/training_data/{cohort}/model_results/model_standard_subset_conservative_weighted_chr_{chr}_NPR_10.joblib`
- Feature importance CSVs
- Test set predictions

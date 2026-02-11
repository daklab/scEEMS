# Step 8: SHAP Analysis

Model interpretability using SHAP (SHapley Additive exPlanations) values.

## Overview

Computes TreeSHAP explanations for high-confidence predictions (pred_prob > 0.95), split by genomic context:
- **Enhancer variants**: >10kb from TSS
- **Promoter variants**: <=10kb from TSS

## Scripts

| Script | Description |
|--------|-------------|
| `shap_analysis.py` | Compute SHAP values per chromosome |
| `aggregate_shap_scores.py` | Combine per-chromosome SHAP to Parquet |
| `summary_cell_types_shap.py` | Cross-cell-type SHAP comparison |
| `run_pipeline.sh` | Run full SHAP pipeline |

## Outputs

- Per-chromosome SHAP CSVs (enhancer/promoter)
- Aggregated SHAP Parquet files
- Cross-cell-type SHAP summary TSV

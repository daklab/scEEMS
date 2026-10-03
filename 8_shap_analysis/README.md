# Step 8: SHAP Analysis

Which feature categories drive the scEEMS predictions, for promoter-like and enhancer-like variants.

## Overview

`shap_analysis.py` computes TreeSHAP values of the scEEMS (`weighted_full`) model for the model's
predicted eQTLs: variant-gene pairs with `pred_prob` above the cell type's tau*, the threshold selected by
the S-LDSC analysis of step 9 (`tau_star.json`). Each chromosome is explained by the LOCO model that held
it out. The |SHAP| values of each variant are summed within 11 categories that partition the 4,840
features:

| Category | Features |
|---|---|
| `ABC` | Activity-by-Contact scores |
| `CRE` | cell type CRE annotations (H3K27ac, H3K4me1, ATAC-seq) |
| `chromBPNet` | BPNet and ChromBPNet scores |
| `Enformer` | Enformer scores |
| `TF` | composite transcription factor scores |
| `abs_gpn` | GPN-STAR score |
| `gene_lof` | GeneBayes gene constraint |
| `variant_type` | length difference and variant type indicators |
| `baseline` | baseline annotations |
| `distance` | log distance to the TSS |
| `gnomad_MAF` | gnomAD minor allele frequency |

Variants within 10 kb of the TSS are promoter-like, the others enhancer-like. `aggregate_shap.py` averages
the summed |SHAP| over the predicted eQTLs of each cell type and region, and `summary_to_wide.py` reshapes
that table for plotting.

## Scripts

| Script | Description |
|--------|-------------|
| `shap_analysis.py` | TreeSHAP for one cell type and chromosome |
| `aggregate_shap.py` | Mean summed \|SHAP\| per category, cell type and region, and its share of the total |
| `summary_to_wide.py` | One row per cell type and region, one column per category |
| `run_shap.sh` | SLURM array, one task per chromosome |
| `run_pipeline.sh` | Everything on one machine |

## Running

```bash
cd 8_shap_analysis
sbatch --export=ALL,cohort=Mic_mega_eQTL run_shap.sh     # each cell type
python aggregate_shap.py                                # when all have finished
python summary_to_wide.py
```

The largest cell types need up to about 100 GB and a few hours for chromosome 1.

## Inputs

- the `weighted_full` models and `feature_cols.pkl` of step 5
- the `weighted_full` predictions of step 7 (`predictions_parquet/weighted_full/predictions.parquet`)
- the per-gene `all_variants` tables of step 3, which are not part of the data release
- `{aggregate_dir}/tau_star.json` from step 9

## Outputs

- `{output_dir}/{cohort}/shap/shap_category_chr{N}.parquet`: one row per predicted eQTL with
  `variant_id, gene_id, distance_TSS, region, tau_star` and the summed |SHAP| of each category
- `{aggregate_dir}/shap/all_cohorts_shap_summary.tsv`: `cohort, region, category, mean_abs_shap,
  share, n_variants`; `share` is the category's fraction of the total within the cell type and region
- `{aggregate_dir}/shap/all_cohorts_shap_summary_wide.tsv`

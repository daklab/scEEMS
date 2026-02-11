# Step 3: Cell Featurization

Create per-gene feature matrices combining variant annotations with cell-type-specific scores.

## Overview

For each gene, this step builds a complete feature matrix by combining:
- Variant-level annotations from Step 2 (Enformer scores, BED annotations, baseline)
- ABC (Activity-By-Contact) enhancer-gene link scores per cell type
- ChromBPNet cell-type-specific variant effect scores
- Aggregated TF binding score summaries (max, min, abs_max per cell type)
- Distance to TSS features

## Scripts

| Script | Description |
|--------|-------------|
| `create_gene_datasets.py` | Build annotated feature matrix for a single gene |
| `data_params.yaml` | Configuration for cell types and metrics |
| `run_pipeline.sh` | Run featurization for all genes in a cohort |

## Inputs

- Annotated variants from Step 2
- PIP data from Step 1
- ABC enhancer-gene predictions
- ChromBPNet variant effect scores
- Gene-to-name mapping file

## Outputs

- `{data_dir}/training_data/{cohort}/all_variants/{gene_id}/annotated_data_{cohort}_{chr}.parquet`

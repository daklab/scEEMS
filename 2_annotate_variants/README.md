# Step 2: Annotate Variants

Annotate the unique variant list with genomic features used as model inputs.

## Overview

This step integrates multiple sources of variant-level annotations:

1. **Enformer variant effect predictions** - Deep learning model predictions of variant effects on gene expression (Avsec et al., 2021)
2. **Brain cell-type epigenomic annotations** - ATAC-seq peaks, promoter/enhancer regions from sorted brain cell types (NeuN+, PU.1+, LHX2+, Olig2+)
3. **Baseline genomic annotations** - Standard functional annotations (e.g., coding regions, conserved elements)
4. **Transcription factor binding predictions** - Cell-type-specific TF binding based on Enformer scores intersected with epigenomic peaks

## Pre-computed Scores Required

The following pre-computed scores must be available before running this step:

### Enformer Variant Effect Predictions
- Source: [Avsec et al., 2021](https://doi.org/10.1038/s41592-021-01252-x)
- Format: Parquet files with diff_32 scores per Enformer output track
- Expected path: `{enformer_dir}/enformer_tensorflow_chr{N}.parquet`

### ChromBPNet Cell-Type-Specific Scores
- Source: [Pampari et al., 2023](https://doi.org/10.1101/2023.11.27.568824)
- Used in Step 3 (cell featurization), not directly in this step
- Scores for: microglia, astrocyte, oligodendrocyte, neuron

### BPNet Scores
- Transcription factor binding predictions from BPNet models

### GPN-MSA Conservation Scores
- Multi-species alignment-based conservation scores

### Brain Epigenomic BED Files
- ENCODE ATAC-seq peaks for brain cell types
- Promoter and enhancer annotations
- Expected location: `{bed_annotations_dir}/`

## Scripts

| Script | Description |
|--------|-------------|
| `annotate_variants.py` | Main annotation script merging all data sources |
| `run_pipeline.sh` | Run annotation for all chromosomes |

## Inputs

- `{data_dir}/susie_vars_pips/variant_list/variant_list_chr{N}.parquet` (from Step 1)
- Pre-computed Enformer scores (Parquet)
- Brain cell-type BED files
- Baseline annotation BED files
- TF binding metadata files

## Outputs

- `{data_dir}/susie_vars_pips/variant_list/annotated_variants/annotated_variants_chr{N}.parquet`

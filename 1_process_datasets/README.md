# Step 1: Process Datasets

Process raw SuSiE fine-mapping RDS files into model-ready variant tables.

## Overview

This step documents how internal raw SuSiE fine-mapping RDS files were converted into:
- variant-level PIP tables
- per-cohort parquet datasets
- merged unique variant lists used downstream

Important: the raw RDS inputs used here are **not publicly distributed** in this repository.
Public users should download preprocessed train/test data from Synapse (see root `README.md` and `download_synapse_data.py`).

## Data Source

Fine-mapping results come from the FunGen-xQTL Analysis Freeze, which includes SuSiE cis-eQTL fine-mapping across multiple brain cell types from snRNA-seq data.

## Scripts

| Script | Description |
|--------|-------------|
| `download_data.py` | Public guardrail script explaining that raw Step 1 inputs are not distributed |
| `create_job_list.py` | Parse downloaded filenames to extract unique gene IDs |
| `get_vars_pips.R` | Extract PIPs and top loci from SuSiE RDS files per gene |
| `create_parquet_files.py` | Convert per-gene CSV files to partitioned Parquet |
| `create_unique_variant_list.py` | Create deduplicated variant list per chromosome |
| `run_pipeline.sh` | Run all scripts in order |

## Inputs

- Internal/private raw RDS files under `{data_dir}/release_04_2024/`
- `Fungen_xQTL.<gene_id>.cis_results_db.export.rds`
- `Fungen_xQTL.<gene_id>.cis_results_db.export_sumstats.rds`

## Outputs

- `{data_dir}/release_04_2024/` - Raw RDS files
- `{data_dir}/susie_vars_pips/{cohort}/PIP_all_parquet/` - All variant PIPs (Parquet)
- `{data_dir}/susie_vars_pips/{cohort}/PIP_top_parquet/` - Top loci PIPs (Parquet)
- `{data_dir}/susie_vars_pips/variant_list/variant_list_chr{N}.parquet` - Unique variant lists

## Cell Types

- `Ast_mega_eQTL` - Astrocytes
- `Exc_mega_eQTL` - Excitatory neurons
- `Inh_mega_eQTL` - Inhibitory neurons
- `Mic_mega_eQTL` - Microglia
- `Oli_mega_eQTL` - Oligodendrocytes
- `OPC_mega_eQTL` - Oligodendrocyte precursor cells

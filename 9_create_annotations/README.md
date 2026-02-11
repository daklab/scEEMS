# Step 9: Create LDSC Annotations and Run Partitioned Heritability

This step creates LDSC-compatible annotations from scEEMS predictions and runs LD score regression to quantify enrichment of Alzheimer's disease heritability in predicted eQTL annotations.

## Overview

1. **Create annotations** (`make_annotations*.py`): Merge scEEMS prediction probabilities with baseline LD annotations at multiple probability thresholds, producing per-chromosome annotation files.
2. **Compute LD scores** (`template_ldscore*.sh`): Use PolyFun's `compute_ldscores.py` to compute LD scores for the new annotations.
3. **Run LDSC regression** (`ldscore_regression/ldscore_template*.py`): Run partitioned heritability analysis using PolyFun's `ldsc.py` to test enrichment of each annotation.
4. **Aggregate results** (`ldscore_regression/aggregate_data.py`): Combine results across model variants and cell types.

## Model Variants

| Variant | Annotation script | Predictions directory | Description |
|---------|------------------|-----------------------|-------------|
| Revised (no control) | `make_annotations.py` | `predictions_parquet_catboost_revised_nocontrol` | Primary model used in paper |
| Original (weighted) | `make_annotations_original.py` | `predictions_parquet_catboost` | Original weighted training |
| Unweighted | `make_annotations_original_unweighted.py` | `predictions_parquet_catboost_unweighted` | Unweighted training baseline |

## Prerequisites

- Completed Steps 1-7 (predictions available)
- Baseline LD annotations (`paths.baseline_annot_dir` in `config.yaml`)
- ABC gene mapping data (`paths.abc_data_dir` in `config.yaml`)
- PLINK reference panel (`paths.plink_ref_dir`)
- PolyFun installation (`paths.ldsc_dir`)
- GWAS summary statistics (`paths.sumstats_file`)
- LD regression weights (`paths.weights_dir`)

## Usage

See `run_pipeline.sh` for the full execution sequence, or run individual steps:

```bash
# 1. Create annotations for one cell type (all chromosomes via SLURM)
sbatch --export=cell_type=Mic_mega_eQTL run_jobs.sh

# 2. Compute LD scores (after annotations complete)
sbatch --export=cell_type=Mic_mega_eQTL template_ldscore.sh

# 3. Run LDSC regression
cd ldscore_regression
sbatch run_jobs.sh

# 4. Aggregate results
python aggregate_data.py
```

## Output

- `MLxQTL_chr{chr}.annot.gz`: Per-chromosome annotation files with binary columns at each threshold
- `MLxQTL_chr{chr}.l2.M`: Annotation count summaries
- `MLxQTL_chr{chr}.l2.ldscore.parquet`: LD scores for annotations
- `ldscore_pareto_data.tsv`: Combined enrichment results across model variants
- `ldscore_pareto_data_all_cell_types.tsv`: Enrichment results across all cell types

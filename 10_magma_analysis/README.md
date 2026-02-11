# Step 10: MAGMA Gene-Set Enrichment Analysis

This step runs MAGMA (Multi-marker Analysis of GenoMic Annotation) gene-set enrichment analysis using scEEMS predictions to define gene-variant associations.

## Overview

1. **Prepare summary statistics** (`run_magma_analysis/make_MAGMA_sumstats*.py`): Format GWAS summary statistics for MAGMA input.
2. **Create MAGMA annotation files** (`make_magma_files_revised.py`, `make_magma_files_knn.py`): Map variants to genes based on scEEMS predictions, ABC enhancer scores, and PLACseq interactions.
3. **Run MAGMA** (`run_magma_analysis/run_magma*.sh`): Execute MAGMA gene-based association tests across multiple GWAS datasets and populations.
4. **Extract significant genes** (`extract_significant_MAGMA_genes*.py`): Identify Bonferroni-significant genes and compare with baseline MAGMA results.

## MAGMA Variants

| Variant | Script | Description |
|---------|--------|-------------|
| Revised | `make_magma_files_revised.py` | Uses cell-type-specific ABC + PLACseq filtering with per-cell-type probability thresholds |
| KNN | `make_magma_files_knn.py` | Selects N closest variants to gene TSS, where N = number of revised predictions per gene |

## GWAS Datasets

| Dataset | Population | Script |
|---------|-----------|--------|
| Bellenguez et al. 2022 | EUR | `run_magma.sh` |
| Kunkle et al. 2019 | EUR | `run_magma_kunkle.sh` |
| Kunkle et al. 2019 | AFR | `run_magma_kunkle_AFR.sh` |
| ADSP | AFR | `run_magma_ADSP_AFR.sh` |
| ADSP | AMR | `run_magma_ADSP_AMR.sh` |
| ADGC | AFR | `run_magma_ADGC_AFR.sh` |
| ADGC | AMR | `run_magma_ADGC_AMR.sh` |
| ADGC | EAS | `run_magma_ADGC_EAS.sh` |

## Prerequisites

- Completed Steps 1-7 (predictions available)
- MAGMA binary installed (`paths.magma_dir` in `config.yaml`)
- PLINK reference panels (`paths.plink_ref_dir`)
- GWAS summary statistics
- Baseline LD annotations (`paths.baseline_annot_dir`)
- ABC enhancer data (`paths.abc_data_dir`)
- PLACseq interaction data (`paths.placseq_file`)
- Gene info file (`paths.gene_info_file`)

## Usage

See `run_pipeline.sh` for the full execution sequence, or run individual steps:

```bash
# 1. Prepare summary statistics
cd run_magma_analysis
python make_MAGMA_sumstats.py
python make_MAGMA_sumstats_v2.py

# 2. Create MAGMA annotation files
cd ..
sbatch --export=cell_type=Mic_mega_eQTL run_jobs.sh

# 3. Run MAGMA (after annotation files complete)
cd run_magma_analysis
sbatch --export=cell_type=Mic_mega_eQTL run_magma.sh

# 4. Extract significant genes
cd ..
python extract_significant_MAGMA_genes.py
```

## Output

- `MAGMA/{cell_type}/prediction/`: MAGMA annotation files based on prediction probabilities
- `MAGMA/{cell_type}/pip/`: MAGMA annotation files based on PIP scores
- `MAGMA/{cell_type}/variant_gene_lists/`: Variant-gene association tables
- `MAGMA/{cell_type}/MAGMA_output/`: MAGMA gene-based test results
- `MAGMA_significant_genes/`: Bonferroni-significant gene lists with/without overlap to baseline

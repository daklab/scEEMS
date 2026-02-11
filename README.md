# scEEMS: Machine Learning-Based Prediction of Cell-Type-Resolved Brain eQTLs

**Version**: 1.0.0

Code for the paper: *Machine Learning-Based Prediction of Cell-type Resolved Brain eQTLs Enhances Discovery of Variants Explaining Alzheimer's Disease Heritability*

Preprint (medRxiv): https://doi.org/10.64898/2025.12.03.25341562

Authors: Chirag M Lakhani, Giacomo Cavalca, Anjing Liu, Rohan Nidumbur, Ru Feng, Towfique Raj, Philip De Jager, The Alzheimer's Disease Functional Genomics Consortium, Gao Wang, David A. Knowles

## Overview

scEEMS (single-cell Expression-mediated Effect Modeling System) is a CatBoost-based framework that predicts whether genetic variants are expression quantitative trait loci (eQTLs) in specific brain cell types. The model integrates variant-level genomic annotations (Enformer, ChromBPNet, GPN-MSA, conservation scores) with cell-type-specific regulatory features (ABC enhancer scores, transcription factor binding) to generate per-variant, per-gene, per-cell-type eQTL probabilities.

These predictions are used for:
- **LDSC partitioned heritability** analysis to quantify enrichment of Alzheimer's disease heritability
- **MAGMA gene-set enrichment** to identify disease-relevant gene sets
- **Informed fine-mapping priors** to improve SuSiE credible set resolution

## System Requirements

### Software Dependencies

- **Operating system (tested)**: Ubuntu 22.04.5 LTS (GNU/Linux 5.15; NYGC cluster login node)
- **Python**: 3.9.18 (conda)
- **R**: 4.3.3
- **Core Python packages**: catboost, numpy, pandas, scipy, scikit-learn, pyarrow, dask, shap, pybedtools, pysam, pyranges, synapseclient, pyyaml
- **Core R packages**: susieR, tidyverse, readr, yaml, pecotmr
- **External tools** (configured via `config.yaml`):
  - [bedtools](https://bedtools.readthedocs.io/) (v2.31.1)
  - [MAGMA](https://ctg.cncr.nl/software/magma) (gene-set analysis binary)
  - [PolyFun/LDSC](https://github.com/omerweisssbrod/polyfun) (LD score regression)
  - PLINK reference panel (for LD score computation)
- **External data** (configured via `config.yaml`): Enformer/ChromBPNet/GPN-MSA variant effect predictions, ABC enhancer scores, baseline LD annotations, gnomAD MAF data, TF binding files, target lists
- **Synapse data for model training**: download required training/test inputs from the Synapse project and use those paths in `config.yaml` (see Additional Information).

### Versions Tested

The pipeline was run on the NYGC cluster with Ubuntu 22.04.5 LTS, Python 3.9.18, and R 4.3.3. The exact conda environment is recorded in `conda_environment_full.txt`; a minimal dependency list is provided in `environment.yml` and `requirements.txt`.

### Hardware

- **No GPU required** for the main pipeline.
- **Full pipeline**: requires a high-memory multi-core environment. Typical runs use **5–20 CPU cores** and **30–200 GB RAM** depending on step.
- **Small-scale/demo runs**: can be executed on a standard desktop (8–32 GB RAM) by restricting to a single cell type and chromosome.

## Pipeline

| Step | Directory | Description |
|------|-----------|-------------|
| 1 | `1_process_datasets/` | Internal preprocessing of raw fine-mapping RDS files (not publicly distributed) |
| 2 | `2_annotate_variants/` | Annotate variants with genomic features (Enformer, ChromBPNet, GPN-MSA) |
| 3 | `3_cell_featurization/` | Create per-gene feature matrices with cell-type-specific annotations |
| 4 | `4_create_training_data/` | Sample positive/negative training examples using PIP thresholds |
| 5 | `5_model_training/` | Train CatBoost classifiers with leave-one-chromosome-out CV |
| 6 | `6_model_inference/` | Score all variants genome-wide using trained models |
| 7 | `7_aggregate_predictions/` | Combine per-gene predictions into chromosome-level files |
| 8 | `8_shap_analysis/` | Compute SHAP values for model interpretability |
| 9 | `9_create_annotations/` | Create LDSC annotations and run partitioned heritability analysis |
| 10 | `10_magma_analysis/` | Run MAGMA gene-set enrichment analysis |
| 11 | `11_fine_mapping_new_priors/` | Re-run SuSiE fine-mapping with scEEMS-informed priors |

Each step directory contains a `README.md` with detailed documentation and a `run_pipeline.sh` example script.
For provenance of how preprocessing files were created (and which upstream inputs are non-public), see `DATA_PROCESSING_NOTES.md`.

## Installation Guide

### 1. Clone and configure

```bash
git clone https://github.com/daklab/scEEMS.git
cd scEEMS
cp config.yaml.example config.yaml
# Edit config.yaml with your paths and credentials
```

### 2. Install dependencies (conda recommended)

```bash
conda env create -f environment.yml
conda activate scEEMS
```

Alternatively, for pip users:
```bash
pip install -r requirements.txt
```

The file `conda_environment_full.txt` records the exact conda environment used during development.

**R packages** (for Steps 1 and 11):
```r
install.packages(c("susieR", "pecotmr", "tidyverse", "readr", "yaml"))
```

**Typical install time**: ~15–30 minutes on a normal desktop (longer if conda needs to download many packages).

## Public Data Download

Publicly downloadable Synapse resources in this repository:

- `model_training` (Project: `syn72248754`): train/test parquet data and supporting files (`columns_dict`, `gene_lof`, `gnomad_MAF`)
- `predictions` (Folder: `syn71338354`): per-cell-type prediction outputs

Download with the helper script:

```bash
cp config.yaml.example config.yaml
python download_synapse_data.py --resource all
```

Or selectively:

```bash
python download_synapse_data.py --resource model_training
python download_synapse_data.py --resource predictions
```

By default files are downloaded to `{paths.data_dir}/synapse_public/`.

## Minimal Training (Chromosome 2 Only)

This public demo runs **one cell type** (Microglia) and **one chromosome** (**chr2 only**) for model training (Step 5).

### 1. Download chr2 demo files

Use `5_model_training/README_minimal_chr2_demo.md` for:
- exact required files
- exact download URLs
- parquet-directory setup (`part.*.parquet` files)

### 2. Configure paths

```bash
cp config.yaml.example config.yaml
# Set:
# - paths.data_dir to your demo root (with chr2 train/test parquet directories)
# - paths.gnomad_maf_dir to where gnomad_MAF_chr2.tsv is stored
# - paths.columns_dict_file to columns_dict.pkl
```

### 3. Run the chr2 demo

```bash
cd 5_model_training
bash run_minimal_training.sh \
  Mic_mega_eQTL \
  2 \
  "<AUX_ROOT>/41588_2024_1820_MOESM4_ESM.xlsx"
```

Equivalent direct command:

```bash
cd 5_model_training
python train_model.py Mic_mega_eQTL 2 \
  --gene_lof_file "<AUX_ROOT>/41588_2024_1820_MOESM4_ESM.xlsx" \
  --yaml_path data_params.yaml \
  --single_chromosome_demo
```

### Expected output

- Model file: `{data_dir}/training_data/Mic_mega_eQTL/model_results/model_standard_subset_conservative_weighted_chr_chr2_NPR_10.joblib`
- Feature importance CSVs
- Test set predictions

### Expected runtime

- **Normal desktop** (8–16 GB RAM, 4–8 cores): ~1–3 hours for a single chromosome (varies by data size).
- **Cluster**: a typical run for this step uses about 5 CPU cores, 50 GB RAM, and up to 5 hours.

For full multi-chromosome training, download the full Synapse `model_training` resource (`syn72248754`) and use the standard training workflow.

## Instructions For Use

### Run the pipeline

Each step can be run independently. See the `README.md` in each directory for specific instructions.

Note: Step 1 uses internal raw inputs that are not publicly distributed.
For public workflows, start from Synapse `model_training` / `predictions` data and run downstream steps.

```bash
cd 5_model_training && bash run_pipeline.sh
cd ../6_model_inference && bash run_pipeline.sh
# ... and so on
```

### Cluster/HPC execution

If you are running on a SLURM cluster, use the `run_jobs.sh` scripts provided in step directories where available and adapt resource requests (paths, partitions, and time limits) to your environment.

## (Optional) Reproduction Instructions

To reproduce all manuscript results exactly, run the full pipeline in order (Steps 1–11) using the same external datasets and annotations described in each step’s `README.md`.

Note: Step 1 raw fine-mapping RDS inputs are not publicly distributed in this repository. Public users can run training/inference workflows from the Synapse `model_training` and `predictions` resources.

## Configuration

All paths and credentials are managed through `config.yaml`. Copy `config.yaml.example` and fill in your values:

```yaml
paths:
  data_dir: "/path/to/data"
  output_dir: "/path/to/output"
  fine_mapping_dir: "/path/to/fine_mapping"
  scratch_dir: "/tmp"

credentials:
  synapse_token: ""

data_sources:
  model_training_synapse_id: "syn72248754"
  predictions_synapse_id: "syn71338354"
```

## Pre-computed Scores

Step 2 requires pre-computed variant effect predictions. See `2_annotate_variants/README.md` for details on obtaining:

- **Enformer** variant effect predictions (Avsec et al., 2021)
- **ChromBPNet** cell-type-specific accessibility predictions (Pampari et al., 2023)
- **GPN-MSA** conservation scores (Benegas et al., 2023)
- **Baseline LD annotations** (Gazal et al., 2017)

## Cell Types

The pipeline supports the following brain cell types:

| Abbreviation | Cell Type |
|---|---|
| Ast | Astrocytes |
| Exc | Excitatory neurons |
| Inh | Inhibitory neurons |
| Mic | Microglia |
| Oli | Oligodendrocytes |
| OPC | Oligodendrocyte precursor cells |

## Additional Information

- **Data access**: Publicly distributed inputs are in Synapse `model_training` (`syn72248754`) and `predictions` (`syn71338354`). Step 1 raw fine-mapping RDS inputs are not publicly distributed in this repo. Step 2 requires large annotation datasets that are not distributed with this repo.
- **Synapse data access**: Data access requires a [Synapse](https://www.synapse.org/) account and may require approval for controlled-access datasets.
- **Synapse model training data** (Project: [syn72248754](https://www.synapse.org/Synapse:syn72248754)):

  **Training data** (per cell type):

  | Cell Type | Synapse ID |
  |-----------|-----------|
  | Astrocytes (Ast) | syn72248805 |
  | Excitatory neurons (Exc) | syn72248803 |
  | Inhibitory neurons (Inh) | syn72248804 |
  | Microglia (Mic) | syn72248802 |
  | OPC | syn72248807 |
  | Oligodendrocytes (Oli) | syn72248806 |

  **Test data** (per cell type):

  | Cell Type | Synapse ID |
  |-----------|-----------|
  | Astrocytes (Ast) | syn72248942 |
  | Excitatory neurons (Exc) | syn72248940 |
  | Inhibitory neurons (Inh) | syn72248941 |
  | Microglia (Mic) | syn72248939 |
  | OPC | syn72248944 |
  | Oligodendrocytes (Oli) | syn72248943 |

  **Supporting files**:

  | File | Synapse ID |
  |------|-----------|
  | README_model_training.txt | syn72249887 |
  | columns_dict | syn72248797 |
  | gene_lof | syn72248800 |
  | gnomad_MAF | syn72248799 |

  **Predictions root folder**: syn71338354

  **Predictions by cell type**:

  | Resource | Synapse ID |
  |------|-----------|
  | Ast_mega_eQTL | syn71338358 |
  | Exc_mega_eQTL | syn71338356 |
  | Inh_mega_eQTL | syn71338357 |
  | Mic_mega_eQTL | syn71338355 |
  | OPC_mega_eQTL | syn71338360 |
  | Oli_mega_eQTL | syn71338359 |
  | README_predictions.txt | syn71338459 |
- **Storage**: Full pipeline outputs are large (multi-terabyte scale depending on cohorts and annotations). Ensure adequate disk space in `paths.output_dir`.
- **Compute**: The pipeline is designed for HPC. For small tests or demos, limit to a single cohort and chromosome to reduce runtime and memory usage.

### Download From Synapse

1. Install and authenticate the Synapse client.
2. Download public `model_training` and/or `predictions` resources by Synapse ID.

Project URLs:
```text
https://www.synapse.org/Synapse:syn72248754
https://www.synapse.org/Synapse:syn71338354
```

Example commands:
```bash
pip install synapseclient
synapse login

# Download model_training project
synapse get syn72248754

# Download predictions folder
synapse get syn71338354

# Or download specific model_training folders
synapse get syn72248797  # columns_dict
synapse get syn72248800  # gene_lof
synapse get syn72248799  # gnomad_MAF
synapse get syn72248798  # train
synapse get syn72248801  # test
```

Raw Step 1 fine-mapping RDS files are not publicly distributed in this repository.

## Citation

If you use this code, please cite:

```bibtex
@article{lakhani2025sceems,
  title={Machine Learning-Based Prediction of Cell-type Resolved Brain eQTLs Enhances Discovery of Variants Explaining Alzheimer's Disease Heritability},
  author={Lakhani, Chirag M and Cavalca, Giacomo and Liu, Anjing and Nidumbur, Rohan and Feng, Ru and Raj, Towfique and De Jager, Philip and The Alzheimer's Disease Functional Genomics Consortium and Wang, Gao and Knowles, David A.},
  journal={medRxiv},
  year={2025},
  doi={10.64898/2025.12.03.25341562},
  url={https://doi.org/10.64898/2025.12.03.25341562}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

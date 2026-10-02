# scEEMS: Machine Learning-Based Prediction of Cell-Type-Resolved Brain eQTLs

**Version**: 2.0.0 (revised manuscript)

Code for the paper: *Machine Learning-Based Prediction of Cell-type Resolved Brain eQTLs Enhances Discovery
of Variants Explaining Alzheimer's Disease Heritability*

Preprint (medRxiv): https://doi.org/10.64898/2025.12.03.25341562

Authors: Chirag M Lakhani, Giacomo Cavalca, Anjing Liu, Rohan Nidumbur, Ru Feng, Towfique Raj, Philip De
Jager, The Alzheimer's Disease Functional Genomics Consortium, Gao Wang, David A. Knowles

## Overview

single-cell Enhanced Expression Modifier Scores (scEEMS) are CatBoost models that predict, for six brain
cell types, the probability that a variant is a causal eQTL for a gene. Each variant-gene pair is described
by 4,840 features: deep learning variant effect predictions (Enformer, BPNet, ChromBPNet and composite
transcription factor scores), the GPN-STAR DNA language model score, Activity-by-Contact (ABC) scores, cell
type cis-regulatory element annotations, baseline genomic annotations, distance to the TSS, variant type,
gnomAD allele frequency and GeneBayes gene constraint. The models are trained on fine-mapped single-cell
eQTLs of six cell types (FunGen-xQTL) with leave-one-chromosome-out cross-validation.

The predictions are used to:
- partition Alzheimer's disease heritability with stratified LD score regression (S-LDSC);
- link variants to genes for eQTL-informed MAGMA gene analysis (eMAGMA), in European and non-European GWAS;
- serve as priors for eQTL fine-mapping, which is then colocalized with AD GWAS.

### Changes in the revision

- **GPN-STAR**: the absolute GPN-STAR log-likelihood ratio is a new feature. The scores of all 15.3 million
  SNVs are in the data release (`model_training/gpn_star/`); `10_magma_analysis/noneur_gpn_star/` shows how
  variants are scored with the GPN-STAR model.
- **Data-driven feature weights**: instead of a fixed tenfold weight on the DL-VEP features, one CatBoost
  feature weight per feature category is selected for each cell type by a multi-objective Bayesian search
  (step 5), on odd and even chromosomes separately so a held-out chromosome never informs its own weights.
- **Comparison models**: Unweighted (Full) (all weights 1) and Weighted (Restricted) (trained only on
  positives with PIP > 0.9 in a credible set) are trained alongside scEEMS, Weighted (Full).
- **Fine-mapping with five priors and colocalization** (step 11): uniform, EMS and the three models' priors,
  colocalized with the Bellenguez et al. AD GWAS; cross-cell-type sharing of credible sets (step 12).
- **Predicted-eQTL threshold**: the probability cut that defines predicted eQTLs (tau*) is chosen per cell
  type by S-LDSC (step 9) and used by MAGMA (step 10) and SHAP (step 8).

## Pipeline

| Step | Directory | Description | Runs from the data release |
|------|-----------|-------------|---|
| 1 | `1_process_datasets/` | SuSiE eQTL fine-mapping results (FunGen-xQTL) to variant tables | no: inputs not distributed |
| 2 | `2_annotate_variants/` | Variant annotations: Enformer, cell type epigenomics, baseline, TF scores | no: large external inputs |
| 3 | `3_cell_featurization/` | Per-gene feature tables (`all_variants`) with ABC, ChromBPNet and distance features | no: steps 1-2 |
| 4 | `4_create_training_data/` | Training and test sets (positives and matched negatives) | no: steps 1-3 (outputs are in the release) |
| 5 | `5_model_training/` | Feature-weight search, LOCO training of the three models, AUPRC | **yes** |
| 6 | `6_model_inference/` | Score every cis-variant of every gene | no: step 3 tables |
| 7 | `7_aggregate_predictions/` | Collect predictions; export or import the released TSVs | **yes** (import/export) |
| 8 | `8_shap_analysis/` | SHAP attribution by feature category | no: step 3 tables |
| 9 | `9_create_annotations/` | S-LDSC heritability of predicted eQTLs; tau* | partly: the scEEMS threshold sweep, tau* and the PIP > 0.10 comparison, with LDSC reference data |
| 10 | `10_magma_analysis/` | eMAGMA gene analysis, European and non-European GWAS | partly (see its README) |
| 11 | `11_finemap_coloc/` | eQTL fine-mapping with five priors; colocalization with AD GWAS | no: controlled-access genotypes (outputs are in the release) |
| 12 | `12_crosscell_coloc/` | Sharing of eQTL credible sets between cell types | no: step 11 fits |

Each step directory has a `README.md` (method, inputs, how to run, outputs) and SLURM scripts. Steps 1-4,
6, 8, 11 and 12 need inputs that are not in the data release (controlled-access genotype and expression
data, or intermediate feature tables of several terabytes); their code documents exactly how the released
data and results were produced.

## Data release

The data release is on Synapse in folder [syn69670587](https://www.synapse.org/Synapse:syn69670587):

| Folder | Content | Used by |
|---|---|---|
| `model_training/` | `train/`, `train_restricted/` and `test/` sets per cell type; `gpn_star/` (GPN-STAR scores); `feature_weights/` (selected weights); `gnomad_MAF/`, `gene_lof/`, `columns_dict/`; `model_features.tsv` | step 5 |
| `predictions/` | scEEMS predictions per cell type and chromosome (tabix-indexed TSV) | steps 9-10, your own analyses |
| `fine_mapping/` | Credible sets of the five fine-mapping priors per cell type and chromosome (tabix-indexed TSV) | your own analyses |

Each folder has a README describing its files and columns. Download with the helper script (Synapse account
required; listing works without one):

```bash
python download_synapse_data.py --dry-run                                   # what would be downloaded
python download_synapse_data.py --resource model_training --cell-type Mic   # microglia only, 3.2 GB
python download_synapse_data.py                                             # everything, ~66 GB
```

Files go to `paths.release_dir` with the release's folder layout, which is where the code looks for them.

## Installation

```bash
git clone https://github.com/daklab/scEEMS.git
cd scEEMS
conda env create -f environment.yml       # Python environment, steps 1-10
conda activate scEEMS
cp config.yaml.example config.yaml        # then set release_dir, output_dir and any inputs you need
```

`environment.yml` pins the versions used for the manuscript (installation takes about 10 minutes);
`conda_environment_full.txt` lists the full environment used to train and score the models. The
fine-mapping and colocalization steps (11-12) use R 4.5 in a second environment, `environment_r.yml`, plus
seven R packages that are not on conda, installed by `install_r_packages.R` (see the header of either file).

### System requirements

- Linux (tested on Ubuntu 22.04), Python 3.9 (`environment.yml`), R 4.5 for steps 11-12 (`environment_r.yml`)
- No GPU, except for the optional GPN-STAR scoring of non-European variants in step 10
- Step 5: 10 CPU cores and 12 GB (microglia) to 60 GB (excitatory neurons) of memory per held-out
  chromosome, 10-30 minutes each
- External tools for specific steps: PolyFun/LDSC (step 9), MAGMA v1.10 (step 10); see each step's README

## Configuration

All paths are set in `config.yaml` (gitignored), or in another file named by the `SCEEMS_CONFIG`
environment variable. Steps 5 and 7 need only `release_dir` and `output_dir`:

```yaml
paths:
  release_dir: "/path/to/scEEMS_data"     # the downloaded data release
  output_dir: "/path/to/output"           # everything the pipeline writes
  data_dir: "/path/to/data"               # steps 1-4 and their outputs (not needed for steps 5 and 7)
```

Path settings are templates that may refer to `{release_dir}`, `{output_dir}`, `{data_dir}`, other path
settings and `{cohort}`; anything not set in `config.yaml` takes the default in `shared/config.py`
(inputs from the data release, outputs under `output_dir`). `config.yaml.example` lists every setting.

## Quick start: reproduce a microglia model

```bash
python download_synapse_data.py --resource model_training --cell-type Mic
cd 5_model_training
python train_loco.py Mic_mega_eQTL 1
```

This trains the three microglia models with chromosome 1 held out (10-20 minutes with 10 CPU cores) and
scores the chromosome 1 test set: 242 variant-gene pairs, AUPRC 0.8015 for scEEMS (`weighted_full`), 0.7555
for Unweighted (Full) and 0.7499 for Weighted (Restricted). Models trained from the data release reproduce
the published models exactly (see `5_model_training/README.md` for the one exception).

## Cell types

| Abbreviation | Cell type |
|---|---|
| Ast | Astrocytes |
| Exc | Excitatory neurons |
| Inh | Inhibitory neurons |
| Mic | Microglia |
| Oli | Oligodendrocytes |
| OPC | Oligodendrocyte precursor cells |

The code refers to each cell type's data as `{Cell}_mega_eQTL` (e.g. `Mic_mega_eQTL`).

## Citation

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

This project is licensed under the MIT License; see [LICENSE](LICENSE).

# Step 4: Create Training Data

Build the training and test sets of step 5: positive variant-gene pairs from the fine-mapping results of
step 1 and, for each positive, matched negatives from the same gene, with their features from step 3.

**This step needs the outputs of steps 1-3, which are not distributed.** Its outputs are the
`model_training/train/`, `train_restricted/` and `test/` folders of the data release, so step 5 can be run
without it.

## Data sets

| Split | Positives | Used for |
|---|---|---|
| `train` | variants with PIP > 0.05 in a 95% credible set whose top variant has PIP > 0.10, and variants with PIP > 0.50 outside any credible set | Weighted (Full) and Unweighted (Full) models |
| `train_restricted` | variants with PIP > 0.90 in a 95% credible set | Weighted (Restricted) model |
| `test` | variants with PIP >= 0.90 | held-out evaluation (AUPRC) |

Positives come from the top loci of the cell type's MEGA analysis (`PIP_top`). For each positive, 10
negatives (`NPR`) are sampled from the same gene: variants with PIP < 0.01 and the same variant type (SNV,
insertion or deletion), with `random_state` 42, with replacement when the gene has too few. Only variants
with features (step 3) are used. Thresholds are in `data_params.yaml`.

A chromosome without positives gets no file (microglia has no chromosome 21 file for `test` and
`train_restricted`); step 5 treats it as a chromosome without data.

## Scripts

| Script | Description |
|--------|-------------|
| `create_training_datasets.py` | One split of one cell type and chromosome |
| `data_params.yaml` | PIP thresholds of the three splits |
| `run_pipeline.sh` | The three splits of every chromosome of one cell type |

## Running

```bash
bash 4_create_training_data/run_pipeline.sh Mic_mega_eQTL
# or one file:
python 4_create_training_data/create_training_datasets.py 22 train_restricted Mic_mega_eQTL 10
```

## Inputs

- `{susie_pips_dir}/PIP_top_parquet/` and `PIP_all_parquet/` of the cell type's MEGA analysis (step 1)
- `{variant_list_dir}/annotated_variants_just_variants/` (step 2)
- `{all_variants_dir}/<gene_id>/` (step 3)

## Outputs

`{training_sets_dir}/<split>_NPR_10_PIP_<positive>_<negative>/annotated_data_<cohort>_chr{N}.parquet`, i.e.
`train_NPR_10_PIP_0.1_0.01/`, `train_restricted_NPR_10_PIP_0.9_0.01/` and `test_NPR_10_PIP_0.9_0.01/`
(default `training_sets_dir`: `{data_dir}/training_data/{cohort}/training_data`). Each file has the
step 3 features of its variant-gene pairs plus `label` (1 positive, 0 negative).

To train step 5 on your own step 4 outputs instead of the data release, set in `config.yaml`:

```yaml
paths:
  train_dir: "{data_dir}/training_data/{cohort}/training_data/train_NPR_10_PIP_0.1_0.01"
  train_restricted_dir: "{data_dir}/training_data/{cohort}/training_data/train_restricted_NPR_10_PIP_0.9_0.01"
  test_dir: "{data_dir}/training_data/{cohort}/training_data/test_NPR_10_PIP_0.9_0.01"
```

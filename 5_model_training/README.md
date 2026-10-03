# Step 5: Model Training

Train the scEEMS CatBoost classifiers with leave-one-chromosome-out (LOCO) cross-validation and measure
their held-out AUPRC, or evaluate the published models on the same test sets without training. Everything
this step reads is in the `model_training/` folder of the data release, so it can be run without steps 1-4.

## Overview

For each cell type and each held-out chromosome, `train_loco.py` trains three models on the other 21
chromosomes:

| Model | Manuscript name | Training data | Feature-category weights |
|---|---|---|---|
| `weighted_full` | Weighted (Full), **scEEMS** | `train/` | selected weights |
| `unweighted_full` | Unweighted (Full) | `train/` | all 1 |
| `weighted_restricted` | Weighted (Restricted) | `train_restricted/` | selected weights |

and scores the held-out chromosome's test set (PIP > 0.9 positives, 10 matched negatives each).
`weighted_full` is the production model: step 6 uses it to score every variant, and its predictions are
the scEEMS predictions of the data release.

- **Features**: 4,840 per variant-gene pair: DL-VEP scores (Enformer, BPNet, ChromBPNet, composite TF
  scores), ABC scores, cell type CRE annotations, baseline annotations, log distance to the TSS, variant
  type, gnomAD MAF, GeneBayes gene constraint and the GPN-STAR score (`shared/featurize.py`).
- **CatBoost settings** (fixed): depth 5, 1,000 iterations, learning rate 0.03, L2 leaf regularization 5,
  at least 10 samples per leaf, Newton leaf estimation (10 iterations), log loss, random seed 9448.
- **Sample weights**: negatives 1; positive *i* gets N_neg x p_i / sum(p), where p is the fine-mapping PIP,
  so both classes carry equal total weight and higher-PIP positives weigh more.
- **Feature-category weights**: CatBoost feature weights multiply the score of every candidate split on a
  feature. One weight per category (CRE, Enformer, chromBPNet, ABC, TF, GPN-STAR, gene constraint,
  variant type; all other features 1) was selected for each cell type by a multi-objective search
  (`search_feature_weights.py`) run separately on the odd and on the even chromosomes. A model for a
  held-out odd chromosome uses the weights selected on the even chromosomes and vice versa, so the weights
  never saw the chromosome they are evaluated on.

## Scripts

| Script | Description |
|--------|-------------|
| `evaluate_published_models.py` | Score the test sets with the published models of the data release and report their pooled AUPRC (no training) |
| `train_loco.py` | Train the three models for one cell type and one held-out chromosome; score its test set |
| `evaluate_auprc.py` | Pooled held-out AUPRC per model, with paired bootstrap CIs of the differences |
| `search_feature_weights.py` | Optional: feature-category weight search on the odd or the even chromosomes |
| `select_feature_weights.py` | Optional: pick the weights from the two searches (`best_configs_{cohort}.json`) |
| `run_train.sh` | SLURM array, one task per held-out chromosome |
| `run_search.sh` | SLURM array, odd and even searches |
| `run_pipeline.sh` | All of the above for one cell type on one machine |

The data release already contains the selected weights (`model_training/feature_weights/`), so the
search is only needed to reproduce the selection itself.

## Inputs

From the data release (`paths.release_dir` in `config.yaml`, see the top-level README):

- `model_training/train/{cohort}/`, `train_restricted/{cohort}/`, `test/{cohort}/`
- `model_training/gpn_star/gpn_star_scores_all.parquet`, `gnomad_MAF/`, `gene_lof/`, `columns_dict/`
- `model_training/feature_weights/best_configs_{cohort}.json`
- `model_training/models/{cohort}/`: the published models (`evaluate_published_models.py` only)

Each `annotated_data_{cohort}_chr{N}.parquet` is a directory of parquet part files; keep the downloaded
directories as they are. Files are read in natural order (chromosome, then part number), the order the
published models were trained in, whatever order your filesystem lists them in.

## Evaluating the published models

```bash
cd 5_model_training
python evaluate_published_models.py Mic_mega_eQTL
```

Each held-out chromosome's test set is scored by the published models that held that chromosome out, and
the pooled AUPRC of each model is printed (2-6 minutes per cell type and under 2 GB of memory, no training). For `weighted_full`
(scEEMS) these are the held-out AUPRCs of the manuscript: astrocytes 0.663, excitatory neurons 0.684,
inhibitory neurons 0.639, microglia 0.693, oligodendrocytes 0.762, OPCs 0.758.

To score variants (step 6) or compute SHAP values (step 8) with the published models instead of your own,
set `model_dir: "{published_models_dir}"` in `config.yaml`.

## Quick start: training microglia, one held-out chromosome

About 3.2 GB of data: `train/`, `train_restricted/` and `test/` for `Mic_mega_eQTL`, plus `gpn_star/`,
`gnomad_MAF/`, `gene_lof/`, `columns_dict/` and `feature_weights/`.

```bash
cd 5_model_training
python train_loco.py Mic_mega_eQTL 1
```

This trains the three models with chromosome 1 held out (10-20 minutes with 10 CPU cores, 12 GB of memory)
and writes
`{output_dir}/Mic_mega_eQTL/models/{model}_chr1.joblib` and
`{output_dir}/Mic_mega_eQTL/test_predictions/test_pred_chr1.parquet`. The chromosome 1 test set has 242
variant-gene pairs (22 positives); the expected AUPRCs on it are 0.8015 for `weighted_full`, 0.7555 for
`unweighted_full` and 0.7499 for `weighted_restricted`.

## Full run

```bash
cd 5_model_training
sbatch --export=ALL,cohort=Mic_mega_eQTL run_train.sh      # 22 tasks
python evaluate_auprc.py Mic_mega_eQTL                     # after all 22 have finished
```

or `bash run_pipeline.sh Mic_mega_eQTL` without SLURM. Each held-out chromosome takes 10-30 minutes with
10 CPU cores and needs 12 GB (microglia) to 60 GB (excitatory neurons) of memory.

## Reproducing the feature-weight search

```bash
cd 5_model_training
sbatch --export=ALL,cohort=Mic_mega_eQTL run_search.sh     # odd and even, 83 trials each
python select_feature_weights.py Mic_mega_eQTL
sbatch --export=ALL,cohort=Mic_mega_eQTL,weights=<output_dir>/Mic_mega_eQTL/feature_weight_search/best_configs_Mic_mega_eQTL.json run_train.sh
```

Each trial trains one model per chromosome of its parity that has a test set, so a search is long (days for
the larger cell types); a requeued or resubmitted task resumes. Trials 0-2 are fixed seed configurations,
3-42 Sobol quasi-random and 43-82 Gaussian-process (Optuna GPSampler). The selected "best" weights are the
highest-AUPRC trial of the Gaussian-process phase.

The published search read GPN-STAR scores from an earlier run of the GPN-STAR model on the training and test
variants. The released scores, which the published models use, come from a later genome-wide run and differ
from those in the last digits for a quarter of the SNVs (median absolute difference 0.001). With the released
scores the trial AUPRCs differ from the published search by about 0.001 (trial 0 of the microglia
odd-chromosome search: 0.6441 vs 0.6427), so a rerun can settle on slightly different weights; the released
`feature_weights/` are the weights of the published models.

## Outputs

Under `paths.output_dir`:

- `{cohort}/models/{weighted_full,unweighted_full,weighted_restricted}_chr{N}.joblib`, `feature_cols.pkl`
  (the feature columns and their order, used by steps 6 and 8)
- `{cohort}/test_predictions/test_pred_chr{N}.parquet`: `chrom, variant_id, gene_id, y_true, is_snv,
  pred_weighted_full, pred_unweighted_full, pred_weighted_restricted`
- `{cohort}/test_predictions/auprc.json`
- `{cohort}/feature_weight_search/optuna_{cohort}_{odd,even}.db`, `best_configs_{cohort}.json`
- `{cohort}/published_models/test_pred_chr{N}.parquet` (same columns) and `auprc.json`, from
  `evaluate_published_models.py`

## Notes

- Microglia has no chromosome 21 test set and no chromosome 21 restricted training file (no microglia
  variant on chromosome 21 reaches PIP > 0.9), so `evaluate_auprc.py` pools 21 chromosomes for microglia.
- The excitatory neuron chromosome 19 training file was added to the data release after the excitatory
  neuron models of the manuscript were trained without it. Models trained from the release therefore
  differ slightly from the published ones, except for the model that holds out chromosome 19; the
  held-out AUPRC changes by less than 0.001.
- `is_insertion` and `is_deletion` are swapped relative to their meaning (`shared/featurize.py`); the
  published models were trained this way, so the definition is kept.

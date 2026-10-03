# Step 9: Heritability Annotations (S-LDSC)

Measure how much Alzheimer's disease GWAS heritability the predicted eQTLs explain, with stratified LD score
regression (S-LDSC, PolyFun's implementation) on top of the baseline-LD model, and select the prediction
threshold tau* that step 8 uses to define predicted eQTLs.

## Overview

Four sets of binary annotations are built on the variants of the reference panel. A variant is 1 if it
is in the annotation for any gene (maximum over genes):

| Set | Builder | Columns | Analysis |
|---|---|---|---|
| `pareto/{cohort}_{model}` | `make_annotations.py` | 20: pred_prob > 0.80, 0.81, ..., 0.99 | threshold sweep per model; tau* |
| `pip/{cohort}` | `make_annotations_pip.py` | 1: fine-mapping PIP > 0.10 | fine-mapped comparator of the predicted eQTLs at tau* |
| `top5000/{cohort}` | `select_topn.py`, `make_annotations_topn.py` | 9: top 5,000 variants of each ranking | size-matched comparison |
| `cs/{cohort}` | `make_annotations_cs.py` | 31: 95% credible-set members | credible sets vs predictions; fine-mapping priors |

For each set, `compute_ldscores.py` computes the LD scores of all its columns (per chromosome), and
`ldscore_regression.py` fits the baseline plus ONE column at a time. Columns of a set are nested or
overlapping, so a joint fit would split one signal across them; comparing columns is a ranking over
separate models.

- **tau\***: `tau_star.py` computes the standardized effect size tau\* (Gazal et al. 2017) of every
  threshold with a block-jackknife standard error, and selects, per cell type, the threshold with the
  largest tau\* for `weighted_full`. Enrichment always increases with the threshold; tau\* has an
  interior maximum. In the manuscript tau\* is 0.98 for astrocytes and excitatory neurons, 0.97 for
  inhibitory neurons and microglia, and 0.99 for oligodendrocytes and OPCs; for OPCs no threshold is
  significant.
- **Size matching**: enrichment grows as an annotation gets smaller, so the top-5,000 set gives every
  ranking exactly 5,000 variants: the three models, the FunGen-xQTL fine-mapping PIP (`pip`), and the
  credible sets of the five fine-mapping priors of step 10 (`pip_{prior}`).
- **Credible sets**: `{cell}_cs` is every member of a 95% credible set of the FunGen-xQTL fine-mapping
  (eQTL and other genes, no PIP cut); `{cell}_{prior}_cs_pip{00,10,...,50}` are the credible-set members
  with PIP above 0 (the whole set), 0.10, ..., 0.50 under each prior of step 10, on the eQTL genes
  fine-mapped under all five priors.

Predictions are joined to the reference variants on position and alleles (BP, A1 = alternative allele,
A2 = reference allele). pred_prob is not on the same scale across models (the weighted_restricted model
predicts a much rarer event), so models are never pooled into one ranking.

## What runs from the data release

The scEEMS (`weighted_full`) part of this step runs from the data release: its threshold sweep, tau\* and
the comparison with the fine-mapped eQTLs (PIP > 0.10; the PIPs are in the released predictions). It also
needs the LDSC reference data listed under Inputs, which are not part of the data release. For one cell
type (chain the jobs with `--dependency=afterok:<jobid>`, or wait for each to finish):

```bash
C=Mic_mega_eQTL
python 7_aggregate_predictions/import_release_predictions.py $C      # from the top of the repository
cd 9_create_annotations
sbatch --export=ALL,cohort=$C,model=weighted_full run_annotations.sh
sbatch --export=ALL,set=pareto/${C}_weighted_full run_ldscores.sh
sbatch --export=ALL,set=pareto/${C}_weighted_full --array=1-20 run_ldscore_regression.sh
sbatch --export=ALL,cohort=$C run_annotations_pip.sh
sbatch --export=ALL,set=pip/$C run_ldscores.sh
sbatch --export=ALL,set=pip/$C --array=1 run_ldscore_regression.sh
python tau_star.py
python aggregate_prediction_vs_pip.py
```

The rest cannot be run from the data release: the sweeps of the two comparison models need their
predictions (steps 6-7, which need the step 3 feature tables), the size-matched top-5,000 comparison and the
credible-set annotations need the fine-mapping results of step 10, which come from controlled-access
data, and the `{cell}_cs` column needs the step 1 fine-mapping exports. Their code documents how the
manuscript's results were computed.

## Scripts

| Script | Description |
|--------|-------------|
| `make_annotations.py` | Threshold-sweep annotation, one chromosome x cell type x model |
| `make_annotations_pip.py` | PIP > 0.10 annotation, one chromosome x cell type |
| `select_topn.py` | Genome-wide top-N variants of each ranking, one cell type |
| `make_annotations_topn.py` | Top-N annotation, one chromosome x cell type |
| `make_annotations_cs.py` | Credible-set annotation, one chromosome x cell type |
| `compute_ldscores.py` | LD scores of one annotation set, one chromosome (PolyFun) |
| `ldscore_regression.py` | S-LDSC of one column of one annotation set (PolyFun) |
| `tau_star.py` | tau\* and jackknife p-values of the sweep; selects tau\* per cell type |
| `aggregate_prediction_vs_pip.py` | Predicted eQTLs at tau\* vs fine-mapped eQTLs (PIP > 0.10) |
| `aggregate_topn_sldsc.py` | Size-matched results |
| `aggregate_cs_sldsc.py` | Credible-set results |
| `run_annotations.sh`, `run_annotations_pip.sh`, `run_select_topn.sh`, `run_annotations_topn.sh`, `run_annotations_cs.sh` | SLURM jobs for the builders |
| `run_ldscores.sh`, `run_ldscore_regression.sh` | SLURM arrays for LD scores (per chromosome) and S-LDSC (per column) of any set |
| `run_pipeline.sh` | Everything for one cell type on one machine |

## Inputs

Settings in `config.yaml` (`paths`):

- `ldsc_dir`: a PolyFun installation (https://github.com/omerwe/polyfun), which provides
  `compute_ldscores.py` and `ldsc.py`
- `baseline_annot_dir`: the baseline-LD annotations of the reference-panel variants with their LD scores,
  `baseline_chr{N}.annot.gz`, `baseline_chr{N}.l2.ldscore.parquet`, `baseline_chr{N}.l2.M` (the manuscript
  uses 79 baseline annotations; Gazal et al. 2017)
- `weights_dir`: regression weights, `weights_chr{N}.l2.ldscore.parquet`
- `ldsc_bfile_prefix`: the PLINK reference panel, as a prefix to which the chromosome number is appended
  (`{prefix}{N}.bed/.bim/.fam`)
- `sumstats_file`: GWAS summary statistics munged with PolyFun's `munge_polyfun_sumstats.py`, on the same
  genome build as the panel (the manuscript uses Bellenguez et al. 2022, GRCh38)
- `susie_pips_dir` (default `{data_dir}/susie_vars_pips/{cohort}`): the step 1 fine-mapping exports
  `PIP_top_parquet/` and, if present, `PIP_top_other_parquet/` (credible-set set only)
- `sldsc_dir` (default `{output_dir}/sldsc`): where annotations, LD scores and results are written

From earlier steps:

- step 7: `{predictions_parquet_dir}/{model}/predictions.parquet` (all three models for the sweep and
  top-N; `weighted_full` for the PIP set)
- step 10: `{aggregate_dir}/finemapping_comparison/{cell}_{prior}_cs_variants.tsv` and
  `gene_consensus.tsv` (top-N and credible-set sets only)

## Environments

The builders, `tau_star.py` and the aggregation scripts run in the `scEEMS` environment.
`compute_ldscores.py` and `ldscore_regression.py` run in the PolyFun environment (`POLYFUN_ENV`, default
`polyfun`), with `pyyaml` added (`pip install pyyaml`) because they read `config.yaml`.

## Running the full analysis

For one cell type on SLURM, with the predictions of all three models (step 7) and the step 10 results
(chain the jobs with `--dependency=afterok:<jobid>`, or wait for each to finish):

```bash
cd 9_create_annotations
C=Mic_mega_eQTL

# 1. threshold sweep of the three models, then tau* (tau_star.py summarizes every cell type with results)
for m in weighted_full unweighted_full weighted_restricted; do
    sbatch --export=ALL,cohort=$C,model=$m run_annotations.sh; done
for m in weighted_full unweighted_full weighted_restricted; do
    sbatch --export=ALL,set=pareto/${C}_$m run_ldscores.sh; done
for m in weighted_full unweighted_full weighted_restricted; do
    sbatch --export=ALL,set=pareto/${C}_$m --array=1-20 run_ldscore_regression.sh; done
python tau_star.py

# 2. fine-mapped comparator (PIP > 0.10)
sbatch --export=ALL,cohort=$C run_annotations_pip.sh
sbatch --export=ALL,set=pip/$C run_ldscores.sh
sbatch --export=ALL,set=pip/$C --array=1 run_ldscore_regression.sh
python aggregate_prediction_vs_pip.py

# 3. size-matched top 5,000 (needs step 10)
sbatch --export=ALL,cohort=$C run_select_topn.sh
sbatch --export=ALL,cohort=$C run_annotations_topn.sh
sbatch --export=ALL,set=top5000/$C run_ldscores.sh
sbatch --export=ALL,set=top5000/$C --array=1-9 run_ldscore_regression.sh
python aggregate_topn_sldsc.py

# 4. credible sets (needs step 10)
sbatch --export=ALL,cohort=$C run_annotations_cs.sh
sbatch --export=ALL,set=cs/$C run_ldscores.sh
sbatch --export=ALL,set=cs/$C --array=1-31 run_ldscore_regression.sh
python aggregate_cs_sldsc.py
```

`run_pipeline.sh` runs the same sequence for one cell type on one machine.

LD scores take from about 75 minutes (chromosome 22, 8 cores) to most of a day per chromosome
(chromosome 6, the MHC, is the slowest) and up to about 16 GB of memory. Each S-LDSC run takes 15-35
minutes and up to 60 GB of memory.

## Outputs

Under `{sldsc_dir}`:

- `{set}/MLxQTL_chr{N}.annot.gz`, `.l2.M`, `.l2.ldscore.parquet` for every annotation set
- `results/{set}/{column}.results` (+ `.log`, `.part_delete`, `.delete`): one S-LDSC run per column
- `top{N}_selection/{cohort}_top{N}.tsv.gz`: the selected variants of each ranking

Under `{aggregate_dir}`:

- `tau_star_jackknife.tsv`: every cell type x model x threshold: `Prop._SNPs`, `Prop._h2`, enrichment,
  LDSC coefficient, tau\* with jackknife SE, z and p, h2 and M
- `tau_star.json`: `{cohort: tau*}`, read by step 8
- `prediction_vs_pip_sldsc.tsv`: the predicted (`arm = prediction`, pred_prob > tau\*) and fine-mapped
  (`arm = pip`, PIP > 0.10) eQTL annotations of each cell type, in LDSC's `.results` columns
- `topn_sldsc.tsv`: the size-matched runs (`arm`, `label`, `ranking`, `gene_set`, LDSC columns)
- `finemap_cs_sldsc.tsv`: the credible-set runs (`arm`, `prior`, `pip_threshold`, `gene_set`, LDSC columns)

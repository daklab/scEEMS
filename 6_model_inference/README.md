# Step 6: Model Inference

Score every cis-variant of every gene with the three models of step 5.

## Overview

For each gene, `model_inference.py` featurizes the gene's `all_variants` table from step 3 exactly as in
training (`shared/featurize.py`, with the feature columns and order saved in `feature_cols.pkl`) and
scores it with the three models that held out the gene's chromosome, so no gene is scored by a model that
saw its chromosome in training. The `weighted_full` predictions are the scEEMS predictions; the other two
models' predictions are used as fine-mapping priors in the manuscript's comparisons (step 11).

Genes are those of `list_genes.csv` (MEGA eQTL genes) and `list_genes_other.csv` (other genes) written by
step 3: about 7,100 (microglia) to 11,500 (inhibitory neurons) per cell type.

This step needs the per-gene `all_variants` tables of step 3, which are not part of the data release
(tens of terabytes). The released predictions (`predictions/`) are this step's output for `weighted_full`,
after step 7.

## Scripts

| Script | Description |
|--------|-------------|
| `split_gpn_by_chr.py` | One-time: split the GPN-STAR scores into one file per chromosome |
| `make_gene_lists.py` | List the genes to score per cell type; prints the SLURM array size |
| `model_inference.py` | Score one gene (one row of the gene list) with the three models |
| `find_missing_genes.py` | List genes not scored yet, for resubmission |
| `run_inference.sh` | SLURM array, one task per gene |
| `run_pipeline.sh` | All genes of one cell type on one machine |

## Running

```bash
cd 6_model_inference
python split_gpn_by_chr.py
python make_gene_lists.py Mic_mega_eQTL                   # prints --array=1-7100%200
sbatch --export=ALL,cohort=Mic_mega_eQTL --array=1-7100%200 run_inference.sh
python find_missing_genes.py Mic_mega_eQTL                # resubmit any rows it lists
```

Each task takes about a minute; the largest genes need up to 30 GB of memory.

## Inputs

- `{model_dir}`: the 66 models of step 5 (`{model}_chr{N}.joblib`) and `feature_cols.pkl`
- `{all_variants_dir}/{gene_id}/`: per-gene feature tables from step 3
- `{gene_list_dir}/list_genes.csv`, `list_genes_other.csv`: from step 3
- the GPN-STAR, gnomAD MAF, gene constraint and columns dictionary files of the data release

## Outputs

- `{output_dir}/{cohort}/list_genes_all.csv`
- `{output_dir}/{cohort}/predictions/{weighted_full,unweighted_full,weighted_restricted}/{gene_id}_predictions.tsv`
  with columns `variant_id, chr, pos, ref, alt, pip, gene_id, pred_prob`. `pip` is the FunGen-xQTL
  fine-mapping PIP of the variant for the gene (step 1); `pred_prob` is the model's probability that the
  variant is a causal eQTL for the gene.

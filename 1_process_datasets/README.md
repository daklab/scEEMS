# Step 1: Process Datasets

Turn the FunGen-xQTL SuSiE eQTL fine-mapping results into the per-variant tables the rest of the pipeline
reads.

**This step needs the FunGen-xQTL fine-mapping exports, which are not distributed with this repository or
the data release.** Its scripts document how the released training data were made; steps 5-11 start from
the data release (see the top-level README).

## Overview

The exports hold, for each gene, the SuSiE fine-mapping of every eQTL analysis that tested it: the MEGA
analysis of each cell type (e.g. `Mic_mega_eQTL`), the individual studies (e.g. `Mic_DeJager_eQTL`,
`Mic_Kellis_eQTL`) and other eQTL analyses. For each analysis, this step writes:

- `PIP_all`: the PIP of every variant of each gene's cis window;
- `PIP_top`: the top loci of each gene, with the variants' credible-set membership at 95%, 70% and 50%
  coverage (`cs_coverage_0.95` etc.: 0 outside any credible set, otherwise the credible set's index).

The scEEMS models are trained on the MEGA analyses. A cell type's **other genes** are the genes its MEGA
analysis did not fine-map but its DeJager or Kellis analysis did. `create_parquet_files_other.py` collects
their results into the MEGA analysis' folder, so that step 3 featurizes them, step 6 scores them and step 9
includes their credible sets.

Finally, `create_unique_variant_list.py` lists the unique variants of each chromosome over all analyses:
the variants step 2 annotates.

## Scripts

| Script | Description |
|--------|-------------|
| `create_job_list.py` | List the genes of the fine-mapping exports |
| `get_vars_pips.R` | Extract one gene's PIPs and top loci for every analysis (R) |
| `create_parquet_files.py` | Per-gene CSV tables of every analysis -> parquet datasets partitioned by chromosome |
| `create_parquet_files_other.py` | Collect a cell type's other genes into its MEGA analysis' folder |
| `create_unique_variant_list.py` | Unique variants of one chromosome over all analyses |
| `download_data.py` | Explains that the fine-mapping exports are not distributed |
| `run_pipeline.sh` | All of the above, in order |

## Running

```bash
bash 1_process_datasets/run_pipeline.sh
```

The Python scripts run in the `scEEMS` environment (`environment.yml`) and `get_vars_pips.R` in the
`scEEMS_R` environment (`environment_r.yml`); `run_pipeline.sh` switches between them. With about 18,000
genes, run `get_vars_pips.R` as an array job (one task per line of `genes_list.txt`) on a cluster.

## Inputs

- `{finemapping_rds_dir}/Fungen_xQTL.<gene_id>.cis_results_db.export.rds` and
  `Fungen_xQTL.<gene_id>.cis_results_db.export_sumstats.rds`, one pair per gene
  (default `finemapping_rds_dir`: `{data_dir}/release_04_2024`)

## Outputs

Under `susie_dir` (default `{data_dir}/susie_vars_pips`), one folder per eQTL analysis
(`susie_pips_dir`, `{susie_dir}/<analysis>`):

- `genes_list.txt`: the genes of the exports; a gene's line number is its index for `get_vars_pips.R`
- `<analysis>/PIP_all/<gene_id>_pips.csv` and `<analysis>/PIP_top/<gene_id>_top_loci.csv`
- `<analysis>/PIP_all_parquet/PIP_all.parquet/` and `<analysis>/PIP_top_parquet/PIP_top.parquet/`,
  partitioned by `chr`; `PIP_all` has columns `variant_id, pos, ref, alt, pip, gene_id`
- `<cell>_mega_eQTL/PIP_all_other_parquet/PIP_all_other.parquet/` and
  `<cell>_mega_eQTL/PIP_top_other_parquet/PIP_top_other.parquet/`: the other genes of each cell type
- `variant_list/variant_list_chr{N}.parquet` (`variant_list_dir`): `variant_id, chr, pos, ref, alt`

## Notes

- `PIP_all_other` holds every other gene (a gene fine-mapped by both DeJager and Kellis contributes both
  studies' rows), while `PIP_top_other` holds only the other genes fine-mapped by exactly one of the two
  studies. This is how the released data were built.
- The DLPFC_Klein analyses in the exports are skipped.

# Step 7: Aggregate Predictions

Collect the per-gene predictions of step 6 into one parquet dataset per cell type and model, and export
the scEEMS predictions in the tabix-indexed format of the data release.

## Scripts

| Script | Description |
|--------|-------------|
| `create_parquet_scored.py` | Per-gene TSVs of one model -> parquet dataset partitioned by chromosome |
| `extract_predictions_tsv.py` | scEEMS (`weighted_full`) predictions of one chromosome -> bgzipped, tabix-indexed TSV |
| `import_release_predictions.py` | The reverse: released TSVs -> the `weighted_full` parquet dataset, to run steps 9-10 from the data release |
| `run_aggregate.sh` | SLURM: the three models of one cell type |
| `run_export.sh` | SLURM array: the 22 chromosomes of one cell type |
| `run_pipeline.sh` | Everything, all cell types, on one machine |

## Running

```bash
cd 7_aggregate_predictions
sbatch --export=ALL,cohort=Mic_mega_eQTL run_aggregate.sh     # after step 6 has scored every gene
sbatch --export=ALL,cohort=Mic_mega_eQTL run_export.sh        # after run_aggregate.sh
```

## Outputs

- `{output_dir}/{cohort}/predictions_parquet/{model}/predictions.parquet/chr=chr{N}/`: columns
  `variant_id, pos, ref, alt, pip, gene_id, pred_prob`, partitioned by `chr`. Steps 8-12 read these.
- `{output_dir}/release/predictions/{cohort}/predictions_{cohort}_{N}.tsv.gz` (+ `.tbi`): columns
  `#CHROM POS ID REF ALT GENE_ID PIP PRED_PROBABILITY`, `CHROM` without the `chr` prefix, sorted by
  position. These are the files of the `predictions/` folder of the data release.

## Starting from the data release

The scEEMS (`weighted_full`) analyses of steps 9 and 10 read the `weighted_full` parquet dataset, which can
be rebuilt from the released predictions without steps 5-7 (the comparison models' predictions are not
released):

```bash
python import_release_predictions.py Mic_mega_eQTL     # reads {release_dir}/predictions/Mic_mega_eQTL/
```

# Step 3: Cell Featurization

Build the feature table of every gene: each variant of the gene's cis window with its fine-mapping PIP and
all model features except those step 5 adds (gene constraint, gnomAD MAF, GPN-STAR, variant type).

**This step needs the outputs of steps 1 and 2 and several external inputs.** Its outputs (`all_variants`,
tens of terabytes over all cell types) are not distributed: step 4 samples the training and test data of
the data release from them, and steps 6 and 8 read them to score and explain every variant.

## Overview

`create_gene_lists.py` lists a cell type's genes: `list_genes.csv`, the genes its MEGA analysis
fine-mapped, and `list_genes_other.csv`, its other genes, fine-mapped only by its DeJager or Kellis
analysis (step 1). `create_gene_datasets.py` then builds one gene's table from the gene's fine-mapped
variants that have annotations (step 2), adding:

- **distance to the TSS**: `distance_TSS` = TSS - position (TSS from `ABC_gene_id_name_mapping.csv`, the
  mean over the gene's entries), `abs_distance_TSS` and `abs_distance_TSS_log`;
- **ABC scores** in microglia, astrocytes, oligodendrocytes and neurons: the highest ABC score of the
  gene's predicted enhancer-gene links within 1,024 bp of the variant (0 if none);
- **ChromBPNet scores** in the four cell types: for each assay and metric (log counts difference, summed
  absolute log probability difference, Jensen-Shannon distance), over the variant's peaks, the largest
  value (at least 0), the smallest value (at most 0) and the largest absolute value;
- **composite TF scores** per cell type: maximum, minimum and maximum absolute Enformer ChIP-seq prediction
  over the cell type's transcription factors;
- all variant annotations of step 2.

A variant that appears more than once for a gene keeps the row with the highest PIP.

## Scripts

| Script | Description |
|--------|-------------|
| `create_gene_lists.py` | List a cell type's MEGA genes and other genes |
| `create_gene_datasets.py` | Build the feature table of one gene (`F`: MEGA gene, `T`: other gene) |
| `run_pipeline.sh` | Both, for every gene of one cell type |

## Running

```bash
bash 3_cell_featurization/run_pipeline.sh Mic_mega_eQTL
```

There are 7,000 to 11,500 genes per cell type. A gene takes about 10-15 minutes with 4 CPU cores and up to
30 GB of memory (the original runs allowed 2 hours per gene); on a cluster, run
`create_gene_datasets.py <row> <cohort> F` (and `T`) as an array job over the rows of the two gene lists.

## Inputs

| Setting | Content |
|---|---|
| `susie_pips_dir` | `PIP_all_parquet/` and `PIP_all_other_parquet/` of the cell type's MEGA analysis (step 1) |
| `variant_list_dir` | `annotated_variants/` and `annotated_variants_just_variants/` (step 2) |
| `abc_data_dir` | `ABC_gene_id_name_mapping.csv` (`gene_id, gene_name, gene_TSS`) and `ABC_results_<cell>_v2/<cell>/Predictions/EnhancerPredictionsAllPutative.tsv.gz` for microglia, astrocyte, oligodendrocyte and neuron |
| `chrombpnet_dir` | `chrombpnet_<cell>_chr{N}_variant_peak_pairs_scored.csv`: ChromBPNet scores of variant-peak pairs per cell type and chromosome |
| `tf_file`, `targets_file` | as in step 2 |

## Outputs

- `{gene_list_dir}/list_genes.csv` and `list_genes_other.csv` (default `gene_list_dir`:
  `{data_dir}/training_data/{cohort}`): `gene_id, chr`; a gene's row number is its index for
  `create_gene_datasets.py`
- `{all_variants_dir}/<gene_id>/annotated_data_<cohort>_<chr>.parquet` (default `all_variants_dir`:
  `{data_dir}/training_data/{cohort}/all_variants`): one row per variant with `pip` and the features

## Notes

- `gene_TSS` in `ABC_gene_id_name_mapping.csv` is the gene's lowest genomic coordinate for genes on either
  strand, so for genes on the minus strand `distance_TSS` is measured from the gene's 3' end. The released
  features and models use this definition.

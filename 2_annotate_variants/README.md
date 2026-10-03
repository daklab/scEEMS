# Step 2: Annotate Variants

Annotate the unique variants of step 1 with the variant-level features of the model.

**This step needs large external inputs** (Enformer predictions for every variant, brain and baseline BED
annotations) that are not distributed with this repository. The features it produces are part of the
training and test data of the data release.

## Overview

For each chromosome, `annotate_variants.py` joins the unique variants of step 1 with:

1. **Enformer variant effect predictions** (Avsec et al. 2021, https://doi.org/10.1038/s41592-021-01252-x):
   the `diff_32_<track>` score of every Enformer output track except the CAGE tracks. Only variants with
   Enformer predictions are kept.
2. **Brain cell type annotations** (`bed_annotations_dir`): ATAC-seq peaks, promoters and enhancers of sorted
   brain nuclei (PU.1+ microglia, NeuN+ neurons, LHX2+ astrocytes, Olig2+ oligodendrocytes), and their
   intersections and unions; 1 if the variant falls in the region, else 0.
3. **Baseline genomic annotations** (`baseline_annotations_dir`): one 0/1 column per BED file.
4. **Composite TF scores**: for each cell type, the maximum Enformer ChIP-seq prediction over the cell
   type's transcription factors (`none-none-<cell>`; `none-none-all` over every TF of `tf_file` and the
   cell types' lists), and its product with the cell type's promoter, enhancer or ATAC annotation
   (`<cell>_enhancer_promoter_union_atac_500`).

`subset_annotated_variants.py` then keeps only the identifier columns, which steps 3 and 4 join to the
fine-mapped variants to keep those with annotations.

The other features are added later: ABC scores, ChromBPNet scores and distance to the TSS per gene in
step 3, and gene constraint, gnomAD MAF and GPN-STAR in step 5 (from the data release).

## Scripts

| Script | Description |
|--------|-------------|
| `annotate_variants.py` | Annotate the variants of one chromosome |
| `subset_annotated_variants.py` | Keep the identifier columns of one chromosome's annotated variants |
| `run_pipeline.sh` | Both, for all chromosomes |

## Inputs

| Setting | Content |
|---|---|
| `variant_list_dir` | `variant_list_chr{N}.parquet` from step 1 |
| `enformer_dir` | `enformer_tensorflow_chr{N}.parquet`: Enformer predictions per variant (`CHR, BP, REF, ALT, SNP` and one `diff_32_<track identifier>` column per track) |
| `bed_annotations_dir` | the brain BED files named in `annotate_variants.py`, e.g. `PU1_optimal_peak_IDR_ENCODE.ATAC.bed`, `PU1_promoter.bed`, `PU1_enhancer.bed` |
| `baseline_annotations_dir` | baseline annotation BED files; every `.bed` file except `*.unmapped.bed` becomes a column |
| `tf_file` | table of TF ChIP-seq experiments with a `Target of assay` column (an ENCODE experiment report; its first line is skipped) |
| `targets_file` | Enformer's `targets_human.txt`, the description of each output track |

## Running

```bash
bash 2_annotate_variants/run_pipeline.sh
```

Run the chromosomes as separate jobs on a cluster. Chromosome 21, one of the smallest, takes about 45
minutes and 62 GB of memory with 4 CPU cores; the largest chromosomes need several times more (the original
runs allowed 10 CPU cores, 150 GB and 14 hours per chromosome).

## Outputs

- `{variant_list_dir}/annotated_variants/annotated_variants_chr{N}.parquet`: one row per variant (index
  `variant_id`), with the Enformer scores, brain and baseline annotations and composite TF scores
- `{variant_list_dir}/annotated_variants_just_variants/annotated_just_variants_chr{N}.parquet`: index
  `variant_id` with `CHR, BP, REF, ALT, SNP`

## Notes

- `pybedtools` writes temporary files to the system temporary directory; set `TMPDIR` to change it.

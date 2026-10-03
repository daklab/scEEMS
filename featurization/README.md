# Featurizing and scoring new variants

Compute the 4,840 features of the scEEMS models for any set of variant-gene pairs and score them with the
published models: for each of the six brain cell types, the probability that the variant is a causal eQTL
variant of the gene. The features are computed as for the training data (steps 2 and 3 of the pipeline),
so the variants need not be in the data release: rare and novel variants, SNVs and indels, can be scored.

## What you need

**1. The data release**: `model_training/` (the published models, gene constraint, columns dictionary and
the GPN-STAR scores of the training variants) and `featurization/` (1.75 GB):

```bash
python download_synapse_data.py --resource featurization model_training
```

| `featurization/` file | Content |
|---|---|
| `chrombpnet_models/` | ChromBPNet (ATAC) and BPNet (H3K27ac, H3K4me3) models, 5 folds per cell type and assay |
| `chrombpnet_peaks.tsv.gz` | the peaks of each cell type and assay that variants are paired with |
| `celltype_annotations.bed.gz` | the 40 cell type ATAC/promoter/enhancer annotations |
| `abc_scores.tsv.gz` | ABC enhancer-gene links of four brain cell types |
| `genes.tsv` | gene ID, name, chromosome and TSS of every gene a variant can be paired with |
| `targets_human.txt` | the Enformer output tracks |

**2. Public resources**, downloaded with their own tools (no account needed):

```bash
# GRCh38 (the ENCODE analysis set, the genome the features were computed on)
curl -L -o GRCh38.fa.gz https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
gunzip GRCh38.fa.gz

# GPN-STAR: the genome-wide scores of the primate model, one file per chromosome you need (0.8-4.8 GB each),
# and the model's calibration table
R=https://huggingface.co/datasets/songlab/gpn-star-scores/resolve/a7b13bbf0d2338d74a7e5f0f8466e41ac0722f50
for c in chr11; do curl -L --create-dirs -o gpn-star/llr_$c.parquet $R/data/gpn-star-hg38-p243-200m/llr/llr_$c.parquet; done
curl -L -o gpn-star/calibration_llr.parquet https://huggingface.co/songlab/gpn-star-hg38-p243-200m/resolve/16773016130e826cd3d72e91ed312a3cb27ba2b5/calibration_table/llr.parquet

# gnomAD v3.1 genomes, one VCF and index per chromosome you need
G=https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes
for c in chr11; do curl -L --create-dirs -o gnomad/gnomad.genomes.v3.1.sites.$c.vcf.bgz $G/gnomad.genomes.v3.1.sites.$c.vcf.bgz
                   curl -L -o gnomad/gnomad.genomes.v3.1.sites.$c.vcf.bgz.tbi $G/gnomad.genomes.v3.1.sites.$c.vcf.bgz.tbi; done
```

The Enformer model (TensorFlow Hub `deepmind/enformer/1`, about 1 GB) is downloaded on first use; set
`TFHUB_CACHE_DIR` to keep it. For a few variants the gnomAD VCFs need not be downloaded: `--vcf` also accepts
their URLs (`$G/gnomad.genomes.v3.1.sites.{chrom}.vcf.bgz`), and only the variants' positions are read.

**3. Three conda environments**: `environment.yml` (scEEMS, steps 1 and 4-7), `environment_enformer.yml`
(step 2) and `environment_chrombpnet.yml` (step 3):

```bash
conda env create -f environment.yml
conda env create -f environment_enformer.yml
conda env create -f environment_chrombpnet.yml
```

Steps 2 and 3 run on a GPU if there is one (TensorFlow 2.15 and PyTorch 2.5 with CUDA 12; tested on NVIDIA
L40S), otherwise on the CPU. If packages installed with `pip install --user` are on your system, set
`PYTHONNOUSERSITE=1` when creating and using these environments.

## Input

A tab-separated file with a header row and one variant per row, in GRCh38 with 1-based positions:

| column | |
|---|---|
| `chrom`, `pos`, `ref`, `alt` | the variant (or one `variant_id` column, `chrom:pos:ref:alt`); SNVs, insertions, deletions and multi-base substitutions; autosomes |
| `name` | optional, e.g. an rsID; carried to the outputs |
| `gene_id` | optional Ensembl gene ID: the gene the variant is scored for (one row per variant-gene pair); a variant without one is paired with every gene whose TSS is within 1 Mb |

The example scores the two PICALM/EED locus variants for both genes ([example/picalm_eed.tsv](example/picalm_eed.tsv)):

```
chrom   pos       ref  alt  name        gene_id
chr11   86156833  A    G    rs10792832  ENSG00000073921
chr11   86156833  A    G    rs10792832  ENSG00000074266
chr11   86157598  T    C    rs3851179   ENSG00000073921
chr11   86157598  T    C    rs3851179   ENSG00000074266
```

## Running

```bash
OUT=out/picalm_eed
conda activate scEEMS
python featurization/1_prepare_variants.py featurization/example/picalm_eed.tsv $OUT --fasta GRCh38.fa
conda activate scEEMS_enformer
python featurization/2_score_enformer.py $OUT --fasta GRCh38.fa
conda activate scEEMS_chrombpnet
python featurization/3_score_chrombpnet.py $OUT --fasta GRCh38.fa
conda activate scEEMS
python featurization/4_score_gpn_star.py $OUT --scores 'gpn-star/llr_{chrom}.parquet' \
    --calibration gpn-star/calibration_llr.parquet --fasta GRCh38.fa
python featurization/5_extract_gnomad_maf.py $OUT --vcf 'gnomad/gnomad.genomes.v3.1.sites.{chrom}.vcf.bgz'
python featurization/6_build_features.py $OUT
python featurization/7_score_variants.py $OUT
```

Steps 2-5 depend only on step 1 and can run in parallel.

| Script | Environment | Output in OUT | |
|---|---|---|---|
| `1_prepare_variants.py` | scEEMS | `variants.tsv`, `pairs.tsv` | checks the variants against GRCh38 (REF/ALT swapped to the genome's orientation if needed) and pairs them with genes |
| `2_score_enformer.py` | scEEMS_enformer | `enformer.parquet` | Enformer `diff_32` of 5,313 tracks per variant (about 1 s per variant on a GPU) |
| `3_score_chrombpnet.py` | scEEMS_chrombpnet | `chrombpnet.tsv` | ChromBPNet scores of every variant-peak pair (peak centre within 1,024 bp) |
| `4_score_gpn_star.py` | scEEMS | `gpn_star.tsv` | GPN-STAR log-likelihood ratio per SNV |
| `5_extract_gnomad_maf.py` | scEEMS | `gnomad_MAF.tsv` | gnomAD allele frequency per variant |
| `6_build_features.py` | scEEMS | `features.parquet` | the features of every variant-gene pair |
| `7_score_variants.py` | scEEMS | `predictions.tsv`, `feature_matrix.parquet` | the scEEMS (Weighted (Full)) prediction of each pair in each cell type, and the model input |

`7_score_variants.py` scores each pair with, for each cell type, the published model that was trained
without the variant's chromosome (as the released predictions were); `--models` adds the Unweighted (Full)
and Weighted (Restricted) models. `feature_matrix.parquet` holds the 4,840 features exactly as the models
receive them, the same for every cell type.

## How the features are computed

| Features | Computed from | As in training |
|---|---|---|
| Enformer (4,675) | reference minus alternate prediction, summed over the central 32 output bins (4,096 bp), of the 393,216 bp sequence centred on the variant; CAGE tracks dropped | step 2 |
| ChromBPNet (36) | per cell type, score and assay, the largest absolute score over the variant's peaks (0 near no peak); scores averaged over the 5 fold models | step 3 |
| TF (4) | per cell type, the largest absolute Enformer ChIP `diff_32` of its transcription factors | step 3 |
| cell type CRE (44) | whether the variant position lies in each annotation; per cell type, its promoter/enhancer/ATAC union times its TF score | step 2 |
| ABC (4) | per cell type, the largest ABC score of the gene's enhancer links within 1,024 bp of the variant | step 3 |
| distance (1) | log of the distance between the variant and the gene's TSS (`genes.tsv`) | step 3 |
| GPN-STAR (1) | absolute log-likelihood ratio of the primate GPN-STAR model; for SNVs only | step 5 |
| variant type (5), gene constraint (1), gnomAD (1) | REF/ALT lengths, GeneBayes constraint of the gene, gnomAD v3.1 genomes INFO/AF | step 5 |
| baseline annotations (68) | 0, their value throughout the training data | |

The sequences of steps 2 and 3 are built as the official Enformer usage notebook builds them (kipoiseq),
reimplemented in `variant_sequence.py` so the environments do not need kipoiseq. GPN-STAR scores of new
SNVs come from GPN-STAR's published genome-wide scores of the model that scored the training variants
(`gpn-star-hg38-p243-200m`), with their mutation-rate calibration removed; training variants keep their
training score. Variants missing from gnomAD (and genes without a constraint estimate) get the median of
the scored pairs, as in training.

## Agreement with the training features

Tested against the training data of the data release:

- **Feature building** (step 6), from the inputs of the training data: all 4,839 features identical for
  the 1,892 held-out microglia test pairs (the TSS distance after conversion to the 32-bit floats the
  models use).
- **Sequences** (`variant_sequence.py`): identical to kipoiseq 0.7.1, including the one-hot encodings, in
  all 65,076 tested windows (5,707 Enformer, 59,369 ChromBPNet): random and real training variants, SNVs,
  insertions and deletions up to 900 bp and multi-base substitutions, REF-mismatched variants and windows
  at chromosome ends.
- **ChromBPNet** (step 3), 59 test variants: the same 99 variant-peak pairs as in training, and the same
  scores up to rounding (log count difference within 3e-5).
- **Enformer** (step 2): the same scores as in training for variants scored by the same procedure, up to
  GPU rounding, which also varies between runs (TensorFlow's default TF32 arithmetic): a difference of up
  to about 0.1% of the 32-bin sums, correlation over the tracks >= 0.998 per variant. Some training
  variants were scored on an earlier version of the genome or with a variant of the procedure; their
  scores differ more, most for indels.
- **GPN-STAR** (step 4): within 0.01 of the training scores of 5,000 chromosome 22 SNVs.
- **gnomAD** (step 5): identical to the training values for all 1,648 test variants in gnomAD (v3.1
  genomes; gnomAD v3.1.2 and v4 differ for some variants).
- **End to end**, the example above: the 24 predictions (4 variant-gene pairs x 6 cell types) differ from
  the released predictions of the same pairs by at most 0.018 (microglia, PICALM: rs10792832 0.977 vs 0.977,
  rs3851179 0.831 vs 0.828).

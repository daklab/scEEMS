# Step 10: eMAGMA Gene Analysis

Test which genes carry Alzheimer's disease GWAS signal through the variants scEEMS links to them, with
MAGMA gene analysis on eQTL-informed gene annotations (eMAGMA), and replicate the results in non-European
GWAS.

Method references:
- MAGMA: de Leeuw CA, Mooij JM, Heskes T, Posthuma D. MAGMA: Generalized Gene-Set Analysis of GWAS Data. *PLoS Computational Biology*. 2015;11(4):e1004219. https://doi.org/10.1371/journal.pcbi.1004219
- E-MAGMA: Gerring ZF, Mina-Vargas A, Gamazon ER, Derks EM. E-MAGMA: an eQTL-informed method to identify risk genes using genome-wide association study summary statistics. *Bioinformatics*. 2021;37(16):2245-2249. https://doi.org/10.1093/bioinformatics/btab115

## Overview

For each cell type, a gene's SNP set is defined in two ways (`make_magma_files.py`):

- **prediction**: variants the scEEMS model calls eQTLs of the gene, pred_prob > tau* (the cell type's
  threshold from the S-LDSC analysis of step 9). Each variant is assigned to a single gene: the gene with
  the most PLAC-seq enhancer-promoter interactions overlapping the variant in the matched brain cell type
  (Nott et al. 2019: PU.1 microglia, NeuN neurons, Olig2 oligodendrocytes), ties broken by the highest ABC
  score (>= 0.005) in that cell type; variants with neither are dropped. Astrocytes have no PLAC-seq data,
  so ABC alone.
- **pip**: variants with fine-mapping PIP > 0.10 for the gene.

Two comparison annotations:

- **TSS-distance control** (`make_magma_files_knn.py`): for every gene, the same number of variants, taken
  as the variants closest to its TSS, so significance through proximity and SNP count alone can be
  separated from the model's choice of variants.
- **positional MAGMA** (`make_positional_annotation.py`): SNPs within gene bodies, the standard MAGMA
  annotation.

The gene analyses (`run_magma.py`, MAGMA v1.10, SNP-wise mean model) use the European AD GWAS of Bellenguez
et al. 2022 with the ADSP European reference panel for LD. The eMAGMA annotations are also tested against
the ADGC African American (AFR), Hispanic/Latino (AMR) and East Asian (EAS) AD GWAS with the ADSP reference
panel of the same population. The prediction annotation therefore also contains variants of a
non-European variant panel (SNVs with rsIDs from the ADSP AFR, AMR and EAS reference panels), scored by the
scEEMS model (`score_noneur.py`).

Genes are significant at P < 0.05 / 18,700 (Bonferroni); replication in a non-European GWAS is P < 0.05.

## Scripts

| Script | Description |
|--------|-------------|
| `make_magma_sumstats.py` | GWAS summary statistics in MAGMA's format (European, ADGC AFR/AMR/EAS) |
| `make_positional_annotation.py` | Positional (gene-body) MAGMA annotation |
| `score_noneur.py` | Score the non-European panel with the scEEMS model (one cell type, one chromosome) |
| `make_magma_files.py` | eMAGMA prediction and pip annotations (one cell type, one chromosome) |
| `make_magma_files_knn.py` | TSS-distance control annotations (one cell type, one chromosome) |
| `run_magma.py` | MAGMA gene analysis of one annotation for one chromosome |
| `aggregate_magma.py` | Significant genes, multi-ancestry table, positional baseline |
| `run_score_noneur.sh`, `run_make_magma.sh`, `run_magma.sh` | SLURM arrays, one task per chromosome |
| `run_pipeline.sh` | Everything, in order, on one machine |
| `noneur_gpn_star/` | Optional: GPN-STAR scores of the non-European panel's SNVs (GPU) |

## Inputs

From earlier steps:
- `{output_dir}/{cohort}/predictions_parquet/weighted_full/predictions.parquet` (step 7, or
  `7_aggregate_predictions/import_release_predictions.py` from the data release)
- `{aggregate_dir}/tau_star.json` (step 9)
- the `weighted_full` models of step 5 (non-European scoring only)

External, set in `config.yaml` (`paths`):

| Setting | Content |
|---|---|
| `magma_binary`, `magma_gene_loc_file` | MAGMA v1.10 executable and its NCBI build 38 gene locations (`NCBI38.gene.loc`) |
| `plink_ref_dir` | ADSP reference panels: `{plink_ref_dir}/{EUR,AFR,AMR,EAS}/plink/ADSP_{POP}_chr{N}.{bed,bim,fam}` |
| `baseline_annot_dir` | S-LDSC baseline annotations `baseline_chr{N}.annot.gz` (European variant set with SNP IDs; as in step 9) |
| `abc_data_dir` | ABC predictions `ABC_results_{cell}_v2/{cell}/Predictions/EnhancerPredictionsAllPutative.tsv.gz` and `ABC_gene_id_name_mapping.csv` (as in steps 3-4) |
| `placseq_file` | PLAC-seq interactions: `gene_id, start, end, interaction_type` (`PU1_`, `Olig2_`, `NeuN_enhancer_interactions`) |
| `sumstats_file` | Bellenguez et al. 2022 summary statistics, munged parquet with `SNP CHR BP P N` (as in step 9) |
| `adgc_sumstats_dir` | ADGC summary statistics `{AFA,HISP,EAS}_common_apoe_adj_p-valueOnly.txt` |
| `gene_info_file` | NCBI `Homo_sapiens.gene_info.gz` (gene symbols, Entrez to Ensembl) |
| `multi_ancestry_predictions_dir` | Non-European panel, per cell type: `MAGMA_features_{cohort}_chr{N}.parquet` (features) and `MAGMA_predictions_{cohort}_chr{N}.tsv.gz` (variant and gene columns) |
| `noneur_gpn_star_file` | GPN-STAR scores of the panel's SNVs (`noneur_gpn_star/`) |
| `noneur_variants_file` | The panel's variants (`CHR BP REF ALT`), input of `noneur_gpn_star/` |
| `magma_dir` | Where this step writes (for example `{output_dir}/magma`) |

## Running

```bash
cd 10_magma_analysis
python make_magma_sumstats.py                     # about 15 GB of memory, 10 minutes
python make_positional_annotation.py
sbatch --export=ALL,annotation=positional run_magma.sh
for c in Ast Exc Inh Mic Oli OPC; do sbatch --export=ALL,cohort=${c}_mega_eQTL run_score_noneur.sh; done
# when those have finished:
for c in Ast Exc Inh Mic Oli OPC; do sbatch --export=ALL,cohort=${c}_mega_eQTL run_make_magma.sh; done
# when those have finished:
for c in Ast Exc Inh Mic Oli OPC; do
    for pop in EUR AFR AMR EAS; do sbatch --export=ALL,annotation=emagma,cohort=${c}_mega_eQTL,pop=${pop} run_magma.sh; done
    sbatch --export=ALL,annotation=knn,cohort=${c}_mega_eQTL run_magma.sh
done
# when all MAGMA jobs have finished:
python aggregate_magma.py
```

The non-European GPN-STAR scores (`noneur_gpn_star/`, needed once before `score_noneur.py`):

```bash
cd noneur_gpn_star
python aggregate_noneur_variants.py               # SNVs per chromosome
sbatch --export=ALL,dir=<magma_dir>/noneur_gpn_star,msa=<GPN-STAR MSA>,model=<GPN-STAR model> run_score_gpn_star.sh
python merge_scores_noneur.py                     # -> noneur_gpn_star_file
```

This needs a GPU and the `gpn` package with the GPN-STAR model and its alignment
(https://github.com/songlab-cal/gpn).

## Outputs

Under `{magma_dir}`:
- `sumstats/{bellenguez,ADGC_AFR,ADGC_AMR,ADGC_EAS}_MAGMA_sumstats.txt`
- `noneur_predictions/{cohort}/MAGMA_predictions_{cohort}_chr{N}.tsv.gz`
- `{cohort}/{prediction,pip}/{cohort}_chr{N}_MAGMA.genes.annot`, `{cohort}/variant_gene_lists/`
- `{cohort}/MAGMA_output{,_ADGC_AFR,_ADGC_AMR,_ADGC_EAS}/{cohort}_chr{N}_MAGMA_{prediction,pip}.genes.out`
- `knn/{cohort}/{prediction,pip,MAGMA_output}/`: the TSS-distance control, same layout
- `positional/bellenguez_MAGMA.genes.annot`, `positional/MAGMA_output/positional_chr{N}_MAGMA.genes.out`

Under `{aggregate_dir}/magma/` (`aggregate_magma.py`):
- `MAGMA_significant_genes_{pred,pip}.txt`: significant (gene, cell type) pairs, European GWAS; with
  `_no_overlap`, without genes significant in positional MAGMA
- `MAGMA_knn_significant_genes_{pred,pip}[_no_overlap].txt`: the same for the TSS-distance control
- `MAGMA_multiancestry_genes_combined{,_pip}.txt`: every eMAGMA gene result with `P_EUR`, `P_AFR`,
  `P_AMR`, `P_EAS` (1 where a gene was not tested), `gene_TSS` and `GENE_NAME`
- `MAGMA_positional_genes_combined.txt`: positional MAGMA results with Ensembl IDs

## Notes

- `run_magma.py` asks MAGMA for adaptive permutation p-values (`--gene-settings adap-permp=10000`); the
  `PERMP` column changes between runs. The `P` column used above is reproducible on a given machine; on
  machines with different processors it can differ in the last printed digit.
- The non-European panel's feature tables contain the same features as the main variant set (step 3) but
  were built by code that is not part of this repository; `score_noneur.py` only rescores them with the
  step 5 models.

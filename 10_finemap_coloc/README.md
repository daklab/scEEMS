# Step 10: Fine-Mapping with scEEMS Priors and Colocalization with AD GWAS

**This step needs controlled-access data.** It fine-maps eQTLs from individual-level ROSMAP genotypes
and single-nucleus pseudobulk expression, available through the AD Knowledge Portal (Synapse) under a
data use agreement, and computes GWAS LD from the ADSP whole-genome sequencing reference panel, available
through NIAGADS. The code is provided so that the analysis can be reproduced by anyone with access to
those data.

**Its results are in the data release:** `fine_mapping/` holds the eQTL credible sets of all six cell
types under the five priors (see its README). Most users should start from those files.

## Overview

For each gene and cell type, SuSiE fine-maps the gene's cis-eQTLs under five priors:

| Prior | Prior weights |
|---|---|
| `uniform` | none (every variant equal) |
| `EMS` | the uniform fit, with PIPs re-weighted within each credible set by the expression modifier scores of Wang et al. (2021), from GTEx whole blood (microglia) or frontal cortex BA9 (other cell types) |
| `scEEMS_Weighted_Full` | scEEMS predictions of the Weighted (Full) model (step 6) |
| `scEEMS_Unweighted_Full` | predictions of the Unweighted (Full) model |
| `scEEMS_Weighted_Restricted` | predictions of the Weighted (Restricted) model |

The scEEMS predictions are recalibrated from the balanced training prior (0.5) to 1.8 / 14,716 per
variant by a logit shift and passed to SuSiE as prior weights; variants without a prediction get
1.8 / 14,716. SuSiE settings: L starts at 5 and grows by 5 (up to 30) while every effect forms a credible
set; scaled prior variance 0.2, prior and residual variances estimated; 95% credible sets with purity
(minimum absolute LD correlation) at least 0.5.

Each eQTL fit is then colocalized (`coloc.susie`) with the AD GWAS of Bellenguez et al. (2022), fine-mapped
over the same window with `susie_rss` on LD from the ADSP European reference panel, under two GWAS
priors: `uniform` and `polyfun`, per-SNP heritabilities estimated by PolyFun from the GWAS and 83
annotations (baseline-LF, microglia CRE, Roadmap E051 H3K27ac and two chromBPNet annotations). Variants are
matched across data sets by a canonical key (chromosome, position, alleles sorted).

`coloc.susie` uses the Bayes factors of the fits, which the EMS re-weighting does not change, so EMS
colocalizations equal the uniform ones; `compute_clpp.R` adds the CLPP (products of eQTL and GWAS PIPs),
which does respond to EMS.

## Inputs (config settings)

| Setting | Contents |
|---|---|
| `eqtl_data_dir` | ROSMAP eQTL data in the FunGen-xQTL layout: `genotype/ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.{N}.{bed,bim,fam}`; per cell type `{cell}/phenotype/` (region list, `phenotype_by_chrom/` expression BED files) and `{cell}/covariate/`; `reference/TADB_enhanced_cis.bed` (the cis window of each gene) |
| `adsp_plink_dir` | ADSP European reference panel, `ADSP_EUR_chr{N}.{bed,bim,fam}` |
| `gwas_sumstats_raw_file` | Bellenguez et al. (2022) stage 1, GWAS Catalog GCST90027158 (GRCh38) |
| `ldsc_dir`, `sumstats_file`, `ldsc_annotation_dir` | PolyFun installation, LDSC-munged GWAS and annotation LD scores, for the PolyFun prior |
| `ems_dir` | EMS tables `ems_top_{tissue}.tsv.bgz` (EMS release of Wang et al. 2021) |
| `predictions_dir` | step 6 per-gene predictions of the three models |
| `gene_mapping_file` | optional: `gene_id, gene_name, gene_TSS`, for gene names in `credset_all.tsv` |

Outputs go to `finemap_dir` (fits, coloc results, EMS vectors), `ld_cache_dir`, `gwas_cache_dir`,
`adsp_backing_dir`, `snpvar_dir` and `gwas_sumstats_file` (prepared GWAS), all under
`{output_dir}/fine_mapping/` by default, and tables to `{aggregate_dir}`.

## Scripts and order

| Order | Script | Job script | What |
|---|---|---|---|
| 1 | `prep_gwas_sumstats.py` | | GWAS summary statistics -> tabix table with z-scores |
| 1 | `compute_snpvar.py`, `snpvar_to_tabix.py` | `run_snpvar.sh` | PolyFun per-SNP priors for the GWAS |
| 1 | `build_backing.R` | `run_build_backing.sh` | reference panel -> bigsnpr backing files |
| 1 | `prepare_ems_priors.py CELL` | | EMS scores per gene |
| 2 | `make_gene_list.R` | | genes to fine-map in the GWAS (13,206 autosomal genes) |
| 3 | `precompute_ld.R` | `run_precompute_ld.sh` | LD of each gene's window |
| 4 | `precompute_gwas.R` | `run_precompute_gwas.sh` | GWAS fits (uniform and PolyFun) per gene |
| 5 | `finemap_and_coloc.R` | `run_finemap_coloc.sh` | eQTL fit (one cell type, prior and gene) + coloc |
| 6 | `compute_clpp.R`, `aggregate_clpp.py` | `run_clpp.sh` | CLPP |
| 6 | `aggregate_finemap.R`, `build_comparison.R`, `aggregate_coloc.R`, `build_coloc_table.R`, `compute_eqtl_marginal.R`, `build_credset_table.R`, `extract_credible_sets_tsv.py` | `run_aggregate.sh` | result tables and the release export |

`common.R` holds what the R scripts share (prior names, canonical keys, eQTL file locations and loading).
Every per-gene job skips genes whose output exists, so a resubmitted array reruns only failed genes.

```bash
cd 10_finemap_coloc
python prep_gwas_sumstats.py
sbatch run_snpvar.sh
sbatch run_build_backing.sh
for c in Mic Ast Exc Inh Oli OPC; do python prepare_ems_priors.py $c; done
Rscript make_gene_list.R                                    # prints the number of genes, N
sbatch --array=1-N%150 run_precompute_ld.sh
sbatch --array=1-N%200 run_precompute_gwas.sh               # after the LD jobs
for c in Mic Ast Exc Inh Oli OPC; do for p in uniform EMS scEEMS_Weighted_Full scEEMS_Unweighted_Full scEEMS_Weighted_Restricted; do for g in uniform polyfun; do
  sbatch --export=ALL,CELL=$c,PRIOR=$p,GWAS_PRIOR=$g --array=1-<regions of $c>%200 run_finemap_coloc.sh
done; done; done
sbatch run_clpp.sh                                          # when fine-mapping is done
python aggregate_clpp.py                                    # when the CLPP array has finished
sbatch run_aggregate.sh
```

The job scripts run in the `scEEMS_R` environment (`environment_r.yml` and `install_r_packages.R` at the
top of the repository: R 4.5.1, susieR 0.14.2, coloc 5.2.3, bigsnpr 1.12.21, pecotmr 0.3.16, htslib, and the
Python packages of the helper scripts). PolyFun runs in its own environment (`POLYFUN_ENV`, default
`polyfun`).

## Outputs

Under `{finemap_dir}/{cell}/`:
- `fine_mapping_{prior}/{gene}.{chr}.univariate_bvsr.rds` (SuSiE fit) and `.cs.tsv` (variant, PIP,
  credible set)
- `coloc_{prior}{tag}/{gene}.{chr}.coloc.tsv`, tag `""` (uniform GWAS prior) or `.polyfun`
- `gwas_susie_{prior}{tag}/`: the GWAS fit used, `ems_vector/`: EMS scores

Under `{aggregate_dir}/finemapping_comparison/`:
- `{cell}_{prior}_{gene,cs,cs_variants,pip010}.tsv`: per-gene, per-credible-set and per-variant tables
- `gene_level_comparison.tsv`, `gene_consensus.tsv`, `cs_matched_comparison.tsv`, `cs_match_summary.tsv`:
  the priors compared with the uniform prior (compare on genes with `in_all_priors`, which the scEEMS
  priors could fine-map)
- `{cell}_{prior}_bellenguez_{gwas_prior}_coloc.tsv`, `coloc_all.tsv`: every tested pair of credible sets
  with PP.H0-H4 and the GWAS p-values of the lead variants; colocalized = PP.H4 > 0.8
- `credset_all.tsv`: both credible sets of every colocalization, with GWAS and eQTL p-values

Under `{aggregate_dir}/clpp/`: `clpp_all.tsv`, `clpp_ems_gene_level.tsv`.

Under `{output_dir}/release/fine_mapping/`: the credible sets in the format of the data release.

# Step 11: Cross-Cell-Type Colocalization

**This step needs the eQTL SuSiE fits of step 10, which are computed from controlled-access ROSMAP
genotype and expression data. It cannot be run from the public data release; the code is provided so the
analysis can be inspected and rerun by those with data access.** The credible sets themselves are in the
`fine_mapping/` folder of the data release.

## Overview

Is an eQTL credible set specific to one cell type, or shared with another? For every gene, prior and cell
type with at least one credible set, `coloc_crosscell.R` runs `coloc.susie` between that cell type's
stored SuSiE fit and the fit of every other cell type with a credible set for the same gene under the same
prior. No model is refitted; variants are aligned on the canonical variant keys stored with each fit.

- **Shared credible set**: a credible set is shared if any credible set of the same gene in another cell
  type, under the same prior, colocalizes with it at PP.H4 > 0.8; otherwise it is not. The denominator is
  every credible set of that cell type under that prior, so a set whose gene has no credible set in any
  other cell type counts as not shared.
- **Shared eGene**: a gene with at least one credible set under the prior, at least one of which is shared.
- **Priors**: uniform, scEEMS Weighted (Full), Unweighted (Full) and Weighted (Restricted). EMS is not
  compared: it re-weights PIPs inside the uniform fit but leaves the credible sets and Bayes factors
  unchanged, so its results would equal uniform's.
- **Protein-coding genes only**: the scEEMS priors exist only for genes with scEEMS predictions, which are
  protein-coding, so a scEEMS prior is never fitted for a lncRNA while the uniform prior is. Without the
  restriction those genes would count as eGenes the scEEMS priors lost although the prior never ran there.

## Scripts

| Script | Description |
|--------|-------------|
| `make_gene_lists.py` | Genes with a step 10 fit (one per array task), and the GENCODE protein-coding gene list |
| `coloc_crosscell.R` | Cross-cell-type colocalization of one gene, all priors and cell types |
| `aggregate_crosscell.py` | Per-credible-set, per-eGene and per-cell-type sharing tables |
| `run_coloc_crosscell.sh` | SLURM array, one task per gene |
| `run_pipeline.sh` | Everything on one machine |

## Running

In the `scEEMS_R` environment (`environment_r.yml` at the top of the repository):

```bash
cd 11_crosscell_coloc
python make_gene_lists.py                                     # prints the array size, e.g. 1-13206
mkdir -p logs && sbatch --array=1-13206%200 run_coloc_crosscell.sh
python aggregate_crosscell.py                                 # after the array has finished
```

Each task takes well under a minute and about 330 MB of memory.

## Inputs

- `{finemap_dir}/{cell}/fine_mapping_{prior}/{gene_id}.{chr}.univariate_bvsr.rds`: the step 10 eQTL fits,
  for cell types Ast, Exc, Inh, Mic, Oli, OPC and priors `uniform`, `scEEMS_Weighted_Full`,
  `scEEMS_Unweighted_Full`, `scEEMS_Weighted_Restricted`
- `{aggregate_dir}/finemapping_comparison/{cell}_{prior}_cs.tsv` (optional): the step 10 credible-set tables, used to add
  credible-set size, top PIP and purity and to check the credible-set counts
- `gencode_gtf_file`: the GENCODE v45 basic annotation, `gencode.v45.basic.annotation.gtf.gz`, from
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/

## Outputs

Under `{crosscell_dir}` (default `{output_dir}/crosscell_coloc`):

- `genes.tsv` (gene_id, chr; no header) and `protein_coding_genes.txt`
- `per_gene/{gene_id}.{chr}.tsv`: one row per colocalization test, `gene_id, chr, prior, cell, idx,
  other_cell, other_idx, note, nsnps, PP.H0.abf ... PP.H4.abf`, where `idx` and `other_idx` are the
  credible-set indices of the two fits. Rows with a `note` record cases with nothing to test:
  `focal credible set` (one per credible set of `cell`; the denominator), `no fit` / `no credible set`
  (for `cell`), `no fit in other cell` / `no credible set in other cell`,
  `fewer than 10 shared variants`, `credible set has no shared variants`, `coloc returned no result`.

Under `{aggregate_dir}`, protein-coding genes only:

- `crosscell_coloc_pairs.tsv`: all rows of the per-gene tables
- `crosscell_cs_status.tsv`: one row per credible set (`prior, cell, gene_id, idx`) with `n_other_tested`
  (cell types it was tested against), `max_H4`, `shared_with` (cell types with PP.H4 > 0.8), `shared`
  (0/1) and, if the step 10 tables are present, `cs_size, min_abs_corr, top_pip`
- `crosscell_egene_status.tsv`: one row per (`prior, cell, gene_id`) with `n_cs`, `n_cs_shared`,
  `shared_with` and `shared` (0/1)
- `crosscell_summary.tsv`: per (`prior, cell`): credible sets and eGenes, the numbers shared, and
  `pct_cs_shared`, `pct_cs_shared_with_partner` (denominator: sets with at least one test) and
  `pct_egenes_shared`

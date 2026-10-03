#!/bin/bash
#SBATCH --job-name=sceems_aggregate_finemap
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4

# Tables of the fine-mapping and colocalization results, after every fine-mapping task has finished:
#   1. per (cell type, prior) fit tables and the comparison with the uniform prior
#   2. per (cell type, eQTL prior, GWAS prior) coloc tables, combined into coloc_all.tsv
#   3. credset_all.tsv, the marginal eQTL statistics of its genes, and credset_all.tsv again with them
#   4. the credible sets in the format of the data release
#   sbatch run_aggregate.sh
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS_R
cd "${SLURM_SUBMIT_DIR}"
CELLS="Mic Ast Exc Inh Oli OPC"
PRIORS="uniform EMS scEEMS_Weighted_Full scEEMS_Unweighted_Full scEEMS_Weighted_Restricted"

for c in $CELLS; do for p in $PRIORS; do
    Rscript --vanilla aggregate_finemap.R "$c" "$p"
done; done
Rscript --vanilla build_comparison.R

for c in $CELLS; do for p in $PRIORS; do for g in uniform polyfun; do
    Rscript --vanilla aggregate_coloc.R "$c" "$p" "$g"
done; done; done
Rscript --vanilla build_coloc_table.R

Rscript --vanilla build_credset_table.R
for c in $CELLS; do Rscript --vanilla compute_eqtl_marginal.R "$c"; done
Rscript --vanilla build_credset_table.R

python -u extract_credible_sets_tsv.py

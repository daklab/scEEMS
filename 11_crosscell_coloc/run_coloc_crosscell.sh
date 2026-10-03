#!/bin/bash
#SBATCH --job-name=sceems_crosscell
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=4G
#SBATCH --time=0:20:00
#SBATCH --cpus-per-task=1
#SBATCH --output=/dev/null
#SBATCH --error=logs/%x_%a.err

# Cross-cell-type colocalization, one array task per gene (a row of {crosscell_dir}/genes.tsv), with all
# priors and cell types inside the task (well under a minute and about 330 MB per gene). Run
# make_gene_lists.py first; it prints the array size:
#   mkdir -p logs && sbatch --array=1-13206%200 run_coloc_crosscell.sh
# A gene whose TSV exists is skipped and a TSV is written even when nothing could be tested, so
# resubmitting fills only the gaps; aggregate_crosscell.py reports genes without a TSV.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS_R
cd "${SLURM_SUBMIT_DIR}"

Rscript --vanilla coloc_crosscell.R "${SLURM_ARRAY_TASK_ID}"

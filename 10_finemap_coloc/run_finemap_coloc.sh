#!/bin/bash
#SBATCH --job-name=sceems_finemap
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=24G
#SBATCH --time=0:15:00
#SBATCH --cpus-per-task=1

# eQTL fine-mapping and colocalization for one cell type, eQTL prior and GWAS prior, one array task per
# row of the cell type's region list (its line count minus the header). Run after run_precompute_gwas.sh;
# the EMS prior also needs prepare_ems_priors.py and the scEEMS priors the step 6 predictions.
#   N=$(( $(wc -l < <eqtl_data_dir>/Mic/phenotype/snuc_pseudo_bulk.Mic.mega.normalized.log2cpm.region_list.txt) - 1 ))
#   sbatch --export=ALL,CELL=Mic,PRIOR=scEEMS_Weighted_Full,GWAS_PRIOR=polyfun --array=1-${N}%200 run_finemap_coloc.sh
# PRIOR: uniform | EMS | scEEMS_Weighted_Full | scEEMS_Unweighted_Full | scEEMS_Weighted_Restricted
# GWAS_PRIOR: uniform | polyfun. The eQTL fit is shared by the two GWAS priors. Finished genes are
# skipped, so resubmitting reruns only failed tasks.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS_R
cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

Rscript --vanilla finemap_and_coloc.R "${SLURM_ARRAY_TASK_ID}" "${CELL}" "${PRIOR}" "${GWAS_PRIOR:-uniform}"

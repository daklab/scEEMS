#!/bin/bash
#SBATCH --job-name=baseline_ldscore                  # Job name
#SBATCH --partition=cpu                      # Partition Name
#SBATCH --mail-type=FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=user@example.com         # Where to send mail
#SBATCH --mem=50G                            # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --time=40:00:00                       # Time limit
#SBATCH --array=1-22%22

# Activate the environment
conda activate scEEMS

# Set paths from config.yaml (adjust these to match your config)
plink_dir=${PLINK_REF_DIR}  # paths.plink_ref_dir in config.yaml
annot_dir=${ANNOT_DIR}      # paths.output_dir/MLxQTL_pareto_nocontrol/${cell_type}
polyfun_dir=${POLYFUN_DIR}  # paths.ldsc_dir in config.yaml

cd ${polyfun_dir}

python compute_ldscores.py	\
  --annot ${annot_dir}/MLxQTL_chr${SLURM_ARRAY_TASK_ID}.annot.gz \
  --bfile ${plink_dir}/ADSP_chr${SLURM_ARRAY_TASK_ID} \
  --out ${annot_dir}/MLxQTL_chr${SLURM_ARRAY_TASK_ID}.l2.ldscore.parquet \
  --allow-missing \
  --ld-wind-kb 1000

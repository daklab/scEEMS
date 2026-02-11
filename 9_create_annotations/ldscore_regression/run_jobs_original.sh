#!/bin/bash
#SBATCH --job-name=ldscore                  # Job name
#SBATCH --partition=cpu                      # Partition Name
#SBATCH --mail-type=FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=user@example.com         # Where to send mail
#SBATCH --mem=150G                            # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --time=0:30:00                       # Time limit
#SBATCH --array=1-19%19

# Activate the environment
conda activate scEEMS

# Run LDSC regression (original weighted model)
python ldscore_template_original.py ${SLURM_ARRAY_TASK_ID} ${cell_type}
# Usage: sbatch --export=cell_type=Mic_mega_eQTL run_jobs_original.sh

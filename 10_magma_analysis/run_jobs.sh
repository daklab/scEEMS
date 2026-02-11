#!/bin/bash
#SBATCH --job-name=magma           # Job name
#SBATCH --partition=cpu                      # Partition Name
#SBATCH --mail-type=FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=user@example.com         # Where to send mail
#SBATCH --mem=50G                            # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --time=1:50:00                       # Time limit
#SBATCH --cpus-per-task=5               # Standard output and error log
#SBATCH --array=1-22%22

# Activate the environment
conda activate scEEMS

# Run the Python script (revised MAGMA files)
python -u make_magma_files_revised.py ${SLURM_ARRAY_TASK_ID} ${cell_type}
# Usage: sbatch --export=cell_type=Mic_mega_eQTL run_jobs.sh

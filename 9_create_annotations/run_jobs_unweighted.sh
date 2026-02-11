#!/bin/bash
#SBATCH --job-name=annotations_unweighted          # Job name
#SBATCH --partition=cpu                      # Partition Name
#SBATCH --mail-type=FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=user@example.com         # Where to send mail
#SBATCH --mem=50G                            # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --time=1:20:00                       # Time limit
#SBATCH --cpus-per-task=5               # Standard output and error log
#SBATCH --array=1-22%22

# Activate the environment
conda activate scEEMS

# Run the Python script (original unweighted model)
python -u make_annotations_original_unweighted.py ${SLURM_ARRAY_TASK_ID} ${cell_type}
# Usage: sbatch --export=cell_type=Mic_mega_eQTL run_jobs_unweighted.sh

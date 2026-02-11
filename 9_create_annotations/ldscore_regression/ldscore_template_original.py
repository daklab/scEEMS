import os
import csv, subprocess
import sys
import pandas as pd
import itertools
import yaml

# Load config
with open(os.path.join(os.path.dirname(__file__), '..', '..', 'config.yaml'), 'r') as f:
    config = yaml.safe_load(f)

# Change to the PolyFun directory
os.chdir(config['paths']['ldsc_dir'])


idx = int(sys.argv[1])

idx = idx - 1

cell_type = sys.argv[2]

# Base configuration
cohort_id = 'bellenguez_omics_dl'
meta_cols = ['CHR', 'SNP', 'BP', 'A1', 'A2']

# Important paths
path_sumstats = config['paths']['sumstats_file']
weight = os.path.join(config['paths']['weights_dir'], 'weights_chr')

# Annotation paths
bl = os.path.join(config['paths']['baseline_annot_dir'], 'baseline_chr')
MLxQTL = os.path.join(config['paths']['output_dir'], 'MLxQTL_pareto', f'{cell_type}_original', 'MLxQTL_chr')

# Read baseline annotations
bl_file = f'{bl}22.annot.gz'
bl_string = pd.read_csv(bl_file, sep='\t').columns.tolist()
bl_string = [col for col in bl_string if col not in meta_cols]

# Read glass lab annotations
MLxQTL_file = f'{MLxQTL}22.annot.gz'
MLxQTL_string = pd.read_csv(MLxQTL_file, sep='\t').columns.tolist()
MLxQTL_string = [col for col in MLxQTL_string if col not in meta_cols]

# Combine all selected annotations
annot_string_array = bl_string + [MLxQTL_string[idx]]
# Remove any columns that should be excluded

annot_string = ','.join(annot_string_array)

# Set up output paths
path_out = config['paths']['ldsc_results_dir']
out_bellenguez = f'{path_out}/MLxQTL_pareto_original/bellenguez_2022/{cell_type}'


file_out = f'{out_bellenguez}/{cell_type}_{cohort_id}_subset_{idx}'

# Create output directory if it doesn't exist
os.makedirs(out_bellenguez, exist_ok=True)

# Set up annotation files string
annot_files = ','.join([bl, MLxQTL])

# Create the submit command
submit_command = (f'python ldsc.py '
                  f'--h2 {path_sumstats} '
                  f'--ref-ld-chr {annot_files} '
                  f'--w-ld-chr {weight} '
                  f' --anno {annot_string} '
                  f'--out {file_out} '
                  f'--overlap-annot '
                  f'--not-M-5-50 '
                  f'--print-coefficients'
                  f' --print-delete-vals')

print(f"Running analysis for subset {idx}")
print("Selected annotations:")
print("\nSubmit command:")
print(submit_command)

# Execute the command
exit_status = subprocess.call(submit_command, shell=True)
if exit_status == 1:
    print(f"Job for subset {idx} failed to submit")
"""
Compare PIPs between uniform and scEEMS-informed fine-mapping.

Compares posterior inclusion probabilities from SuSiE runs with uniform
priors vs. scEEMS prediction-based priors across cell types.

Usage:
    python compare_pips.py

Requires config.yaml with paths configured.
"""

import os
import pandas as pd
from dask import dataframe as dd
from dask.diagnostics import ProgressBar
import yaml

ProgressBar().register()

# Load configuration
with open('../config.yaml', 'r') as f:
    config = yaml.safe_load(f)

data_dir = config['paths']['data_dir']
fine_mapping_dir = config['paths']['fine_mapping_dir']
out_dir = os.path.join(data_dir, 'training_data/aggregate_results/finemapping_comparison')
os.makedirs(out_dir, exist_ok=True)

cell_types = ['Mic', 'Ast']

for cell_type in cell_types:
    print(f"Processing {cell_type}...")

    # Load uniform fine-mapping results
    uniform_dir = os.path.join(fine_mapping_dir, cell_type, 'fine_mapping_top_uniform')
    uniform_df = dd.read_csv(
        f'{uniform_dir}/*univariate_bvsr_top_loci.tsv', sep='\t').compute()

    # Load predicted fine-mapping results
    predicted_dir = os.path.join(fine_mapping_dir, cell_type, 'fine_mapping_top_predicted')
    predicted_df = dd.read_csv(
        f'{predicted_dir}/*univariate_bvsr_top_loci.tsv', sep='\t').compute()

    # Parse variant IDs
    for df in [uniform_df, predicted_df]:
        df[['chr_val', 'pos', 'ref', 'alt']] = df['variant_id'].str.split(':', expand=True)
        df.drop(columns=['chr_val'], inplace=True)

    # Save comparison files
    uniform_df.to_csv(os.path.join(out_dir, f'{cell_type}_uniform_finemapping.tsv'),
                      sep='\t', index=False)
    predicted_df.to_csv(os.path.join(out_dir, f'{cell_type}_predictions_finemapping.tsv'),
                        sep='\t', index=False)

    print(f"  Uniform: {len(uniform_df)} variants in credible sets")
    print(f"  Predicted: {len(predicted_df)} variants in credible sets")

print("Comparison files saved")

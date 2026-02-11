"""
Aggregate per-chromosome SHAP values into combined Parquet files.

Usage:
    python aggregate_shap_scores.py <cohort>

Arguments:
    cohort: Cell type cohort (e.g., Mic_mega_eQTL)
"""

import os
import sys
import pandas as pd
from dask import dataframe as dd
from dask.diagnostics import ProgressBar
import yaml

ProgressBar().register()

# Load configuration
with open('../config.yaml', 'r') as f:
    config = yaml.safe_load(f)

data_dir = config['paths']['data_dir']
cohort = sys.argv[1]

training_data_dir = os.path.join(data_dir, f'training_data/{cohort}')
shap_dir = os.path.join(training_data_dir, 'model_results/shap_values')
parquet_dir = os.path.join(shap_dir, 'parquet')
os.makedirs(parquet_dir, exist_ok=True)

enhancer_shap = dd.read_csv(f'{shap_dir}/enhancer_shap_values_chr*.csv', sep=',').compute()
promoter_shap = dd.read_csv(f'{shap_dir}/promoter_shap_values_chr*.csv', sep=',').compute()

print(f'{cohort} enhancer shap shape: {enhancer_shap.shape}')
print(f'{cohort} promoter shap shape: {promoter_shap.shape}')

enhancer_parquet_file = os.path.join(parquet_dir, 'enhancer_shap_values.parquet')
promoter_parquet_file = os.path.join(parquet_dir, 'promoter_shap_values.parquet')

if os.path.exists(enhancer_parquet_file):
    os.remove(enhancer_parquet_file)
if os.path.exists(promoter_parquet_file):
    os.remove(promoter_parquet_file)

enhancer_shap.to_parquet(enhancer_parquet_file, index=False)
promoter_shap.to_parquet(promoter_parquet_file, index=False)

print(f"Aggregated SHAP values saved to {parquet_dir}")

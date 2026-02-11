"""
Aggregate per-gene predictions into chromosome-level Parquet files.

Combines individual gene prediction TSV files into a single Parquet dataset
partitioned by chromosome for efficient downstream access.

Usage:
    python aggregate_predictions.py <cohort>

Arguments:
    cohort: Cell type cohort (e.g., Mic_mega_eQTL)

Requires config.yaml with paths configured.
"""

import os
import sys
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
prediction_dir = os.path.join(training_data_dir, 'predictions_catboost')
prediction_dir_parquet = os.path.join(training_data_dir, 'predictions_parquet_catboost')
predictions_parquet_file = os.path.join(prediction_dir_parquet, 'predictions.parquet')

os.makedirs(prediction_dir_parquet, exist_ok=True)

# Read all per-gene prediction files
predictions_df = dd.read_csv(f'{prediction_dir}/*_predictions.tsv', sep='\t')
predictions_df = predictions_df.repartition(partition_size='100MB')
predictions_df = predictions_df.persist()

# Write to Parquet partitioned by chromosome
predictions_df.to_parquet(predictions_parquet_file, engine='pyarrow',
                          compression='snappy', write_index=False,
                          partition_on=['chr'], overwrite=True)

print(f"Aggregated predictions saved to {predictions_parquet_file}")

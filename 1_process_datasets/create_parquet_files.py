"""
Convert the per-gene CSV tables of get_vars_pips.R (PIP_all and PIP_top) of every eQTL analysis into
parquet datasets partitioned by chromosome.

Usage:
    python create_parquet_files.py

Input:   {susie_dir}/<analysis>/PIP_all/*.csv and PIP_top/*.csv
Output:  {susie_dir}/<analysis>/PIP_all_parquet/PIP_all.parquet and PIP_top_parquet/PIP_top.parquet
"""

import os
import sys

from dask import dataframe as dd
from dask.diagnostics import ProgressBar

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

ProgressBar().register()

data_directory = path("susie_dir")

# One folder per eQTL analysis written by get_vars_pips.R; other folders (e.g. variant_list) are skipped
subdirectories = [f.name for f in os.scandir(data_directory)
                  if f.is_dir() and os.path.isdir(os.path.join(f.path, 'PIP_all'))]

# Convert PIP_all CSVs to Parquet
for subdirectory in subdirectories:
    print(subdirectory)
    path = f'{data_directory}/{subdirectory}/PIP_all'
    try:
        df = dd.read_csv(path + '/*.csv', sep=' ', dtype={'pip': 'float64'})
        write_path = f'{data_directory}/{subdirectory}/PIP_all_parquet'
        os.makedirs(write_path, exist_ok=True)
        df = df.repartition(partition_size="100MB")
        df.to_parquet(f'{write_path}/PIP_all.parquet', engine='pyarrow',
                      compression='snappy', partition_on='chr', overwrite=True)
    except ValueError as e:
        if "Length mismatch" in str(e):
            print(f"Skipping subdirectory {subdirectory} due to error: {str(e)}")
        else:
            raise e

# Column type mapping for PIP_top files
dtypes_mapping = {
    'chr': 'string',
    'pos': 'int64',
    'ref': 'string',
    'alt': 'string',
    'betahat': 'double',
    'sebetahat': 'double',
    'maf': 'double',
    'pip': 'float64',
    'cs_coverage_0.95': 'int64',
    'cs_coverage_0.7': 'int64',
    'cs_coverage_0.5': 'int64',
    'cs_coverage_0.95_min_corr': 'int64',
    'cs_coverage_0.7_min_corr': 'int64',
    'cs_coverage_0.5_min_corr': 'int64',
    'gene_id': 'string'
}

# Convert PIP_top CSVs to Parquet
for subdirectory in subdirectories:
    print(subdirectory)
    path = f'{data_directory}/{subdirectory}/PIP_top'
    try:
        df = dd.read_csv(path + '/*.csv', sep=' ', dtype=dtypes_mapping)
        write_path = f'{data_directory}/{subdirectory}/PIP_top_parquet'
        os.makedirs(write_path, exist_ok=True)
        df = df.repartition(partition_size="100MB")
        df.to_parquet(f'{write_path}/PIP_top.parquet', engine='pyarrow',
                      compression='snappy', partition_on='chr', overwrite=True)
    except ValueError as e:
        if "Length mismatch" in str(e):
            print(f"Skipping subdirectory {subdirectory} due to error: {str(e)}")
        else:
            raise e

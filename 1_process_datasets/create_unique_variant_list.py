"""
Create a deduplicated variant list across all cohorts for a given chromosome.

Reads PIP_all parquet files from all cohorts and produces a unique variant list
per chromosome for downstream annotation.

Usage:
    python create_unique_variant_list.py <chromosome_number>

Arguments:
    chromosome_number: Integer chromosome number (1-22)

Requires config.yaml with paths.data_dir set.
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

data_directory = os.path.join(config['paths']['data_dir'], 'susie_vars_pips')
variant_dir = os.path.join(data_directory, 'variant_list')
os.makedirs(variant_dir, exist_ok=True)

# Get all subdirectory names (cohorts)
subdirectories = [f.name for f in os.scandir(data_directory) if f.is_dir()]
subdirectories = [x for x in subdirectories if x != 'variant_list']

i = sys.argv[1]
print(i)

dfs = []
for subdirectory in subdirectories:
    print(subdirectory)
    chromosome_val = f'chr{i}'
    parquet_file = f'{data_directory}/{subdirectory}/PIP_all_parquet/PIP_all.parquet'
    df = dd.read_parquet(parquet_file, engine='pyarrow')
    df = df[df['chr'] == chromosome_val]
    df = df[['variant_id', 'chr', 'pos', 'ref', 'alt']].drop_duplicates()
    dfs.append(df)

merged_df = dd.concat(dfs)
unique_df = merged_df.drop_duplicates()
unique_df = unique_df.repartition(partition_size="100MB")
unique_df.to_parquet(f'{variant_dir}/variant_list_chr{i}.parquet',
                     engine='pyarrow', compression='snappy', overwrite=True)

print(f"Saved unique variant list for chr{i}")

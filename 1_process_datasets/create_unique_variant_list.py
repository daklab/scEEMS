"""
List the unique variants of one chromosome over the fine-mapping results of all eQTL analyses: the
variants step 2 annotates.

Usage:
    python create_unique_variant_list.py <chromosome_number>        (1-22)

Input:   {susie_dir}/<analysis>/PIP_all_parquet/PIP_all.parquet (create_parquet_files.py)
Output:  {variant_list_dir}/variant_list_chr{N}.parquet
"""

import os
import sys

from dask import dataframe as dd
from dask.diagnostics import ProgressBar

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

ProgressBar().register()

data_directory = path("susie_dir")
variant_dir = path("variant_list_dir")
os.makedirs(variant_dir, exist_ok=True)

# One folder per eQTL analysis with a PIP_all parquet dataset; other folders (e.g. variant_list) are skipped
subdirectories = [f.name for f in os.scandir(data_directory)
                  if f.is_dir() and os.path.isdir(os.path.join(f.path, 'PIP_all_parquet', 'PIP_all.parquet'))]

i = sys.argv[1]
print(i)

dfs = []
for subdirectory in subdirectories:
    print(subdirectory)
    chromosome_val = f'chr{i}'
    parquet_file = f'{data_directory}/{subdirectory}/PIP_all_parquet/PIP_all.parquet'
    # read only this chromosome's partition (the dataset is partitioned by chr)
    df = dd.read_parquet(parquet_file, engine='pyarrow', filters=[('chr', '==', chromosome_val)])
    df = df[df['chr'] == chromosome_val]
    df = df[['variant_id', 'chr', 'pos', 'ref', 'alt']].drop_duplicates()
    dfs.append(df)

merged_df = dd.concat(dfs)
unique_df = merged_df.drop_duplicates()
unique_df = unique_df.repartition(partition_size="100MB")
unique_df.to_parquet(f'{variant_dir}/variant_list_chr{i}.parquet',
                     engine='pyarrow', compression='snappy', overwrite=True)

print(f"Saved unique variant list for chr{i}")

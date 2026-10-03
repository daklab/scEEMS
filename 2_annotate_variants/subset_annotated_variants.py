#!/usr/bin/env python
"""
Keep only the identifier columns (CHR, BP, REF, ALT, SNP; index variant_id) of one chromosome's annotated
variants (annotate_variants.py). Steps 3 and 4 join the fine-mapped variants to this small table to keep
the variants that have annotations, instead of reading the full annotation table.

Usage:
    python subset_annotated_variants.py <chromosome_number>        (1-22)

Input:   {variant_list_dir}/annotated_variants/annotated_variants_chr{N}.parquet
Output:  {variant_list_dir}/annotated_variants_just_variants/annotated_just_variants_chr{N}.parquet
"""
import os
import sys

import dask
from dask import dataframe as dd
from dask.diagnostics import ProgressBar

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

ProgressBar().register()
dask.config.set({'temporary_directory': path("scratch_dir")})

chr_num = sys.argv[1]
variant_list = path("variant_list_dir")
annotation_file = f'{variant_list}/annotated_variants/annotated_variants_chr{chr_num}.parquet'

annotation_df = dd.read_parquet(annotation_file, engine='pyarrow', columns=['CHR', 'BP', 'REF', 'ALT', 'SNP'])
annotation_df = annotation_df.repartition(partition_size="100MB")

write_path = f'{variant_list}/annotated_variants_just_variants'
os.makedirs(write_path, exist_ok=True)
write_file = f'{write_path}/annotated_just_variants_chr{chr_num}.parquet'
annotation_df.to_parquet(write_file, engine='pyarrow', compression='snappy', overwrite=True)
print(f"Saved the variant identifiers of chr{chr_num}")

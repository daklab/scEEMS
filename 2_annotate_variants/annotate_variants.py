"""
Annotate the unique variants of one chromosome (step 1) with the variant-level features.

Merges the variant list with:
- Enformer variant effect predictions (CAGE tracks dropped)
- brain cell type ATAC-seq, promoter and enhancer BED annotations (bed_annotations_dir)
- baseline genomic annotations (baseline_annotations_dir)
- composite transcription factor scores: the maximum Enformer ChIP score over each cell type's TFs, and
  its product with the cell type's promoter/enhancer/ATAC annotation

Usage:
    python annotate_variants.py <chromosome_number>        (1-22)

Inputs:  {variant_list_dir}/variant_list_chr{N}.parquet, {enformer_dir}/enformer_tensorflow_chr{N}.parquet,
         the BED files, tf_file and targets_file (see README.md)
Output:  {variant_list_dir}/annotated_variants/annotated_variants_chr{N}.parquet
"""

import os
import sys
import pandas as pd
from dask import dataframe as dd
from dask.diagnostics import ProgressBar
from tqdm import tqdm
import numpy as np
from pybedtools import BedTool
import dask

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

ProgressBar().register()

scratch_dir = path("scratch_dir")
dask.config.set({'temporary_directory': scratch_dir})


##############################################################################################################
# Functions

def anti_join(x, y, on):
    """Return rows in x which are not present in y."""
    ans = pd.merge(left=x, right=y, how='left', indicator=True, on=on)
    ans = ans.loc[ans._merge == 'left_only', :].drop(columns='_merge')
    return ans


def merge_beds(snp_list_df, bed_list, bed_directory):
    """Intersect variant positions with BED file annotations."""
    iter_bim = [[str(x1), str(x2 - 1), str(x2)]
                for (x1, x2) in np.array(snp_list_df[['CHR', 'BP']])]
    full_df = snp_list_df
    bimbed = BedTool(iter_bim)
    for col_name, file_name in bed_list.items():
        bed_for_annot = BedTool(os.path.join(bed_directory, file_name))
        annotbed = bimbed.intersect(bed_for_annot)
        bp = [x.start + 1 for x in annotbed]
        df_int = pd.DataFrame({'BP': bp, 'ANNOT': 1})
        df_int_unique = df_int.groupby('BP').first().reset_index()
        df_annot = pd.merge(snp_list_df, df_int_unique, how='left', on='BP')
        df_annot.fillna(0, inplace=True)
        df_annot[['ANNOT']] = df_annot[['ANNOT']].astype(int)
        df_annot = df_annot[['BP', 'ANNOT']]
        df_annot = df_annot.rename(columns={'ANNOT': col_name})
        df_annot_unique = df_annot.groupby('BP').first().reset_index()
        full_df = pd.merge(full_df, df_annot_unique, how='left', on='BP')
    return full_df


def merge_beds_baseline(snp_list_df, bed_list, bed_directory):
    """Intersect variant positions with baseline BED annotations."""
    iter_bim = [['chr' + str(x1), str(x2 - 1), str(x2)]
                for (x1, x2) in np.array(snp_list_df[['CHR', 'BP']])]
    full_df = snp_list_df
    bimbed = BedTool(iter_bim)
    for row in tqdm(bed_list):
        col_name = row['column_name']
        file_name = row['file_name']
        bed_for_annot = BedTool(f'{bed_directory}/{file_name}')
        annotbed = bimbed.intersect(bed_for_annot)
        bp = [x.start + 1 for x in annotbed]
        df_int = pd.DataFrame({'BP': bp, 'ANNOT': 1})
        df_int_unique = df_int.groupby('BP').first().reset_index()
        df_annot = pd.merge(snp_list_df, df_int_unique, how='left', on='BP')
        df_annot.fillna(0, inplace=True)
        df_annot[['ANNOT']] = df_annot[['ANNOT']].astype(int)
        df_annot = df_annot[['BP', 'ANNOT']]
        df_annot = df_annot.rename(columns={'ANNOT': col_name})
        df_annot_unique = df_annot.groupby('BP').first().reset_index()
        full_df = pd.merge(full_df, df_annot_unique, how='left', on='BP')
    return full_df


def multiply_row(row, col1, col2):
    return int(row[col1] * row[col2])


##############################################################################################################
# TF list setup

TF_file = path("tf_file")
TF_df = pd.read_csv(TF_file, sep='\t', skiprows=1)

targets_txt = path("targets_file")
df_targets = pd.read_csv(targets_txt, sep='\t')

df_targets_cage = df_targets[df_targets['description'].str.contains("CAGE")]
df_targets_subset = df_targets[df_targets['sum_stat'].str.contains("mean")]
df_targets_subset_TF = df_targets_subset[df_targets_subset['description'].str.contains("CHIP:")]
df_targets_subset_TF[['CHIP', 'TF', 'cell_type']] = df_targets_subset_TF.description.str.split(":", expand=True)

TF_microglia = ['ETS2', 'FLI1', 'SPI1', 'IRF8', 'PRDM1', 'CEBPA', 'CEBPB', 'CEBPD', 'CEBPE',
                'KLF11', 'KLF2', 'TFEC', 'BACH1', 'BATF', 'BATF3', 'RUNX1', 'RUNX2', 'RUNX3',
                'LYL1', 'TAL1', 'ELK3', 'MEF2C', 'CREB3L2', 'ATF4', 'MAFB', 'SALL1', 'eGFP-SALL1',
                'SMAD5', 'MEF2A', 'MEF2B', 'SMAD2', 'USF1', 'STAT3', 'eGFP-MAFG', 'NFYB', 'NRF1',
                'CREB1', 'IRF1', 'SOX9']

TF_neuronal = ['ASCL1', 'ASCL5', 'BHLHA15', 'BHLHE22', 'NEUROD1', 'NEUROD2', 'NEUROD6',
               'TWIST2', 'EGR4', 'SP8', 'SP9', 'MEIS3', 'PKNOX2', 'MSC', 'KLF5', 'KLF8',
               'TBR1', 'HLF']

TF_astrocyte = ['MEIS2', 'NFIB', 'TGIF1', 'EMX2', 'LHX2', 'RFX2', 'RFX4', 'RORA', 'RORB',
                'SOX1', 'SOX2', 'SOX21', 'SOX9', 'NFATC4', 'SOX5', 'POU3F2', 'POU3F3',
                'POU3F4', 'FOXG1', 'FOXO1', 'SP5']

TF_oligodendrocyte = ['SOX10', 'SOX13', 'SOX21', 'SOX3', 'SOX6', 'SOX8', 'NHLH2', 'NFIX',
                      'SP7', 'NFE2', 'E2F1', 'CREB5', 'POU3F3', 'MYCN']

microglia_TF_columns = df_targets_subset_TF[df_targets_subset_TF['TF'].isin(TF_microglia)].identifier.tolist()
neuronal_TF_columns = df_targets_subset_TF[df_targets_subset_TF['TF'].isin(TF_neuronal)].identifier.tolist()
astrocyte_TF_columns = df_targets_subset_TF[df_targets_subset_TF['TF'].isin(TF_astrocyte)].identifier.tolist()
oligodendrocyte_TF_columns = df_targets_subset_TF[df_targets_subset_TF['TF'].isin(TF_oligodendrocyte)].identifier.tolist()

TF_all_list = TF_df['Target of assay'].tolist() + TF_microglia + TF_neuronal + TF_astrocyte + TF_oligodendrocyte
TF_all_list = list(set(TF_all_list))
all_TF_columns = df_targets_subset_TF[df_targets_subset_TF['TF'].isin(TF_all_list)].identifier.tolist()

microglia_TF_columns_diff_32 = [f'diff_32_{i}' for i in microglia_TF_columns]
neuronal_TF_columns_diff_32 = [f'diff_32_{i}' for i in neuronal_TF_columns]
astrocyte_TF_columns_diff_32 = [f'diff_32_{i}' for i in astrocyte_TF_columns]
oligodendrocyte_TF_columns_diff_32 = [f'diff_32_{i}' for i in oligodendrocyte_TF_columns]

cage_columns = df_targets_cage.identifier.tolist()
cage_columns_diff_32 = [f'diff_32_{i}' for i in cage_columns]

all_TF_columns_diff_32 = list(set([f'diff_32_{i}' for i in all_TF_columns]))

meta_columns = ['CHR', 'SNP', 'BP', 'REF', 'ALT']
columns_all_TF = meta_columns + all_TF_columns_diff_32

##############################################################################################################
# Brain ATAC-seq/histone BED annotations

beds_dir = path("bed_annotations_dir")

bed_list = {
    'astrocyte_atac': 'LHX2_optimal_peak_IDR_ENCODE.ATAC.bed',
    'neuron_atac': 'NeuN_optimal_peak_IDR_ENCODE.ATAC.bed',
    'oligodendrocyte_atac': 'Olig2_optimal_peak_IDR_ENCODE.ATAC.bed',
    'microglia_atac': 'PU1_optimal_peak_IDR_ENCODE.ATAC.bed',
    'astrocyte_promoter': 'LHX2_promoter.bed',
    'neuron_promoter': 'NeuN_promoter.bed',
    'oligodendrocyte_promoter': 'Olig2_promoter.bed',
    'microglia_promoter': 'PU1_promoter.bed',
    'astrocyte_enhancer': 'LHX2_enhancer.bed',
    'neuron_enhancer': 'NeuN_enhancer.bed',
    'oligodendrocyte_enhancer': 'Olig2_enhancer.bed',
    'microglia_enhancer': 'PU1_enhancer.bed',
    'astrocyte_promoter_intersect_atac': 'LHX2_promoter_atac.bed',
    'neuron_promoter_intersect_atac': 'NeuN_promoter_atac.bed',
    'oligodendrocyte_promoter_intersect_atac': 'Olig2_promoter_atac.bed',
    'microglia_promoter_intersect_atac': 'PU1_promoter_atac.bed',
    'astrocyte_enhancer_intersect_atac': 'LHX2_enhancer_atac.bed',
    'neuron_enhancer_intersect_atac': 'NeuN_enhancer_atac.bed',
    'oligodendrocyte_enhancer_intersect_atac': 'Olig2_enhancer_atac.bed',
    'microglia_enhancer_intersect_atac': 'PU1_enhancer_atac.bed',
    'astrocyte_enhancer_promoter': 'LHX2_promoter_enhancer.bed',
    'neuron_enhancer_promoter': 'NeuN_promoter_enhancer.bed',
    'oligodendrocyte_enhancer_promoter': 'Olig2_promoter_enhancer.bed',
    'microglia_enhancer_promoter': 'PU1_promoter_enhancer.bed',
    'astrocyte_enhancer_promoter_intersect_atac': 'LHX2_promoter_enhancer_atac.bed',
    'neuron_enhancer_promoter_intersect_atac': 'NeuN_promoter_enhancer_atac.bed',
    'oligodendrocyte_enhancer_promoter_intersect_atac': 'Olig2_promoter_enhancer_atac.bed',
    'microglia_enhancer_promoter_intersect_atac': 'PU1_promoter_enhancer_atac.bed',
    'astrocyte_enhancer_promoter_union_atac': 'LHX2_promoter_enhancer_plus_atac.bed',
    'neuron_enhancer_promoter_union_atac': 'NeuN_promoter_enhancer_plus_atac.bed',
    'oligodendrocyte_enhancer_promoter_union_atac': 'Olig2_promoter_enhancer_plus_atac.bed',
    'microglia_enhancer_promoter_union_atac': 'PU1_promoter_enhancer_plus_atac.bed',
    'astrocyte_enhancer_promoter_500': 'LHX2_promoter_enhancer_500.bed',
    'neuron_enhancer_promoter_500': 'NeuN_promoter_enhancer_500.bed',
    'oligodendrocyte_enhancer_promoter_500': 'Olig2_promoter_enhancer_500.bed',
    'microglia_enhancer_promoter_500': 'PU1_promoter_enhancer_500.bed',
    'astrocyte_enhancer_promoter_union_atac_500': 'LHX2_promoter_enhancer_plus_atac_500.bed',
    'neuron_enhancer_promoter_union_atac_500': 'NeuN_promoter_enhancer_plus_atac_500.bed',
    'oligodendrocyte_enhancer_promoter_union_atac_500': 'Olig2_promoter_enhancer_plus_atac_500.bed',
    'microglia_enhancer_promoter_union_atac_500': 'PU1_promoter_enhancer_plus_atac_500.bed',
}

##############################################################################################################
# Main annotation logic

enformer_dir = path("enformer_dir")
variant_list = path("variant_list_dir")

chr_num = sys.argv[1]

parquet_file = f'{variant_list}/variant_list_chr{chr_num}.parquet'
enformer_file = f'{enformer_dir}/enformer_tensorflow_chr{chr_num}.parquet'

variant_df = dd.read_parquet(parquet_file, engine='pyarrow')

# Remap column names
remapping_dict = {'pos': 'BP', 'alt': 'ALT', 'ref': 'REF', 'chr': 'CHR'}
variant_df = variant_df.rename(columns=remapping_dict)
variant_df['CHR'] = variant_df['CHR'].astype(str)

enformer_df = dd.read_parquet(enformer_file)

# Remove CAGE columns from Enformer
enformer_df = enformer_df.drop(cage_columns_diff_32, axis=1)

# Merge variant list with Enformer scores
annotation_df = variant_df.merge(enformer_df, on=['CHR', 'BP', 'REF', 'ALT'], how='inner')
annotation_df = annotation_df.compute()

##############################################################################################################
# BED annotations

annotation_df_variants_only = annotation_df[['CHR', 'BP', 'REF', 'ALT', 'SNP']]
brain_atac_df = merge_beds(annotation_df_variants_only, bed_list, beds_dir)

##############################################################################################################
# Baseline annotations

baseline_beds_dir = path("baseline_annotations_dir")
baseline_files = [f for f in os.listdir(baseline_beds_dir) if not f.endswith('.unmapped.bed')]
baseline_files_dict = [{'column_name': f.replace('.bed', ''), 'file_name': f} for f in baseline_files]
baseline_df = merge_beds_baseline(annotation_df_variants_only, baseline_files_dict, baseline_beds_dir)

##############################################################################################################
# TF score processing

df_merged = annotation_df[columns_all_TF]
df_merged = df_merged.merge(brain_atac_df, on=meta_columns, how='inner')

df_merged['none-none-all'] = df_merged[all_TF_columns_diff_32].max(axis=1).astype('int')
df_merged['none-none-microglia'] = df_merged[microglia_TF_columns_diff_32].max(axis=1).astype('int')
df_merged['none-none-neuron'] = df_merged[neuronal_TF_columns_diff_32].max(axis=1).astype('int')
df_merged['none-none-astrocyte'] = df_merged[astrocyte_TF_columns_diff_32].max(axis=1).astype('int')
df_merged['none-none-oligodendrocyte'] = df_merged[oligodendrocyte_TF_columns_diff_32].max(axis=1).astype('int')
df_merged.drop(all_TF_columns_diff_32, axis=1, inplace=True)

for cell_EP in ['microglia', 'neuron', 'astrocyte', 'oligodendrocyte']:
    for cell_TF in ['microglia', 'neuron', 'astrocyte', 'oligodendrocyte']:
        if cell_EP == cell_TF or cell_TF == 'all':
            df_merged[f'{cell_EP}-enhancer_promoter_union_atac_500-{cell_TF}'] = df_merged.apply(
                lambda row: multiply_row(row, f'{cell_EP}_enhancer_promoter_union_atac_500',
                                         f'none-none-{cell_TF}'), axis=1)

##############################################################################################################
# Final merge and save

annotation_df = annotation_df.merge(df_merged, on=meta_columns, how='inner')
annotation_df = annotation_df.merge(baseline_df, on=meta_columns, how='inner')

annotation_df = dd.from_pandas(annotation_df, npartitions=100)
annotation_df = annotation_df.repartition(partition_size="100MB")
annotation_df = annotation_df.set_index('variant_id')
annotation_df = annotation_df.repartition(partition_size="100MB")

write_path = f'{variant_list}/annotated_variants'
os.makedirs(write_path, exist_ok=True)

write_file = f'{write_path}/annotated_variants_chr{chr_num}.parquet'
annotation_df.to_parquet(write_file, engine='pyarrow', compression='snappy', overwrite=True)

print(f"Saved annotated variants for chr{chr_num}")

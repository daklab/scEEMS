"""
Create per-gene feature matrices with cell-type-specific annotations.

For each gene, combines variant annotations with:
- ABC (Activity-By-Contact) enhancer-gene scores
- ChromBPNet cell-type-specific variant effect scores
- TF binding score summaries
- Distance to TSS

Usage:
    python create_gene_datasets.py <gene_index> <cohort> <other_flag>

Arguments:
    gene_index: 1-based index into gene list
    cohort: Cell type cohort (e.g., Mic_mega_eQTL)
    other_flag: "F" for primary genes, "T" for other-ethnicity genes

Requires config.yaml with paths configured.
"""

import os
import sys
import pandas as pd
from dask import dataframe as dd
from dask.diagnostics import ProgressBar
from tqdm import tqdm
import numpy as np
import dask
import yaml
import shutil

ProgressBar().register()
tqdm.pandas()

# Load configuration
with open('../config.yaml', 'r') as f:
    config = yaml.safe_load(f)

gene_idx = int(sys.argv[1]) - 1
cohort = sys.argv[2]
other = sys.argv[3]

data_dir = config['paths']['data_dir']
scratch_dir = config['paths'].get('scratch_dir', '/tmp')

# Set up temporary directory for dask
temp_dir = os.path.join(scratch_dir, f'{cohort}_{other}_gene_{gene_idx}')
os.makedirs(temp_dir, exist_ok=True)
dask.config.set({'temporary_directory': temp_dir})

##############################################################################################################
# TF list setup

TF_file = config['paths']['tf_file']
TF_df = pd.read_csv(TF_file, sep='\t', skiprows=1)

targets_txt = config['paths']['targets_file']
df_targets = pd.read_csv(targets_txt, sep='\t')

df_targets_subset = df_targets[df_targets['sum_stat'].str.contains("mean")]
df_targets_subset_TF = df_targets_subset[df_targets_subset['description'].str.contains("CHIP:")]
df_targets_subset_TF[['CHIP', 'TF', 'cell_type']] = df_targets_subset_TF.description.str.split(":", expand=True)

TF_microglia = [
    'ETS2', 'FLI1', 'SPI1', 'IRF8', 'PRDM1', 'CEBPA', 'CEBPB', 'CEBPD', 'CEBPE',
    'KLF11', 'KLF2', 'TFEC', 'BACH1', 'BATF', 'BATF3', 'RUNX1', 'RUNX2', 'RUNX3',
    'LYL1', 'TAL1', 'ELK3', 'MEF2C', 'CREB3L2', 'ATF4', 'MAFB', 'SALL1', 'eGFP-SALL1',
    'SMAD5', 'MEF2A', 'MEF2B', 'SMAD2', 'USF1', 'STAT3', 'eGFP-MAFG', 'NFYB', 'NRF1',
    'CREB1', 'IRF1', 'SOX9',
    'SP1', 'MAX', 'ETS1', 'HINFP', 'ZNF416', 'MEF2D', 'NFIA'
]

TF_neuronal = [
    'ASCL1', 'ASCL5', 'BHLHA15', 'BHLHE22', 'NEUROD1', 'NEUROD2', 'NEUROD6',
    'TWIST2', 'EGR4', 'SP8', 'SP9', 'MEIS3', 'PKNOX2', 'MSC', 'KLF5', 'KLF8',
    'TBR1', 'HLF',
    'MEF2B', 'SP1', 'NRF1', 'GATA3', 'MEF2C', 'EGR2', 'GFI1B', 'TEAD1'
]

TF_astrocyte = [
    'MEIS2', 'NFIB', 'TGIF1', 'EMX2', 'LHX2', 'RFX2', 'RFX4', 'RORA', 'RORB',
    'SOX1', 'SOX2', 'SOX21', 'SOX9', 'NFATC4', 'SOX5', 'POU3F2', 'POU3F3',
    'POU3F4', 'FOXG1', 'FOXO1', 'SP5',
    'RFX1', 'MYB', 'NFIC', 'FOXK2', 'ATF3', 'SP1', 'ATF1', 'ZNF416', 'MEF2D'
]

TF_oligodendrocyte = [
    'SOX10', 'SOX13', 'SOX21', 'SOX3', 'SOX6', 'SOX8', 'NHLH2', 'NFIX',
    'SP7', 'NFE2', 'E2F1', 'CREB5', 'POU3F3', 'MYCN',
    'CTCF', 'SP1', 'RFX1', 'FOSL2', 'MYB', 'E2F7', 'ATF7', 'RUNX1', 'USF1',
    'YY2', 'GFI1B', 'FOXJ2', 'ETV1'
]

microglia_TF_columns = df_targets_subset_TF[df_targets_subset_TF['TF'].isin(TF_microglia)].identifier.tolist()
neuronal_TF_columns = df_targets_subset_TF[df_targets_subset_TF['TF'].isin(TF_neuronal)].identifier.tolist()
astrocyte_TF_columns = df_targets_subset_TF[df_targets_subset_TF['TF'].isin(TF_astrocyte)].identifier.tolist()
oligodendrocyte_TF_columns = df_targets_subset_TF[df_targets_subset_TF['TF'].isin(TF_oligodendrocyte)].identifier.tolist()

TF_all_list = list(set(
    TF_df['Target of assay'].tolist() + TF_microglia + TF_neuronal + TF_astrocyte + TF_oligodendrocyte
))
all_TF_columns = df_targets_subset_TF[df_targets_subset_TF['TF'].isin(TF_all_list)].identifier.tolist()

microglia_TF_columns_diff_32 = [f'diff_32_{i}' for i in microglia_TF_columns]
neuronal_TF_columns_diff_32 = [f'diff_32_{i}' for i in neuronal_TF_columns]
astrocyte_TF_columns_diff_32 = [f'diff_32_{i}' for i in astrocyte_TF_columns]
oligodendrocyte_TF_columns_diff_32 = [f'diff_32_{i}' for i in oligodendrocyte_TF_columns]
all_TF_columns_diff_32 = list(set([f'diff_32_{i}' for i in all_TF_columns]))

##############################################################################################################
# Load gene information

susie_dir = os.path.join(data_dir, 'susie_vars_pips')
variant_list = os.path.join(susie_dir, 'variant_list')

if other == "F":
    gene_list_path = os.path.join(data_dir, f'training_data/{cohort}/list_genes.csv')
elif other == "T":
    gene_list_path = os.path.join(data_dir, f'training_data/{cohort}/list_genes_other.csv')

genes_df = pd.read_csv(gene_list_path, sep="\t")
gene = genes_df.iloc[gene_idx]["gene_id"]
chr_val = genes_df.iloc[gene_idx]["chr"]

abc_data_dir = config['paths']['abc_data_dir']
gene_id_gene_name_df = pd.read_csv(f'{abc_data_dir}/ABC_gene_id_name_mapping.csv', sep=',')
gene_id_gene_name_df = gene_id_gene_name_df[gene_id_gene_name_df['gene_id'] == gene]
gene_list_gene_id = gene_id_gene_name_df['gene_name'].tolist()

##############################################################################################################
# Load PIP data

if other == "F":
    all_parquet_dir = f'{susie_dir}/{cohort}/PIP_all_parquet'
    all_parquet_df = dd.read_parquet(f'{all_parquet_dir}/PIP_all.parquet/chr={chr_val}', engine='pyarrow')
elif other == "T":
    all_parquet_dir = f'{susie_dir}/{cohort}/PIP_all_other_parquet'
    all_parquet_df = dd.read_parquet(f'{all_parquet_dir}/PIP_all_other.parquet/chr={chr_val}', engine='pyarrow')

all_parquet_df = all_parquet_df[['variant_id', 'chr', 'pos', 'ref', 'alt', 'gene_id', 'pip']].reset_index(drop=True)
all_parquet_df = all_parquet_df.set_index('variant_id')
all_parquet_df = all_parquet_df[all_parquet_df['gene_id'] == gene]
all_parquet_df = all_parquet_df.persist()

# Load variant annotations
annotation_path_top_pip = f'{variant_list}/annotated_variants_just_variants'
annotation_file_top_pip = f'{annotation_path_top_pip}/annotated_just_variants_{chr_val}.parquet'
annotation_df_top_pip = dd.read_parquet(annotation_file_top_pip, engine='pyarrow')

all_pip_df_annotated = all_parquet_df.merge(annotation_df_top_pip,
                                            left_index=True, right_index=True, how='inner')
all_pip_df_annotated = all_pip_df_annotated.compute()
all_pip_df_annotated = all_pip_df_annotated.reset_index()
all_pip_df_annotated_variants = all_pip_df_annotated[['variant_id', 'chr', 'pos', 'ref', 'alt', 'gene_id', 'pip']]

pos_list_filter = all_pip_df_annotated_variants['pos'].tolist()
annotation_full_file = f'{variant_list}/annotated_variants/annotated_variants_{chr_val}.parquet'
annotation_df_full = dd.read_parquet(annotation_full_file, engine='pyarrow',
                                     filters=[('BP', 'in', pos_list_filter)])
annotation_df_full = annotation_df_full.persist()

##############################################################################################################
# Process gene data

training_data_df = all_pip_df_annotated[all_pip_df_annotated["gene_id"] == gene].copy()
training_data_df = training_data_df.merge(gene_id_gene_name_df, on='gene_id', how='left')

mean_TSS = gene_id_gene_name_df['gene_TSS'].mean()
training_data_df['distance_TSS'] = mean_TSS - training_data_df['BP']
training_data_df['abs_distance_TSS'] = training_data_df['distance_TSS'].abs()
training_data_df['abs_distance_TSS_log'] = np.log(training_data_df['abs_distance_TSS'])

training_data_df = training_data_df.sort_values(by='pip', ascending=False).groupby('variant_id').first().reset_index()
training_data_df = training_data_df.drop(columns=['index', 'gene_id', 'gene_name'], errors='ignore')

##############################################################################################################
# ABC scores

abc_names = ['microglia', 'astrocyte', 'oligodendrocyte', 'neuron']
abc_files = [f'{abc_data_dir}/ABC_results_{f}_v2/{f}/Predictions/EnhancerPredictionsAllPutative.tsv.gz'
             for f in abc_names]
abc_dfs = [pd.read_csv(f, sep='\t') for f in abc_files]

chrombpnet_dir = config['paths']['chrombpnet_dir']
chrombpnet_dfs_dict = {
    f: pd.read_csv(f'{chrombpnet_dir}/chrombpnet_{f}_{chr_val}_variant_peak_pairs_scored.csv', sep='\t')
    for f in abc_names
}

# Filter ABC for target genes
for i in range(len(abc_dfs)):
    abc_dfs[i] = abc_dfs[i][abc_dfs[i]['TargetGene'].isin(gene_list_gene_id)]


def get_max_abc_score_context_length(chr_val, bp, abc_df, context_length=1024):
    """Get max ABC score within a context window around a position."""
    start_pos = bp - context_length
    end_pos = bp + context_length
    mask = (abc_df['chr'] == chr_val) & (
        ((abc_df['start'] >= start_pos) & (abc_df['start'] <= end_pos)) |
        ((abc_df['end'] >= start_pos) & (abc_df['end'] <= end_pos)) |
        ((abc_df['start'] <= start_pos) & (abc_df['end'] >= end_pos))
    )
    if mask.any():
        return abc_df.loc[mask, 'ABC.Score'].max()
    return 0.0


def process_abc_scores(training_data_df, abc_dfs, abc_names):
    """Add ABC scores for each cell type."""
    for abc_column in abc_names:
        training_data_df[f'ABC_score_{abc_column}'] = 0.0
    for abc_df, abc_name in zip(abc_dfs, abc_names):
        for idx, row in tqdm(training_data_df.iterrows()):
            max_abc_score = get_max_abc_score_context_length(row['CHR'], row['BP'], abc_df)
            training_data_df.loc[idx, f'ABC_score_{abc_name}'] = max_abc_score
    training_data_df.drop(columns=['gene_TSS'], inplace=True, errors='ignore')
    return training_data_df


def pivot_chrombpnet_scores(chrombpnet_df, cell_name):
    """Process and pivot chrombpnet scores to wide format."""
    unique_variants = chrombpnet_df['variant_id'].unique()
    final_df = pd.DataFrame({'variant_id': unique_variants})
    metrics = [col for col in ['log_counts_diff_chrombpnet', 'log_probs_diff_abs_sum_chrombpnet',
                               'probs_jsd_diff_chrombpnet']
               if col in chrombpnet_df.columns]
    for metric in metrics:
        assays = chrombpnet_df['assay'].unique()
        for assay in assays:
            mask = (chrombpnet_df['assay'] == assay) & chrombpnet_df[metric].notna()
            assay_data = chrombpnet_df[mask][['variant_id', metric]]
            if len(assay_data) == 0:
                continue
            max_values, min_values, abs_max_values = {}, {}, {}
            for _, row in assay_data.iterrows():
                var_id = row['variant_id']
                value = row[metric]
                if pd.isna(value):
                    continue
                if var_id not in max_values or (value > 0 and value > max_values[var_id]):
                    max_values[var_id] = max(0, value)
                if var_id not in min_values or (value < 0 and value < min_values[var_id]):
                    min_values[var_id] = min(0, value)
                abs_value = abs(value)
                if var_id not in abs_max_values or abs_value > abs_max_values[var_id]:
                    abs_max_values[var_id] = abs_value
            final_df[f'{cell_name}-{metric}_max_{assay}'] = final_df['variant_id'].map(
                lambda x: max_values.get(x, 0))
            final_df[f'{cell_name}-{metric}_min_{assay}'] = final_df['variant_id'].map(
                lambda x: min_values.get(x, 0))
            final_df[f'{cell_name}-{metric}_abs_max_{assay}'] = final_df['variant_id'].map(
                lambda x: abs_max_values.get(x, 0))
    final_df = final_df.fillna(0)
    return final_df


def process_chrombpnet_scores(training_data_df, chrombpnet_dfs_dict, abc_names):
    """Process chrombpnet scores for each cell type and merge with training data."""
    merged_training_data_df = training_data_df.copy(deep=True)
    for cell_name in abc_names:
        print(f"Processing chrombpnet scores for {cell_name}...")
        if cell_name not in chrombpnet_dfs_dict:
            continue
        chrombpnet_df = chrombpnet_dfs_dict[cell_name]
        available_metrics = [col for col in ['log_counts_diff_chrombpnet',
                                             'log_probs_diff_abs_sum_chrombpnet',
                                             'probs_jsd_diff_chrombpnet']
                             if col in chrombpnet_df.columns]
        if not available_metrics:
            print(f"Warning: No required metrics found for {cell_name}")
            continue
        chrombpnet_df_pivot = pivot_chrombpnet_scores(chrombpnet_df, cell_name=cell_name)
        merged_training_data_df = merged_training_data_df.merge(
            chrombpnet_df_pivot, on='variant_id', how='left')
        new_cols = [col for col in merged_training_data_df.columns
                    if col.startswith(f'{cell_name}-') and col not in training_data_df.columns]
        if new_cols:
            merged_training_data_df[new_cols] = merged_training_data_df[new_cols].fillna(0)
    return merged_training_data_df


training_data_df = process_abc_scores(training_data_df, abc_dfs, abc_names)
training_data_df = process_chrombpnet_scores(training_data_df, chrombpnet_dfs_dict, abc_names)

##############################################################################################################
# TF score summaries

training_data_df = training_data_df.drop_duplicates(subset=['variant_id'])
training_data_df = dd.from_pandas(training_data_df, npartitions=10)
training_data_df = training_data_df.persist()
training_data_df = training_data_df.set_index('variant_id')
training_data_df = training_data_df.drop(columns=['CHR', 'BP', 'REF', 'ALT', 'SNP'])
training_data_df = training_data_df.merge(annotation_df_full,
                                          left_index=True, right_index=True, how='left')

# Calculate TF score summaries per cell type
for tf_name, tf_cols in [('microglia', microglia_TF_columns_diff_32),
                          ('astrocyte', astrocyte_TF_columns_diff_32),
                          ('oligodendrocyte', oligodendrocyte_TF_columns_diff_32),
                          ('neuron', neuronal_TF_columns_diff_32)]:
    training_data_df[f'max_{tf_name}_TF_score'] = training_data_df[tf_cols].max(axis=1)
    training_data_df[f'min_{tf_name}_TF_score'] = training_data_df[tf_cols].min(axis=1)
    training_data_df[f'abs_max_{tf_name}_TF_score'] = training_data_df[tf_cols].apply(
        lambda x: max(x.abs().max(), 0), axis=1)

training_data_df = training_data_df.reset_index()
training_data_df = training_data_df.repartition(partition_size="100MB")
training_data_df = training_data_df.persist()

##############################################################################################################
# Save output

output_dir = os.path.join(data_dir, f'training_data/{cohort}/all_variants')
os.makedirs(output_dir, exist_ok=True)

write_dir = os.path.join(output_dir, gene)
os.makedirs(write_dir, exist_ok=True)

training_data_df.to_parquet(f'{write_dir}/annotated_data_{cohort}_{chr_val}.parquet',
                            engine='pyarrow', compression='snappy', overwrite=True)

print(f"Finished writing {gene} data for {cohort} and {chr_val}")

# Clean up temp directory
shutil.rmtree(temp_dir, ignore_errors=True)

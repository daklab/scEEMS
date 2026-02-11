import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import sys
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
import dask
import yaml

import pyranges as pr


def analyze_spaces_in_line(file_path, line_number):
    try:
        # Read the specified line from the file
        with open(file_path, 'r') as file:
            for i, line in enumerate(file, 1):
                if i == line_number:
                    target_line = line.rstrip('\n')
                    break
            else:
                return f"Line {line_number} not found in the file."
        # Split the line by any whitespace to get all tokens
        tokens = target_line.split()
        # Initialize list to store space counts
        space_counts = []
        # Find the original positions of each token in the line
        positions = []
        start_pos = 0
        for token in tokens:
            pos = target_line.find(token, start_pos)
            positions.append(pos)
            start_pos = pos + len(token)
        # Calculate spaces between consecutive tokens
        for i in range(len(positions) - 1):
            current_token = tokens[i]
            current_pos = positions[i]
            next_pos = positions[i + 1]
            # Calculate the number of spaces
            spaces = next_pos - (current_pos + len(current_token))
            space_counts.append((tokens[i], tokens[i + 1], spaces))
        # Print results
        print(f"Line {line_number}: {target_line}")
        print("\nSpaces between consecutive strings:")
        for prev, next, count in space_counts:
            print(f"'{prev}' to '{next}': {count} spaces")
        return space_counts
    except FileNotFoundError:
        return f"File '{file_path}' not found."
    except Exception as e:
        return f"An error occurred: {str(e)}"


#abc_name = 'microglia'

# Set up progress bar
tqdm.pandas()
pbar = ProgressBar(dt=1)
pbar.register()

# Configure dask temporary directory
import datetime

def overlap_placseq_simple_fix(df, placseq_df):
    """
    Simple, robust version that avoids indexing issues.
    Works directly with the necessary data without complex merges.
    """
    # Extract what we need as arrays - no copying of large dataframe
    variant_ids = df['variant_id'].values
    gene_ids = df['gene_id'].values
    n_rows = len(df)
    # Parse positions efficiently
    positions = np.array([int(vid.split(':')[1]) for vid in variant_ids])
    # Initialize result arrays
    microglia_counts = np.zeros(n_rows, dtype=np.int32)
    oligodendrocyte_counts = np.zeros(n_rows, dtype=np.int32)
    neuron_counts = np.zeros(n_rows, dtype=np.int32)
    # Cell type mapping
    cell_type_dict = {
        'microglia': ('PU1_enhancer_interactions', microglia_counts),
        'oligodendrocyte': ('Olig2_enhancer_interactions', oligodendrocyte_counts),
        'neuron': ('NeuN_enhancer_interactions', neuron_counts)
    }
    # Create placseq lookup by gene and interaction type
    placseq_lookup = {}
    for interaction_type in ['PU1_enhancer_interactions', 'Olig2_enhancer_interactions', 'NeuN_enhancer_interactions']:
        type_data = placseq_df[placseq_df['interaction_type'] == interaction_type]
        if not type_data.empty:
            # Group by gene_id and store intervals
            grouped = type_data.groupby('gene_id')[['start', 'end']].apply(
                lambda x: (x['start'].values, x['end'].values)
            ).to_dict()
            placseq_lookup[interaction_type] = grouped
    # Process each variant
    for i in tqdm(range(n_rows), desc="Processing variants"):
        gene_id = gene_ids[i]
        pos = positions[i]
        # Check each cell type
        for cell_type, (interaction_type, count_array) in cell_type_dict.items():
            if interaction_type in placseq_lookup and gene_id in placseq_lookup[interaction_type]:
                starts, ends = placseq_lookup[interaction_type][gene_id]
                # Count overlaps
                overlaps = (starts <= pos) & (ends >= pos)
                count_array[i] = np.sum(overlaps)
    # Add results to original dataframe
    df['microglia_placseq_count'] = microglia_counts
    df['oligodendrocyte_placseq_count'] = oligodendrocyte_counts
    df['neuron_placseq_count'] = neuron_counts
    return df


def get_max_abc_score(chr, bp, abc_df):
    # Find overlapping regions
    mask = (abc_df['chr'] == chr) & \
           (abc_df['start'] <= bp) & \
           (abc_df['end'] >= bp
            )
    if mask.any():
        return abc_df.loc[mask, 'ABC.Score'].max()
    return 0.0

def process_abc_scores(training_data_df, abc_dfs, abc_names):
    # Initialize ABC score columns with zeros
    for abc_column in abc_names:
        training_data_df[f'ABC_score_{abc_column}'] = 0.0
        #training_data_df[f'ABC_context_length_{abc_column}'] = 0.0
    # Process each ABC dataset
    for abc_df, abc_name in zip(abc_dfs, abc_names):
        # For each variant in training data
        for idx, row in tqdm(training_data_df.iterrows()):
            gene_id_val = row['gene_id']
            abc_df_subset = abc_df[abc_df['TargetGeneEnsembl_ID'] == gene_id_val]
            max_abc_score = get_max_abc_score(row['CHR'], row['BP'], abc_df_subset)
            #max_abc_context_length = get_max_abc_score_context_length(row['CHR'], row['BP'], abc_df)
            training_data_df.loc[idx, f'ABC_score_{abc_name}'] = max_abc_score
            #training_data_df.loc[idx, f'ABC_context_length_{abc_name}'] = max_abc_context_length
    # Drop unnecessary columns
    #training_data_df.drop(columns=['gene_TSS'], inplace=True)
    return training_data_df


import pandas as pd
import numpy as np


def sort_and_filter_by_cell_type(pred_merged, cell_type, cell_type_abc_plac_dict):
    """
    Sort and filter dataframe based on cell type configuration.

    Parameters:
    - pred_merged: DataFrame to process
    - cell_type: Cell type key (e.g., 'Mic_mega_eQTL')
    - cell_type_abc_plac_dict: Dictionary with cell type configurations

    Returns:
    - Filtered and sorted DataFrame
    """
    # Get configuration for the specific cell type
    config = cell_type_abc_plac_dict[cell_type]
    placseq_col = config['placseq']
    abc_col = config['abc']
    # Create a copy to avoid modifying original
    df = pred_merged.copy()
    # Handle missing placseq column case
    if placseq_col is None or placseq_col not in df.columns:
        # Only use ABC score for sorting
        df['sort_priority'] = df[abc_col]
        df['has_placseq'] = False
    else:
        # Use placseq as primary sort, ABC as secondary
        df['sort_priority'] = df[placseq_col]
        df['has_placseq'] = True
        # For ties in placseq (including zeros), use ABC score as tiebreaker
        df['secondary_sort'] = df[abc_col]
    # Remove rows where both values are zero (or only ABC is zero if no placseq)
    if placseq_col is None or placseq_col not in df.columns:
        # Only ABC column exists
        df_filtered = df[df[abc_col] != 0].copy()
    else:
        # Both columns exist - remove rows where both are zero
        df_filtered = df[~((df[placseq_col] == 0) & (df[abc_col] == 0))].copy()
    # Sort by variant_id and then by priority scores
    if placseq_col is None or placseq_col not in df_filtered.columns:
        # Sort by variant_id, then by ABC score (descending)
        df_sorted = df_filtered.sort_values(
            ['variant_id', 'sort_priority'],
            ascending=[True, False]
        )
    else:
        # Sort by variant_id, then by placseq (descending), then by ABC (descending)
        df_sorted = df_filtered.sort_values(
            ['variant_id', 'sort_priority', 'secondary_sort'],
            ascending=[True, False, False]
        )
    # Take the top row for each variant_id (highest priority)
    result = df_sorted.groupby('variant_id').first().reset_index()
    # Clean up helper columns
    columns_to_drop = ['sort_priority', 'has_placseq']
    if 'secondary_sort' in result.columns:
        columns_to_drop.append('secondary_sort')
    result = result.drop(columns=columns_to_drop)
    return result





# Load config
with open(os.path.join(os.path.dirname(__file__), '..', 'config.yaml'), 'r') as f:
    config = yaml.safe_load(f)

dask.config.set({'temporary_directory': config['paths']['scratch_dir']})

chr = sys.argv[1]

cell_type = sys.argv[2]



cell_type_abc_plac_dict = {
   'Mic_mega_eQTL': {'placseq': 'microglia_placseq_count', 'abc': 'ABC_score_microglia', 'abc_name': 'microglia'},
   'Ast_mega_eQTL': {'placseq': None, 'abc': 'ABC_score_astrocyte', 'abc_name': 'astrocyte'},
    'Exc_mega_eQTL': {'placseq': 'neuron_placseq_count', 'abc': 'ABC_score_neuron', 'abc_name': 'neuron'},
    'Inh_mega_eQTL': {'placseq': 'neuron_placseq_count', 'abc': 'ABC_score_neuron', 'abc_name': 'neuron'},
    'Oli_mega_eQTL': {'placseq': 'oligodendrocyte_placseq_count', 'abc': 'ABC_score_oligodendrocyte', 'abc_name': 'oligodendrocyte'},
    'OPC_mega_eQTL': {'placseq': 'oligodendrocyte_placseq_count', 'abc': 'ABC_score_oligodendrocyte', 'abc_name': 'oligodendrocyte'},
}


cell_type_abc_plac_dict_cell_type = cell_type_abc_plac_dict[cell_type]

baseline_dir = config['paths']['baseline_annot_dir']

baseline_file = f'{baseline_dir}/baseline_chr{chr}.annot.gz'

baseline_df = pd.read_csv(baseline_file, sep='\t',  usecols = ['CHR', 'SNP', 'BP', 'A1', 'A2'])

annotation_df = baseline_df.copy(deep=True)

#rename A2 to ref and A1 to alt



annotation_df = annotation_df.rename(columns={'A1': 'alt', 'A2': 'ref'})

#append chr to CHR

#pred_prob_threshold = 0.90
pred_prob_threshold_dict = {'Mic_mega_eQTL': 0.95,
                            'Ast_mega_eQTL': 0.97,
                            'Exc_mega_eQTL': 0.94,
                            'Inh_mega_eQTL': 0.94,
                            'OPC_mega_eQTL': 0.91,
                            'Oli_mega_eQTL': 0.98
                            }

pip_threshold = 0.10

gene_id_gene_name_dir = config['paths']['abc_data_dir']
gene_id_gene_name_df = pd.read_csv(f'{gene_id_gene_name_dir}/ABC_gene_id_name_mapping.csv', sep=',')



write_MAGMA_dir = os.path.join(config['paths']['output_dir'], 'aggregate_results', 'MAGMA')





# Create the output directory if it doesn't exist
if not os.path.exists(write_MAGMA_dir):
    os.makedirs(write_MAGMA_dir)

cohorts = ['Ast_mega_eQTL',
              'Exc_mega_eQTL',
              'Inh_mega_eQTL',
              'Mic_mega_eQTL' ,
              'Oli_mega_eQTL',
              'OPC_mega_eQTL']

# for cell_type in tqdm(cohorts):
# Load the predictions dataframe with high probability variants


write_MAGMA_dir_cell_type = f'{write_MAGMA_dir}/{cell_type}'

if not os.path.exists(write_MAGMA_dir_cell_type):
    os.makedirs(write_MAGMA_dir_cell_type)

# Create prediction and pip subdirectories
write_MAGMA_dir_pred = f'{write_MAGMA_dir_cell_type}/prediction'
write_MAGMA_dir_pip = f'{write_MAGMA_dir_cell_type}/pip'

write_MAGMA_dir_variant_gene = f'{write_MAGMA_dir_cell_type}/variant_gene_lists'


if not os.path.exists(write_MAGMA_dir_pred):
    os.makedirs(write_MAGMA_dir_pred)

if not os.path.exists(write_MAGMA_dir_pip):
    os.makedirs(write_MAGMA_dir_pip)


if not os.path.exists(write_MAGMA_dir_variant_gene):
    os.makedirs(write_MAGMA_dir_variant_gene)


# Set up the MAGMA output file paths for both prediction and pip
MAGMA_file_pred = f'{write_MAGMA_dir_pred}/{cell_type}_chr{chr}_MAGMA.genes.annot'
MAGMA_file_pip = f'{write_MAGMA_dir_pip}/{cell_type}_chr{chr}_MAGMA.genes.annot'


#delete if os.path.exists(MAGMA_file_pred):
if os.path.exists(MAGMA_file_pred):
    os.remove(MAGMA_file_pred)


if os.path.exists(MAGMA_file_pip):
    os.remove(MAGMA_file_pip)


predictions_dir = os.path.join(config['paths']['data_dir'], cell_type)
predictions_df = dd.read_parquet(f'{predictions_dir}/predictions_parquet_catboost/predictions.parquet')
predictions_df = predictions_df[predictions_df['chr'] == f'chr{chr}']
predictions_df = predictions_df.persist()
predictions_df_top_pred = predictions_df[predictions_df['pred_prob'] > pred_prob_threshold_dict[f'{cell_type}']].compute()
predictions_df_top_pip = predictions_df[predictions_df['pip'] > pip_threshold].compute()
# rename pos to BP chr to CHR
predictions_df_top_pred = predictions_df_top_pred.rename(columns={'pos': 'BP', 'chr': 'CHR'})
predictions_df_top_pred['CHR'] = predictions_df_top_pred['CHR'].str.replace('chr', '')
predictions_df_top_pred['CHR'] = predictions_df_top_pred['CHR'].astype(int)
predictions_df_top_pip = predictions_df_top_pip.rename(columns={'pos': 'BP', 'chr': 'CHR'})
predictions_df_top_pip['CHR'] = predictions_df_top_pip['CHR'].str.replace('chr', '')
predictions_df_top_pip['CHR'] = predictions_df_top_pip['CHR'].astype(int)
pred_merged = pd.merge(annotation_df, predictions_df_top_pred, on=['CHR', 'BP', 'ref', 'alt'], how='inner')
pip_merged = pd.merge(annotation_df, predictions_df_top_pip, on=['CHR', 'BP', 'ref', 'alt'], how='inner')



AFR_AMR_variants = os.path.join(config['paths']['multi_ancestry_predictions_dir'], cell_type)

pred_merged_AFR_AMR = pd.read_csv(f'{AFR_AMR_variants}/MAGMA_predictions_{cell_type}_chr{chr}.tsv.gz', sep='\t')

pred_merged_AFR_AMR['CHR'] = pred_merged_AFR_AMR['CHR'].str.replace('chr', '')

pred_merged_AFR_AMR['variant_id'] = pred_merged_AFR_AMR['CHR'].astype(str) + ':' + pred_merged_AFR_AMR['BP'].astype(str) + ':' + pred_merged_AFR_AMR['REF'] + ':' + pred_merged_AFR_AMR['ALT']
#rename REF and ALT to ref and alt

pred_merged_AFR_AMR = pred_merged_AFR_AMR.rename(columns={'REF': 'ref', 'ALT': 'alt'})


pred_merged_AFR_AMR = pred_merged_AFR_AMR[pred_merged_AFR_AMR['pred_prob'] > pred_prob_threshold_dict[f'{cell_type}']]

pred_merged = pd.concat([pred_merged, pred_merged_AFR_AMR], ignore_index=True)



# get first row grouped by CHR          SNP        BP ref alt and gene_id
pred_merged = pred_merged.groupby(['CHR', 'SNP', 'BP', 'ref', 'alt', 'gene_id']).first().reset_index()





abc_names = ['microglia', 'astrocyte', 'oligodendrocyte', 'neuron']

abc_dir = config['paths']['abc_data_dir']

abc_files = [f'{abc_dir}/ABC_results_{f}_v2/{f}/Predictions/EnhancerPredictionsAllPutative.tsv.gz' for f in abc_names]

abc_dfs = [pd.read_csv(f, sep='\t') for f in abc_files]


if cell_type == 'Mic_mega_eQTL':
    abc_df_single = abc_dfs[0]
elif cell_type == 'Ast_mega_eQTL':
    abc_df_single = abc_dfs[1]
elif cell_type in ['Exc_mega_eQTL', 'Inh_mega_eQTL']:
    abc_df_single = abc_dfs[3]
elif cell_type in ['Oli_mega_eQTL', 'OPC_mega_eQTL']:
    abc_df_single = abc_dfs[2]


abc_df_single = abc_df_single[abc_df_single['ABC.Score'] >= 0.005]

abc_df_single = abc_df_single[abc_df_single['chr'].str.replace('chr', '') == str(chr)]

#replace chr in chr column and make it int
abc_df_single['chr'] = abc_df_single['chr'].str.replace('chr', '')
abc_df_single['chr'] = abc_df_single['chr'].astype(int)


placseq_file = config['paths']['placseq_file']

placseq_df = pd.read_csv(placseq_file, sep='\t')



pred_merged = process_abc_scores(pred_merged, [abc_df_single], [cell_type_abc_plac_dict_cell_type['abc_name']])

pred_merged = overlap_placseq_simple_fix(pred_merged, placseq_df)


pred_merged = sort_and_filter_by_cell_type(pred_merged, cell_type, cell_type_abc_plac_dict)


#filter by BP == 125149195
#pred_merged[pred_merged['BP'] == 86156833]
#pred_merged[pred_merged['BP'] == 207693]
#pred_merged[pred_merged['BP'] == 130180927]
#pred_merged[pred_merged['BP'] == 134254132]


#by

#find rows where BP is duplicated
#pred_merged[pred_merged.duplicated(subset=['BP'], keep=False)].sort_values(by='BP')

# number of rows grouped by gene_id in pred_merged
pred_merged.groupby('gene_id').size().reset_index(name='counts').sort_values(by='counts', ascending=False)
pip_merged.groupby('gene_id').size().reset_index(name='counts').sort_values(by='counts', ascending=False)
pred_merged_genes = pred_merged['gene_id'].unique()
pip_merged_genes = pip_merged['gene_id'].unique()


pred_merged.to_csv(f'{write_MAGMA_dir_variant_gene}/MAGMA_variant_gene_list_{cell_type}_chr{chr}_pred.tsv', sep='\t', index=False)
pip_merged.to_csv(f'{write_MAGMA_dir_variant_gene}/MAGMA_variant_gene_list_{cell_type}_chr{chr}_pip.tsv', sep='\t', index=False)

#pred_merged[pred_merged['gene_id'] == 'ENSG00000073921']
#pred_merged[pred_merged['gene_id'] == 'ENSG00000074266']


# Create and write to the pred MAGMA file
first_line = '# window_up = 0'
second_line = '# window_down = 0'

# Write the prediction MAGMA file
with open(MAGMA_file_pred, 'w') as f:
    # Write the header lines
    f.write(f"{first_line}\n")
    f.write(f"{second_line}\n")
    # Iterate through genes and write each line
    for gene in pred_merged_genes:
        pred_merged_gene = pred_merged[pred_merged['gene_id'] == gene]
        start_window = pred_merged_gene['BP'].min() - 5
        end_window = pred_merged_gene['BP'].max() + 5
        window_string = f'{chr}:{start_window}:{end_window}'
        rsid_list = pred_merged_gene['SNP'].tolist()
        # Create rsid string with 4 spaces between each rsid
        rsid_string = '    '.join(rsid_list)
        # Create the line with proper spacing:
        # 1 space after gene_id, 4 spaces after window_string
        write_line = f'{gene} {window_string}    {rsid_string}'
        f.write(f"{write_line}\n")

print(f"Prediction MAGMA file created at: {MAGMA_file_pred}")
print(f"Total prediction genes written: {len(pred_merged_genes)}")

# Write the pip MAGMA file
with open(MAGMA_file_pip, 'w') as f:
    # Write the header lines
    f.write(f"{first_line}\n")
    f.write(f"{second_line}\n")
    # Iterate through genes and write each line
    for gene in pip_merged_genes:
        pip_merged_gene = pip_merged[pip_merged['gene_id'] == gene]
        start_window = pip_merged_gene['BP'].min() - 5
        end_window = pip_merged_gene['BP'].max() + 5
        window_string = f'{chr}:{start_window}:{end_window}'
        rsid_list = pip_merged_gene['SNP'].tolist()
        # Create rsid string with 4 spaces between each rsid
        rsid_string = '    '.join(rsid_list)
        # Create the line with proper spacing:
        # 1 space after gene_id, 4 spaces after window_string
        write_line = f'{gene} {window_string}    {rsid_string}'
        f.write(f"{write_line}\n")

print(f"PIP MAGMA file created at: {MAGMA_file_pip}")
print(f"Total PIP genes written: {len(pip_merged_genes)}")





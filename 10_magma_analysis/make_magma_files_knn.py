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


# Load config
with open(os.path.join(os.path.dirname(__file__), '..', 'config.yaml'), 'r') as f:
    config = yaml.safe_load(f)

dask.config.set({'temporary_directory': config['paths']['scratch_dir']})

chr = sys.argv[1]

cell_type = sys.argv[2]


baseline_dir = config['paths']['baseline_annot_dir']

baseline_file = f'{baseline_dir}/baseline_chr{chr}.annot.gz'

baseline_df = pd.read_csv(baseline_file, sep='\t',  usecols = ['CHR', 'SNP', 'BP', 'A1', 'A2'])

annotation_df = baseline_df.copy(deep=True)

#rename A2 to ref and A1 to alt

annotation_df = annotation_df.rename(columns={'A1': 'alt', 'A2': 'ref'})

#append chr to CHR

pred_prob_threshold = 0.95
pip_threshold = 0.10

gene_id_gene_name_dir = config['paths']['abc_data_dir']
gene_id_gene_name_df = pd.read_csv(f'{gene_id_gene_name_dir}/ABC_gene_id_name_mapping.csv', sep=',')


gene_id_gene_name_df_list = gene_id_gene_name_df['gene_id'].unique().tolist()

#extract row with gene_id == ENSG00000112096
#gene_id = 'ENSG00000112096'
#gene_id_gene_name_df[gene_id_gene_name_df['gene_id'] == gene_id]


write_MAGMA_dir = os.path.join(config['paths']['output_dir'], 'aggregate_results', 'MAGMA_knn')

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

if not os.path.exists(write_MAGMA_dir_pred):
    os.makedirs(write_MAGMA_dir_pred)

if not os.path.exists(write_MAGMA_dir_pip):
    os.makedirs(write_MAGMA_dir_pip)

# Set up the MAGMA output file paths for both prediction and pip
MAGMA_file_pred = f'{write_MAGMA_dir_pred}/{cell_type}_chr{chr}_MAGMA.genes.annot'
MAGMA_file_pip = f'{write_MAGMA_dir_pip}/{cell_type}_chr{chr}_MAGMA.genes.annot'

#delete if os.path.exists(MAGMA_file_pred):
if os.path.exists(MAGMA_file_pred):
    os.remove(MAGMA_file_pred)


if os.path.exists(MAGMA_file_pip):
    os.remove(MAGMA_file_pip)





MAGMA_original_dir = os.path.join(config['paths']['output_dir'], 'aggregate_results', 'MAGMA', cell_type, 'variant_gene_lists')

pred_merged = pd.read_csv(f'{MAGMA_original_dir}/MAGMA_variant_gene_list_{cell_type}_chr{chr}_pred.tsv', sep='\t')

pip_merged = pd.read_csv(f'{MAGMA_original_dir}/MAGMA_variant_gene_list_{cell_type}_chr{chr}_pip.tsv', sep='\t')


# number of rows grouped by gene_id in pred_merged
pred_merged.groupby('gene_id').size().reset_index(name='counts').sort_values(by='counts', ascending=False)
pip_merged.groupby('gene_id').size().reset_index(name='counts').sort_values(by='counts', ascending=False)

num_pred_gene_variants = pred_merged.groupby('gene_id').size().reset_index(name='counts').sort_values(by='counts', ascending=False)
num_pip_gene_variants = pip_merged.groupby('gene_id').size().reset_index(name='counts').sort_values(by='counts', ascending=False)


pred_merged_genes = pred_merged['gene_id'].unique()
pip_merged_genes = pip_merged['gene_id'].unique()


#keep only genes that are in the gene_id_gene_name_df_list
pred_merged_genes = [gene for gene in pred_merged_genes if gene in gene_id_gene_name_df_list]
pip_merged_genes = [gene for gene in pip_merged_genes if gene in gene_id_gene_name_df_list]


# Create and write to the pred MAGMA file
first_line = '# window_up = 0'
second_line = '# window_down = 0'

# Write the prediction MAGMA file
with open(MAGMA_file_pred, 'w') as f:
    # Write the header lines
    f.write(f"{first_line}\n")
    f.write(f"{second_line}\n")
    # Iterate through genes and write each line
    for gene in tqdm(pred_merged_genes):
        print(gene)
        num_variants = num_pred_gene_variants[num_pred_gene_variants['gene_id'] == gene]['counts'].values[0]
        TSS_position = gene_id_gene_name_df[gene_id_gene_name_df['gene_id'] == gene]['gene_TSS'].values[0]
        #get the n closest variants to TSS_position in annotation_df using the column BP
        variant_subset = annotation_df.copy(deep=True)
        variant_subset['distance_to_TSS'] = np.abs(variant_subset['BP'] - TSS_position)
        closest_variants = variant_subset.nsmallest(num_variants, 'distance_to_TSS')
        start_window = closest_variants['BP'].min() - 5
        end_window = closest_variants['BP'].max() + 5
        window_string = f'{chr}:{start_window}:{end_window}'
        print(window_string)
        rsid_list = closest_variants['SNP'].tolist()
        # Create rsid string with 4 spaces between each rsid
        rsid_string = '    '.join(rsid_list)
        # Create the line with proper spacing:
        # 1 space after gene_id, 4 spaces after window_string
        write_line = f'{gene} {window_string}    {rsid_string}'
        f.write(f"{write_line}\n")
        print('\n')

print(f"Prediction MAGMA file created at: {MAGMA_file_pred}")
print(f"Total prediction genes written: {len(pred_merged_genes)}")

# Write the pip MAGMA file
with open(MAGMA_file_pip, 'w') as f:
    # Write the header lines
    f.write(f"{first_line}\n")
    f.write(f"{second_line}\n")
    # Iterate through genes and write each line
    for gene in tqdm(pip_merged_genes):
        num_variants = num_pip_gene_variants[num_pip_gene_variants['gene_id'] == gene]['counts'].values[0]
        TSS_position = gene_id_gene_name_df[gene_id_gene_name_df['gene_id'] == gene]['gene_TSS'].values[0]
        #get the n closest variants to TSS_position in annotation_df using the column BP
        variant_subset = annotation_df.copy(deep=True)
        variant_subset['distance_to_TSS'] = np.abs(variant_subset['BP'] - TSS_position)
        closest_variants = variant_subset.nsmallest(num_variants, 'distance_to_TSS')
        start_window = closest_variants['BP'].min() - 5
        end_window = closest_variants['BP'].max() + 5
        window_string = f'{chr}:{start_window}:{end_window}'
        rsid_list = closest_variants['SNP'].tolist()
        # Create rsid string with 4 spaces between each rsid
        rsid_string = '    '.join(rsid_list)
        # Create the line with proper spacing:
        # 1 space after gene_id, 4 spaces after window_string
        write_line = f'{gene} {window_string}    {rsid_string}'
        f.write(f"{write_line}\n")

print(f"PIP MAGMA file created at: {MAGMA_file_pip}")
print(f"Total PIP genes written: {len(pip_merged_genes)}")





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


def anti_join(x, y, on):
    """Return rows in x which are not present in y"""
    ans = pd.merge(left=x, right=y, how='left', indicator=True, on=on)
    ans = ans.loc[ans._merge == 'left_only', :].drop(columns='_merge')
    return ans


import pandas as pd
import requests
import json
from time import sleep



def get_gene_names_batch(ensembl_ids, batch_size=200):
    """
    Get gene names for Ensembl IDs using the Ensembl REST API
    """
    base_url = "https://rest.ensembl.org"
    # Remove duplicates and convert to list
    unique_ids = list(set(ensembl_ids))
    # Filter out None/NaN values
    unique_ids = [id for id in unique_ids if pd.notna(id) and id is not None]
    gene_name_map = {}
    # Process in batches to avoid overwhelming the API
    for i in range(0, len(unique_ids), batch_size):
        batch = unique_ids[i:i + batch_size]
        print(f"Processing batch {i // batch_size + 1}/{(len(unique_ids) - 1) // batch_size + 1}: {len(batch)} genes")
        # Prepare the POST request
        headers = {
            "Content-Type": "application/json",
            "Accept": "application/json"
        }
        # Create the request data
        request_data = {
            "ids": batch
        }
        try:
            # Make the API request
            response = requests.post(
                f"{base_url}/lookup/id",
                headers=headers,
                data=json.dumps(request_data),
                timeout=30  # Add timeout
            )
            print(f"Response status code: {response.status_code}")
            if response.status_code == 200:
                # Check if response has content
                if response.text.strip():
                    try:
                        data = response.json()
                        # Extract gene names
                        for gene_id, info in data.items():
                            if isinstance(info, dict):  # Make sure info is a dictionary
                                if 'display_name' in info:
                                    gene_name_map[gene_id] = info['display_name']
                                elif 'external_name' in info:
                                    gene_name_map[gene_id] = info['external_name']
                                else:
                                    gene_name_map[gene_id] = gene_id  # fallback to ID
                            else:
                                # Handle case where API returns error for specific gene
                                gene_name_map[gene_id] = gene_id
                    except json.JSONDecodeError as e:
                        print(f"JSON decode error: {e}")
                        print(f"Response content: {response.text[:200]}...")
                        # For failed batch, use gene IDs as names
                        for gene_id in batch:
                            gene_name_map[gene_id] = gene_id
                else:
                    print("Empty response received")
                    # For failed batch, use gene IDs as names
                    for gene_id in batch:
                        gene_name_map[gene_id] = gene_id
            else:
                print(f"API request failed with status code: {response.status_code}")
                print(f"Response content: {response.text[:200]}...")
                # For failed batch, use gene IDs as names
                for gene_id in batch:
                    gene_name_map[gene_id] = gene_id
        except requests.exceptions.RequestException as e:
            print(f"Request error: {e}")
            # For failed batch, use gene IDs as names
            for gene_id in batch:
                gene_name_map[gene_id] = gene_id
        except Exception as e:
            print(f"Unexpected error processing batch: {e}")
            # For failed batch, use gene IDs as names
            for gene_id in batch:
                gene_name_map[gene_id] = gene_id
        # Be respectful to the API - longer sleep for large batches
        sleep(0.5)
    return gene_name_map


def add_gene_names(df):
    """
    Add gene names to your dataframe
    """
    # Get unique Ensembl IDs, filtering out NaN values
    ensembl_ids = df['GENE'].dropna().unique()
    if len(ensembl_ids) == 0:
        print("No valid gene IDs found")
        df['GENE_NAME'] = df['GENE']  # Use gene ID as name
        return df
    # Get gene names
    print(f"Fetching gene names for {len(ensembl_ids)} unique genes...")
    gene_name_map = get_gene_names_batch(ensembl_ids)
    # Add gene names to dataframe
    df['GENE_NAME'] = df['GENE'].map(gene_name_map).fillna(df['GENE'])
    return df





# Load config
with open(os.path.join(os.path.dirname(__file__), '..', 'config.yaml'), 'r') as f:
    config = yaml.safe_load(f)

cohorts = ['Ast_mega_eQTL',
              'Exc_mega_eQTL',
              'Inh_mega_eQTL',
              'Mic_mega_eQTL' ,
              'Oli_mega_eQTL',
              'OPC_mega_eQTL']



list_genes_pred = []
list_genes_pip = []

num_genes_dict = {}

num_genes_total = 18700


bonferroni_val = 0.05 / num_genes_total

for cell_type in cohorts:
    print(f'Processing {cell_type}')
    predictions_dir = os.path.join(config['paths']['data_dir'], cell_type)
    predictions_df = dd.read_parquet(f'{predictions_dir}/predictions_parquet_catboost/predictions.parquet')
    num_genes = len(predictions_df['gene_id'].unique())
    print(f'Number of genes in {cell_type}: {num_genes}')
    MAGMA_dir = os.path.join(config['paths']['output_dir'], 'aggregate_results', 'MAGMA_knn', cell_type, 'MAGMA_output')
    print(f'Processing {cell_type}')
    pred_MAGMA = dd.read_csv(f'{MAGMA_dir}/{cell_type}_chr*_MAGMA_prediction.genes.out', sep='\s+').compute()
    pip_MAGMA = dd.read_csv(f'{MAGMA_dir}/{cell_type}_chr*_MAGMA_pip.genes.out', sep='\s+').compute()
    pred_MAGMA['num_genes'] = num_genes
    pip_MAGMA['num_genes'] = num_genes
    #calculate bonferroni correction using P column
    pred_MAGMA['bonferroni'] = 0.05 / num_genes_total
    pip_MAGMA['bonferroni'] = 0.05 / num_genes_total
    # Filter the dataframe to keep only the significant genes
    pred_MAGMA_significant = pred_MAGMA[pred_MAGMA['P'] < pred_MAGMA['bonferroni']]
    pip_MAGMA_significant = pip_MAGMA[pip_MAGMA['P'] < pip_MAGMA['bonferroni']]
    pred_MAGMA_significant['cell_type'] = cell_type
    pip_MAGMA_significant['cell_type'] = cell_type
    # Append the significant genes to the list
    list_genes_pred.append(pred_MAGMA_significant)
    list_genes_pip.append(pip_MAGMA_significant)
    num_genes_dict[cell_type] = num_genes




# Concatenate all the dataframes in the list into a single dataframe
pred_MAGMA_significant_df = pd.concat(list_genes_pred, ignore_index=True)
pip_MAGMA_significant_df = pd.concat(list_genes_pip, ignore_index=True)


pred_MAGMA_significant_df = add_gene_names(pred_MAGMA_significant_df)

pip_MAGMA_significant_df = add_gene_names(pip_MAGMA_significant_df)



aggregate_dir = os.path.join(config['paths']['output_dir'], 'aggregate_results', 'MAGMA_significant_genes')

if not os.path.exists(aggregate_dir):
    os.makedirs(aggregate_dir)





pred_MAGMA_significant_df.to_csv(f'{aggregate_dir}/MAGMA_knn_significant_genes_pred.txt', sep='\t', index=False)

pip_MAGMA_significant_df.to_csv(f'{aggregate_dir}/MAGMA_knn_significant_genes_pip.txt', sep='\t', index=False)


MAGMA_original_dir = os.path.join(config['paths']['output_dir'], 'aggregate_results', 'MAGMA', 'MAGMA_output_original', 'original_chr*_MAGMA_prediction.genes.out')


MAGMA_original_df = dd.read_csv(MAGMA_original_dir, sep='\s+').compute()

num_genes_unique = len(MAGMA_original_df['GENE'].unique())

MAGMA_original_df['bonferroni'] = 0.05 / num_genes_unique

MAGMA_original_df = MAGMA_original_df[MAGMA_original_df['P'] < MAGMA_original_df['bonferroni']]



gene_conversion_file = config['paths']['gene_info_file']

gene_conversion_df = pd.read_csv(gene_conversion_file, sep= '\t')




# Apply the conversion
gene_conversion_df['gene_id'] = gene_conversion_df['dbXrefs'].str.extract(r'Ensembl:(ENSG\d+)', expand=False)



gene_mapping = dict(zip(gene_conversion_df['GeneID'], gene_conversion_df['gene_id']))

MAGMA_original_df['ENSEMBL_GENE'] = MAGMA_original_df['GENE'].map(gene_mapping)

#rename ENSEMBL_GENE to GENE

MAGMA_original_df.drop(columns='GENE', inplace=True)

MAGMA_original_df.rename(columns = {'ENSEMBL_GENE':'GENE'}, inplace=True)



MAGMA_original_genes = MAGMA_original_df['GENE'].tolist()


pred_MAGMA_significant_df_not_overlap = pred_MAGMA_significant_df.loc[~pred_MAGMA_significant_df.GENE.isin(MAGMA_original_genes)]

pip_MAGMA_significant_df_not_overlap = pip_MAGMA_significant_df.loc[~pip_MAGMA_significant_df.GENE.isin(MAGMA_original_genes)]



pred_MAGMA_significant_df_not_overlap.to_csv(f'{aggregate_dir}/MAGMA_knn_significant_genes_pred_no_overlap.txt', sep='\t', index=False)

pip_MAGMA_significant_df_not_overlap.to_csv(f'{aggregate_dir}/MAGMA_knn_significant_genes_pip_no_overlap.txt', sep='\t', index=False)


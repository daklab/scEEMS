"""
List the genes of the FunGen-xQTL SuSiE fine-mapping exports (one Fungen_xQTL.<gene_id>.cis_results_db.export.rds
file per gene), the genes get_vars_pips.R processes.

The exports are not distributed with this repository or the data release; this step documents how the
released training data were made.

Usage:
    python create_job_list.py

Input:   {finemapping_rds_dir}/Fungen_xQTL.<gene_id>.cis_results_db.export.rds
Output:  {susie_dir}/genes_list.txt, one gene ID per line (sorted); a gene's line number is its index for
         get_vars_pips.R
"""

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

input_dir = path("finemapping_rds_dir")
out_dir = path("susie_dir")
os.makedirs(out_dir, exist_ok=True)

files = []
for root, _, filenames in os.walk(input_dir):
    for filename in filenames:
        if filename.endswith('cis_results_db.export.rds') and filename.startswith('Fungen_xQTL'):
            files.append(filename)

if not files:
    raise ValueError(
        f"No Fungen_xQTL.<gene_id>.cis_results_db.export.rds files under {input_dir} (finemapping_rds_dir). "
        "These fine-mapping exports are not part of the data release."
    )

genes = []
for filename in files:
    parts = filename.split('.')
    if len(parts) > 2:
        genes.append(parts[1])

unique_genes = pd.Series(sorted(set(genes)))

out_file = os.path.join(out_dir, 'genes_list.txt')
unique_genes.to_csv(out_file, index=False, header=None)

print(f"Found {len(unique_genes)} unique genes. Saved to {out_file}")

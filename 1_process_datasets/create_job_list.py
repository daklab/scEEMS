"""
Generate a list of unique genes from downloaded SuSiE fine-mapping results.

Parses filenames of downloaded RDS files to extract gene identifiers
and saves a deduplicated gene list for downstream processing.

Note:
    This script expects internal raw RDS files under release_04_2024.
    Those Step 1 raw inputs are not publicly distributed in this repository.

Usage:
    python create_job_list.py

Requires config.yaml with paths.data_dir set.
"""

import os
import pandas as pd
import yaml


# Load configuration
with open('../config.yaml', 'r') as f:
    config = yaml.safe_load(f)

input_dir = os.path.join(config['paths']['data_dir'], 'release_04_2024')
out_dir = os.path.dirname(os.path.abspath(__file__))

files = []
for root, _, filenames in os.walk(input_dir):
    for filename in filenames:
        if filename.endswith('cis_results_db.export.rds') and filename.startswith('Fungen_xQTL'):
            files.append(filename)

if not files:
    raise ValueError(
        f"No fine-mapping RDS files found under {input_dir}. "
        "Run download_data.py first and confirm Synapse data was downloaded."
    )

genes = []
for filename in files:
    parts = filename.split('.')
    if len(parts) > 2:
        genes.append(parts[1])

unique_genes = pd.Series(sorted(set(genes)))

# Save unique genes to file
unique_genes.to_csv(os.path.join(out_dir, 'genes_list.txt'), index=False, header=None)

print(f"Found {len(unique_genes)} unique genes. Saved to genes_list.txt")

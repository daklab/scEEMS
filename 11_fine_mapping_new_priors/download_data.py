"""
Download reference panel data from Synapse for fine-mapping.

Downloads genotype reference panel files needed for SuSiE fine-mapping
with scEEMS-informed priors.

Usage:
    python download_data.py

Requires config.yaml with synapse_token set.
"""

import os
import synapseclient
import synapseutils
import yaml

# Load configuration
with open('../config.yaml', 'r') as f:
    config = yaml.safe_load(f)

synapse_token = config['credentials']['synapse_token']
if not synapse_token:
    raise ValueError("Synapse token not set in config.yaml. "
                     "Get a personal access token from https://www.synapse.org/#!PersonalAccessTokens:")

fine_mapping_dir = config['paths']['fine_mapping_dir']
out_dir = os.path.join(fine_mapping_dir, 'reference')
os.makedirs(out_dir, exist_ok=True)

# Connect to Synapse
syn = synapseclient.Synapse()
syn.login(authToken=synapse_token)

# Download reference panel
files = synapseutils.syncFromSynapse(syn, 'syn53171227', path=out_dir)

print(f"Reference panel downloaded to {out_dir}")

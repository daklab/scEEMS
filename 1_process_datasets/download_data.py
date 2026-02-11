"""
Step 1 data download status for public repository.

This repository does not expose the raw fine-mapping RDS inputs used in Step 1.
The Step 1 scripts are retained to document how public train/test parquet files were
created internally, but these raw source files are not distributed.

For public downloads, use:
    python ../download_synapse_data.py --resource model_training
    python ../download_synapse_data.py --resource predictions
"""

raise RuntimeError(
    "Step 1 raw RDS inputs are not publicly available in this repository. "
    "Use download_synapse_data.py to fetch public model_training/predictions data from Synapse."
)

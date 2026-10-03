"""
Step 1 reads the FunGen-xQTL SuSiE fine-mapping exports, which are not distributed with this repository
or the scEEMS data release. The step 1 scripts document how the released training data were made.

The data release, which steps 5-11 start from, is downloaded with the script at the top of the repository:
    python ../download_synapse_data.py --resource model_training predictions fine_mapping
"""

raise RuntimeError(
    "The step 1 fine-mapping exports are not publicly available. "
    "Use download_synapse_data.py at the top of the repository to download the scEEMS data release."
)

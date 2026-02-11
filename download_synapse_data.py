"""
Download public scEEMS datasets from Synapse.

Public resources:
- model_training (syn72248754): train/test parquet files + supporting files
- predictions (syn71338354): per-cell-type prediction outputs

Usage examples:
  python download_synapse_data.py --resource all
  python download_synapse_data.py --resource model_training
  python download_synapse_data.py --resource predictions --output_dir /path/to/data

Authentication:
- If credentials.synapse_token is set in config.yaml, token login is used.
- Otherwise, Synapse interactive login is used.
"""

import argparse
import os
import synapseclient
import synapseutils
import yaml


DEFAULT_MODEL_TRAINING_ID = "syn72248754"
DEFAULT_PREDICTIONS_ID = "syn71338354"


def load_config():
    with open("config.yaml", "r") as f:
        return yaml.safe_load(f)


def login_synapse(config):
    syn = synapseclient.Synapse()
    token = config.get("credentials", {}).get("synapse_token", "")
    if token:
        syn.login(authToken=token)
    else:
        syn.login()
    return syn


def sync_entity(syn, entity_id, output_dir, label):
    target_dir = os.path.join(output_dir, label)
    os.makedirs(target_dir, exist_ok=True)
    print(f"Syncing {label} ({entity_id}) -> {target_dir}")
    synapseutils.syncFromSynapse(syn, entity_id, path=target_dir)


def main():
    parser = argparse.ArgumentParser(description="Download public scEEMS data from Synapse")
    parser.add_argument(
        "--resource",
        choices=["model_training", "predictions", "all"],
        default="all",
        help="Which public resource to download",
    )
    parser.add_argument(
        "--output_dir",
        default=None,
        help="Output directory (default: paths.data_dir/synapse_public)",
    )
    args = parser.parse_args()

    config = load_config()
    data_sources = config.get("data_sources", {})

    model_training_id = data_sources.get("model_training_synapse_id", DEFAULT_MODEL_TRAINING_ID)
    predictions_id = data_sources.get("predictions_synapse_id", DEFAULT_PREDICTIONS_ID)

    base_data_dir = config.get("paths", {}).get("data_dir", ".")
    output_dir = args.output_dir or os.path.join(base_data_dir, "synapse_public")
    os.makedirs(output_dir, exist_ok=True)

    syn = login_synapse(config)

    if args.resource in ["model_training", "all"]:
        sync_entity(syn, model_training_id, output_dir, "model_training")

    if args.resource in ["predictions", "all"]:
        sync_entity(syn, predictions_id, output_dir, "predictions")

    print("Synapse download complete.")


if __name__ == "__main__":
    main()

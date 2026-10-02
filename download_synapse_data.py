"""
Download the scEEMS data release from Synapse into paths.release_dir (config.yaml), keeping the release's
folder layout, which is where steps 5-12 look for it.

The release folder (default syn69670587) contains:
  model_training/  training and test data, GPN-STAR scores, selected feature weights, supporting files
  predictions/     scEEMS predictions, one tabix-indexed TSV per cell type and chromosome
  fine_mapping/    credible sets of the five fine-mapping priors
Folders are found by name, so the script does not depend on the Synapse IDs of the subfolders.

Usage:
  python download_synapse_data.py --dry-run                                  # list what would be downloaded
  python download_synapse_data.py                                            # everything
  python download_synapse_data.py --resource model_training --cell-type Mic  # step 5 quick start
  python download_synapse_data.py --resource predictions --cell-type Mic Ast

Authentication: a Synapse personal access token from config.yaml (credentials.synapse_token) or the
SYNAPSE_AUTH_TOKEN environment variable; otherwise synapseclient's own login (e.g. ~/.synapseConfig).
Listing (--dry-run) works without logging in.
"""
import argparse
import os
import re

import synapseclient
import yaml

ROOT_ID = "syn69670587"
RESOURCES = ["model_training", "predictions", "fine_mapping"]
CELL_DIR = re.compile(r"^(Ast|Exc|Inh|Mic|Oli|OPC)_mega_eQTL$")


def load_config():
    path = os.environ.get("SCEEMS_CONFIG", os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.yaml"))
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def walk(syn, folder_id, rel, cells):
    """Yield (relative path, file entity id) under a folder, skipping cell type folders not in `cells`."""
    for ch in syn.getChildren(folder_id, includeTypes=["folder", "file"]):
        name, path = ch["name"], os.path.join(rel, ch["name"])
        if ch["type"].endswith("Folder"):
            m = CELL_DIR.match(name)
            if cells and m and m.group(1) not in cells:
                continue
            yield from walk(syn, ch["id"], path, cells)
        else:
            yield path, ch["id"]


def main():
    ap = argparse.ArgumentParser(description="Download the scEEMS data release from Synapse")
    ap.add_argument("--resource", nargs="+", choices=RESOURCES, default=RESOURCES)
    ap.add_argument("--cell-type", nargs="+", choices=["Ast", "Exc", "Inh", "Mic", "Oli", "OPC"],
                    help="only these cell types (shared files are always included)")
    ap.add_argument("--root", default=ROOT_ID, help=f"Synapse ID of the release folder (default {ROOT_ID})")
    ap.add_argument("--output-dir", help="destination (default: paths.release_dir in config.yaml)")
    ap.add_argument("--dry-run", action="store_true", help="list the files, download nothing")
    args = ap.parse_args()

    config = {} if (args.output_dir and args.dry_run) else load_config()
    out_dir = args.output_dir or config.get("paths", {}).get("release_dir")
    if not out_dir:
        raise SystemExit("set paths.release_dir in config.yaml or pass --output-dir")

    syn = synapseclient.Synapse(silent=True)
    token = os.environ.get("SYNAPSE_AUTH_TOKEN") or (config.get("credentials") or {}).get("synapse_token")
    if not args.dry_run:
        syn.login(authToken=token) if token else syn.login()

    folders = {ch["name"]: ch["id"] for ch in syn.getChildren(args.root, includeTypes=["folder"])}
    for res in args.resource:
        if res not in folders:
            print(f"{res}: not found under {args.root} (folders: {', '.join(sorted(folders))})")
            continue
        n = 0
        for rel, eid in walk(syn, folders[res], res, set(args.cell_type or [])):
            n += 1
            if args.dry_run:
                print(rel)
                continue
            dest = os.path.join(out_dir, os.path.dirname(rel))
            os.makedirs(dest, exist_ok=True)
            syn.get(eid, downloadLocation=dest, ifcollision="overwrite.local")
        print(f"{res}: {n} files {'listed' if args.dry_run else 'downloaded to ' + os.path.join(out_dir, res)}")


if __name__ == "__main__":
    main()

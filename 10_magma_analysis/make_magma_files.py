#!/usr/bin/env python
"""
Build the eMAGMA gene annotations (.genes.annot) of one cell type and one chromosome. Each gene's SNP set
is the set of variants linked to it in one of two ways:

  prediction  variants the scEEMS model calls eQTLs of the gene (pred_prob > tau*, the cell type's
              threshold from the S-LDSC analysis of step 9), from the main variant set and from the
              non-European panel (score_noneur.py). Each variant is then assigned to ONE gene: the gene
              with the most cell type-matched PLAC-seq enhancer-promoter interactions overlapping the
              variant, ties broken by the highest ABC score (ABC >= 0.005) in the matched cell type;
              variants supported by neither are dropped. Astrocytes have no PLAC-seq data, so ABC alone.
  pip         variants with fine-mapping PIP > 0.10 for the gene (all genes kept)

Only variants of the European S-LDSC baseline annotation (which carries the SNP IDs of the reference
panel) are used from the main variant set.

usage:   python make_magma_files.py CHR COHORT
outputs: {magma_dir}/{cohort}/{prediction,pip}/{cohort}_chr{CHR}_MAGMA.genes.annot
         {magma_dir}/{cohort}/variant_gene_lists/MAGMA_variant_gene_list_{cohort}_chr{CHR}_{pred,pip}.tsv
"""
import json
import os
import sys

import dask
import dask.dataframe as dd
import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import cohort_name, path

tqdm.pandas()
dask.config.set({"temporary_directory": path("scratch_dir")})


def overlap_placseq_simple_fix(df, placseq_df):
    """Count PLAC-seq enhancer-promoter interactions per (variant, gene) by cell-type interaction column."""
    variant_ids = df["variant_id"].values
    gene_ids = df["gene_id"].values
    n_rows = len(df)
    positions = np.array([int(vid.split(":")[1]) for vid in variant_ids])
    microglia_counts = np.zeros(n_rows, dtype=np.int32)
    oligodendrocyte_counts = np.zeros(n_rows, dtype=np.int32)
    neuron_counts = np.zeros(n_rows, dtype=np.int32)
    cell_type_dict = {
        "microglia": ("PU1_enhancer_interactions", microglia_counts),
        "oligodendrocyte": ("Olig2_enhancer_interactions", oligodendrocyte_counts),
        "neuron": ("NeuN_enhancer_interactions", neuron_counts),
    }
    placseq_lookup = {}
    for interaction_type in ["PU1_enhancer_interactions", "Olig2_enhancer_interactions", "NeuN_enhancer_interactions"]:
        type_data = placseq_df[placseq_df["interaction_type"] == interaction_type]
        if not type_data.empty:
            grouped = type_data.groupby("gene_id")[["start", "end"]].apply(
                lambda x: (x["start"].values, x["end"].values)).to_dict()
            placseq_lookup[interaction_type] = grouped
    for i in tqdm(range(n_rows), desc="Processing variants"):
        gene_id = gene_ids[i]
        pos = positions[i]
        for cell_type, (interaction_type, count_array) in cell_type_dict.items():
            if interaction_type in placseq_lookup and gene_id in placseq_lookup[interaction_type]:
                starts, ends = placseq_lookup[interaction_type][gene_id]
                overlaps = (starts <= pos) & (ends >= pos)
                count_array[i] = np.sum(overlaps)
    df["microglia_placseq_count"] = microglia_counts
    df["oligodendrocyte_placseq_count"] = oligodendrocyte_counts
    df["neuron_placseq_count"] = neuron_counts
    return df


def get_max_abc_score(chr, bp, abc_df):
    mask = (abc_df["chr"] == chr) & (abc_df["start"] <= bp) & (abc_df["end"] >= bp)
    if mask.any():
        return abc_df.loc[mask, "ABC.Score"].max()
    return 0.0


def process_abc_scores(training_data_df, abc_dfs, abc_names):
    for abc_column in abc_names:
        training_data_df[f"ABC_score_{abc_column}"] = 0.0
    for abc_df, abc_name in zip(abc_dfs, abc_names):
        for idx, row in tqdm(training_data_df.iterrows()):
            gene_id_val = row["gene_id"]
            abc_df_subset = abc_df[abc_df["TargetGeneEnsembl_ID"] == gene_id_val]
            max_abc_score = get_max_abc_score(row["CHR"], row["BP"], abc_df_subset)
            training_data_df.loc[idx, f"ABC_score_{abc_name}"] = max_abc_score
    return training_data_df


def sort_and_filter_by_cell_type(pred_merged, cell_type, cell_type_abc_plac_dict):
    """Keep the top gene per variant by PLAC-seq (primary) then ABC (secondary); drop unsupported variants."""
    config = cell_type_abc_plac_dict[cell_type]
    placseq_col = config["placseq"]
    abc_col = config["abc"]
    df = pred_merged.copy()
    if placseq_col is None or placseq_col not in df.columns:
        df["sort_priority"] = df[abc_col]
        df["has_placseq"] = False
    else:
        df["sort_priority"] = df[placseq_col]
        df["has_placseq"] = True
        df["secondary_sort"] = df[abc_col]
    if placseq_col is None or placseq_col not in df.columns:
        df_filtered = df[df[abc_col] != 0].copy()
    else:
        df_filtered = df[~((df[placseq_col] == 0) & (df[abc_col] == 0))].copy()
    if placseq_col is None or placseq_col not in df_filtered.columns:
        df_sorted = df_filtered.sort_values(["variant_id", "sort_priority"], ascending=[True, False])
    else:
        df_sorted = df_filtered.sort_values(
            ["variant_id", "sort_priority", "secondary_sort"], ascending=[True, False, False])
    result = df_sorted.groupby("variant_id").first().reset_index()
    columns_to_drop = ["sort_priority", "has_placseq"]
    if "secondary_sort" in result.columns:
        columns_to_drop.append("secondary_sort")
    return result.drop(columns=columns_to_drop)


def write_annot(path, genes, merged, chr):
    """MAGMA .genes.annot: header + one line per gene 'gene chr:start:end<4sp>rsid<4sp>rsid...'."""
    with open(path, "w") as f:
        f.write("# window_up = 0\n")
        f.write("# window_down = 0\n")
        for gene in genes:
            g = merged[merged["gene_id"] == gene]
            window_string = f"{chr}:{g['BP'].min() - 5}:{g['BP'].max() + 5}"
            rsid_string = "    ".join(g["SNP"].tolist())
            f.write(f"{gene} {window_string}    {rsid_string}\n")


chr = sys.argv[1]
cell_type = cohort_name(sys.argv[2])

# cell type -> PLAC-seq interaction count column and ABC cell type used to assign each variant to one gene
cell_type_abc_plac_dict = {
    "Mic_mega_eQTL": {"placseq": "microglia_placseq_count", "abc": "ABC_score_microglia", "abc_name": "microglia"},
    "Ast_mega_eQTL": {"placseq": None, "abc": "ABC_score_astrocyte", "abc_name": "astrocyte"},
    "Exc_mega_eQTL": {"placseq": "neuron_placseq_count", "abc": "ABC_score_neuron", "abc_name": "neuron"},
    "Inh_mega_eQTL": {"placseq": "neuron_placseq_count", "abc": "ABC_score_neuron", "abc_name": "neuron"},
    "Oli_mega_eQTL": {"placseq": "oligodendrocyte_placseq_count", "abc": "ABC_score_oligodendrocyte", "abc_name": "oligodendrocyte"},
    "OPC_mega_eQTL": {"placseq": "oligodendrocyte_placseq_count", "abc": "ABC_score_oligodendrocyte", "abc_name": "oligodendrocyte"},
}
cell_cfg = cell_type_abc_plac_dict[cell_type]

with open(f"{path('aggregate_dir')}/tau_star.json") as fh:
    pred_prob_threshold = json.load(fh)[cell_type]                # tau* (step 9)
pip_threshold = 0.10

annotation_df = pd.read_csv(f"{path('baseline_annot_dir')}/baseline_chr{chr}.annot.gz", sep="\t",
                            usecols=["CHR", "SNP", "BP", "A1", "A2"])
annotation_df = annotation_df.rename(columns={"A1": "alt", "A2": "ref"})

write_MAGMA_dir = f"{path('magma_dir')}/{cell_type}"
write_MAGMA_dir_pred = f"{write_MAGMA_dir}/prediction"
write_MAGMA_dir_pip = f"{write_MAGMA_dir}/pip"
write_MAGMA_dir_variant_gene = f"{write_MAGMA_dir}/variant_gene_lists"
for d in (write_MAGMA_dir_pred, write_MAGMA_dir_pip, write_MAGMA_dir_variant_gene):
    os.makedirs(d, exist_ok=True)
MAGMA_file_pred = f"{write_MAGMA_dir_pred}/{cell_type}_chr{chr}_MAGMA.genes.annot"
MAGMA_file_pip = f"{write_MAGMA_dir_pip}/{cell_type}_chr{chr}_MAGMA.genes.annot"
for p in (MAGMA_file_pred, MAGMA_file_pip):
    if os.path.exists(p):
        os.remove(p)

# --- scEEMS predictions for this chromosome ---
predictions_df = dd.read_parquet(
    f"{path('predictions_parquet_dir', cohort=cell_type)}/weighted_full/predictions.parquet")
predictions_df = predictions_df[predictions_df["chr"] == f"chr{chr}"].persist()
top_pred = predictions_df[predictions_df["pred_prob"] > pred_prob_threshold].compute()
top_pip = predictions_df[predictions_df["pip"] > pip_threshold].compute()
for d in (top_pred, top_pip):
    d.rename(columns={"pos": "BP", "chr": "CHR"}, inplace=True)
    d["CHR"] = d["CHR"].str.replace("chr", "").astype(int)
pred_merged = pd.merge(annotation_df, top_pred, on=["CHR", "BP", "ref", "alt"], how="inner")
pip_merged = pd.merge(annotation_df, top_pip, on=["CHR", "BP", "ref", "alt"], how="inner")

# --- non-European panel variants passing tau*, added to the prediction links ---
noneur_file = (f"{path('magma_dir')}/noneur_predictions/{cell_type}/"
               f"MAGMA_predictions_{cell_type}_chr{chr}.tsv.gz")
pred_noneur = pd.read_csv(noneur_file, sep="\t")
pred_noneur["CHR"] = pred_noneur["CHR"].str.replace("chr", "")
pred_noneur["variant_id"] = (pred_noneur["CHR"].astype(str) + ":" + pred_noneur["BP"].astype(str)
                             + ":" + pred_noneur["REF"] + ":" + pred_noneur["ALT"])
pred_noneur = pred_noneur.rename(columns={"REF": "ref", "ALT": "alt"})
pred_noneur = pred_noneur[pred_noneur["pred_prob"] > pred_prob_threshold]
pred_merged = pd.concat([pred_merged, pred_noneur], ignore_index=True)

pred_merged = pred_merged.groupby(["CHR", "SNP", "BP", "ref", "alt", "gene_id"]).first().reset_index()

# --- ABC (>= 0.005) + PLAC-seq: one gene per predicted variant ---
abc_name = cell_cfg["abc_name"]
abc_df_single = pd.read_csv(f"{path('abc_data_dir')}/ABC_results_{abc_name}_v2/{abc_name}/Predictions/"
                            f"EnhancerPredictionsAllPutative.tsv.gz", sep="\t")
abc_df_single = abc_df_single[abc_df_single["ABC.Score"] >= 0.005]
abc_df_single = abc_df_single[abc_df_single["chr"].str.replace("chr", "") == str(chr)]
abc_df_single["chr"] = abc_df_single["chr"].str.replace("chr", "").astype(int)

placseq_df = pd.read_csv(path("placseq_file"), sep="\t")

pred_merged = process_abc_scores(pred_merged, [abc_df_single], [abc_name])
pred_merged = overlap_placseq_simple_fix(pred_merged, placseq_df)
pred_merged = sort_and_filter_by_cell_type(pred_merged, cell_type, cell_type_abc_plac_dict)

pred_merged.to_csv(f"{write_MAGMA_dir_variant_gene}/MAGMA_variant_gene_list_{cell_type}_chr{chr}_pred.tsv",
                   sep="\t", index=False)
pip_merged.to_csv(f"{write_MAGMA_dir_variant_gene}/MAGMA_variant_gene_list_{cell_type}_chr{chr}_pip.tsv",
                  sep="\t", index=False)

write_annot(MAGMA_file_pred, pred_merged["gene_id"].unique(), pred_merged, chr)
print(f"Prediction MAGMA file: {MAGMA_file_pred} ({pred_merged['gene_id'].nunique()} genes)")
write_annot(MAGMA_file_pip, pip_merged["gene_id"].unique(), pip_merged, chr)
print(f"PIP MAGMA file: {MAGMA_file_pip} ({pip_merged['gene_id'].nunique()} genes)")

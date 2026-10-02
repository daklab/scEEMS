#!/usr/bin/env python
"""
Distance-matched control annotations for eMAGMA, for one cell type and one chromosome: for every gene of
the prediction and pip annotations of make_magma_files.py, the same number of variants, chosen as the
variants of the European S-LDSC baseline annotation closest to the gene's TSS instead of the variants the
model or fine-mapping selected. Genes significant with these annotations are significant through proximity
to the TSS and SNP count alone, which separates that from the model's choice of variants.

usage:   python make_magma_files_knn.py CHR COHORT        (after make_magma_files.py for the same CHR)
output:  {magma_dir}/knn/{cohort}/{prediction,pip}/{cohort}_chr{CHR}_MAGMA.genes.annot
"""
import os
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import cohort_name, path

tqdm.pandas()

chr = sys.argv[1]
cell_type = cohort_name(sys.argv[2])

annotation_df = pd.read_csv(f"{path('baseline_annot_dir')}/baseline_chr{chr}.annot.gz", sep="\t",
                            usecols=["CHR", "SNP", "BP", "A1", "A2"])
annotation_df = annotation_df.rename(columns={"A1": "alt", "A2": "ref"})

gene_id_gene_name_df = pd.read_csv(f"{path('abc_data_dir')}/ABC_gene_id_name_mapping.csv", sep=",")
gene_id_gene_name_df_list = gene_id_gene_name_df["gene_id"].unique().tolist()

write_MAGMA_dir = f"{path('magma_dir')}/knn/{cell_type}"
write_MAGMA_dir_pred = f"{write_MAGMA_dir}/prediction"
write_MAGMA_dir_pip = f"{write_MAGMA_dir}/pip"
for d in (write_MAGMA_dir_pred, write_MAGMA_dir_pip):
    os.makedirs(d, exist_ok=True)
MAGMA_file_pred = f"{write_MAGMA_dir_pred}/{cell_type}_chr{chr}_MAGMA.genes.annot"
MAGMA_file_pip = f"{write_MAGMA_dir_pip}/{cell_type}_chr{chr}_MAGMA.genes.annot"
for p in (MAGMA_file_pred, MAGMA_file_pip):
    if os.path.exists(p):
        os.remove(p)

vg_dir = f"{path('magma_dir')}/{cell_type}/variant_gene_lists"
pred_merged = pd.read_csv(f"{vg_dir}/MAGMA_variant_gene_list_{cell_type}_chr{chr}_pred.tsv", sep="\t")
pip_merged = pd.read_csv(f"{vg_dir}/MAGMA_variant_gene_list_{cell_type}_chr{chr}_pip.tsv", sep="\t")


def write_knn_annot(path, merged):
    counts = merged.groupby("gene_id").size().reset_index(name="counts")
    genes = [g for g in merged["gene_id"].unique() if g in gene_id_gene_name_df_list]
    with open(path, "w") as f:
        f.write("# window_up = 0\n")
        f.write("# window_down = 0\n")
        for gene in tqdm(genes):
            num_variants = counts.loc[counts["gene_id"] == gene, "counts"].values[0]
            tss = gene_id_gene_name_df.loc[gene_id_gene_name_df["gene_id"] == gene, "gene_TSS"].values[0]
            vs = annotation_df.copy()
            vs["distance_to_TSS"] = np.abs(vs["BP"] - tss)
            closest = vs.nsmallest(num_variants, "distance_to_TSS")
            window_string = f"{chr}:{closest['BP'].min() - 5}:{closest['BP'].max() + 5}"
            rsid_string = "    ".join(closest["SNP"].tolist())
            f.write(f"{gene} {window_string}    {rsid_string}\n")
    return len(genes)


n_pred = write_knn_annot(MAGMA_file_pred, pred_merged)
print(f"kNN prediction MAGMA file: {MAGMA_file_pred} ({n_pred} genes)")
n_pip = write_knn_annot(MAGMA_file_pip, pip_merged)
print(f"kNN PIP MAGMA file: {MAGMA_file_pip} ({n_pip} genes)")

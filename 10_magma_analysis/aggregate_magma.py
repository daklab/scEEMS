#!/usr/bin/env python
"""
Collect the MAGMA gene results (run_magma.py) of all cell types into the tables of the manuscript.

  MAGMA_significant_genes_{pred,pip}.txt          eMAGMA, European GWAS: (gene, cell type) pairs with
                                                  P < 0.05 / 18,700 (Bonferroni over ~18,700 genes)
  MAGMA_significant_genes_{pred,pip}_no_overlap.txt   the same, without genes significant in the
                                                  positional MAGMA baseline
  MAGMA_knn_significant_genes_{pred,pip}[_no_overlap].txt   the same for the TSS-distance control
  MAGMA_multiancestry_genes_combined{,_pip}.txt   every eMAGMA gene result with its European (P_EUR) and
                                                  ADGC African (P_AFR), admixed American (P_AMR) and East
                                                  Asian (P_EAS) p-values (1 where a gene was not tested)
                                                  and its TSS
  MAGMA_positional_genes_combined.txt             positional MAGMA (gene bodies), genes as Ensembl IDs

Gene symbols come from the NCBI gene_info file (current HGNC symbols); the ABC gene mapping is used only
for Ensembl IDs NCBI does not cover, because it carries some obsolete aliases (e.g. CLTH for PICALM).
The NCBI file also maps the positional results' Entrez IDs to Ensembl IDs.

usage:   python aggregate_magma.py              (after every run_magma.py job has finished)
outputs: {aggregate_dir}/magma/*.txt
"""
import os
import sys

import dask.dataframe as dd
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import COHORTS, path

MAGMA = path("magma_dir")
OUT = f"{path('aggregate_dir')}/magma"
NUM_GENES_TOTAL = 18700
BONF = 0.05 / NUM_GENES_TOTAL              # 2.67e-6
os.makedirs(OUT, exist_ok=True)

gmap = pd.read_csv(f"{path('abc_data_dir')}/ABC_gene_id_name_mapping.csv", sep=",")
GENE_TSS = gmap[["gene_id", "gene_TSS"]].drop_duplicates().reset_index(drop=True)

conv = pd.read_csv(path("gene_info_file"), sep="\t")
conv["gene_id"] = conv["dbXrefs"].str.extract(r"Ensembl:(ENSG\d+)", expand=False)
entrez2ens = dict(zip(conv["GeneID"], conv["gene_id"]))
_ncbi = conv.dropna(subset=["gene_id"]).drop_duplicates("gene_id")
GENE_NAME = dict(zip(gmap["gene_id"], gmap.get("gene_name", gmap["gene_id"])))   # ABC fallback
GENE_NAME.update(dict(zip(_ncbi["gene_id"], _ncbi["Symbol"])))                   # NCBI takes priority


def load_out(subdir, cohort, arm, root=MAGMA):
    """All chromosomes' .genes.out of one cell type and annotation (arm = 'prediction' or 'pip')."""
    df = dd.read_csv(f"{root}/{cohort}/{subdir}/{cohort}_chr*_MAGMA_{arm}.genes.out", sep=r"\s+").compute()
    df["cell_type"] = cohort
    return df


def significant(root):
    """(prediction, pip) -> all gene results of the six cell types, flagged at the Bonferroni threshold."""
    out = {}
    for arm in ("prediction", "pip"):
        rows = []
        for coh in COHORTS:
            df = load_out("MAGMA_output", coh, arm, root=root)
            df["bonferroni"] = BONF
            df["significant"] = np.where(df["P"] < BONF, 1, 0)
            rows.append(df)
        tab = pd.concat(rows, ignore_index=True)
        tab["GENE_NAME"] = tab["GENE"].map(GENE_NAME).fillna(tab["GENE"])
        out[arm] = tab
    return out


# ---- positional MAGMA: Entrez -> Ensembl; its significant genes define the no-overlap filter ----
pos = dd.read_csv(f"{MAGMA}/positional/MAGMA_output/positional_chr*_MAGMA.genes.out", sep=r"\s+").compute()
pos["GENE"] = pos["GENE"].map(entrez2ens)
pos["bonferroni"] = 0.05 / pos["GENE"].nunique()
pos_sig_genes = pos.loc[pos["P"] < pos["bonferroni"], "GENE"].dropna().tolist()
pos["GENE_NAME"] = pos["GENE"].map(GENE_NAME).fillna(pos["GENE"])
pos.to_csv(f"{OUT}/MAGMA_positional_genes_combined.txt", sep="\t", index=False)

# ---- eMAGMA and the TSS-distance control, European GWAS ----
for root, prefix in ((MAGMA, "MAGMA"), (f"{MAGMA}/knn", "MAGMA_knn")):
    for arm, tab in significant(root).items():
        tag = "pred" if arm == "prediction" else "pip"
        sig = tab[tab["significant"] == 1].copy()
        sig.to_csv(f"{OUT}/{prefix}_significant_genes_{tag}.txt", sep="\t", index=False)
        sig[~sig["GENE"].isin(pos_sig_genes)].to_csv(f"{OUT}/{prefix}_significant_genes_{tag}_no_overlap.txt",
                                                    sep="\t", index=False)
        print(f"{prefix} {tag}: {len(sig)} significant (gene, cell type) pairs, "
              f"{(~sig['GENE'].isin(pos_sig_genes)).sum()} not significant in positional MAGMA")


# ---- eMAGMA, European and ADGC AFR/AMR/EAS p-values per gene ----
def combine(arm):
    def per_pop(subdir, pcol):
        return pd.concat([load_out(subdir, coh, arm)[["GENE", "CHR", "P", "cell_type"]].rename(columns={"P": pcol})
                          for coh in COHORTS], ignore_index=True)

    out = per_pop("MAGMA_output", "P_EUR")
    for pop in ("AFR", "AMR", "EAS"):
        out = out.merge(per_pop(f"MAGMA_output_ADGC_{pop}", f"P_{pop}"), on=["GENE", "CHR", "cell_type"], how="left")
    out = out.merge(GENE_TSS, left_on="GENE", right_on="gene_id", how="left")
    out["GENE_NAME"] = out["GENE"].map(GENE_NAME).fillna(out["GENE"])
    for c in ("P_EUR", "P_AFR", "P_AMR", "P_EAS"):
        out[c] = out[c].fillna(1.0)
    return out


combine("prediction").to_csv(f"{OUT}/MAGMA_multiancestry_genes_combined.txt", sep="\t", index=False)
combine("pip").to_csv(f"{OUT}/MAGMA_multiancestry_genes_combined_pip.txt", sep="\t", index=False)
print(f"wrote {OUT}")

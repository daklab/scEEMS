#!/usr/bin/env python
"""
Aggregate the cross-cell-type credible-set colocalizations of coloc_crosscell.R over all genes.

A credible set is shared (1) if any credible set of the same gene in any other cell type, under the same
prior, colocalizes with it at PP.H4 > 0.8; otherwise 0. The denominator is every credible set of that cell
type under that prior (the "focal credible set" rows), so a set whose gene has no credible set in any
other cell type counts as not shared. pct_cs_shared_with_partner restricts the denominator to sets with at
least one coloc test, for reference.

Protein-coding genes only (protein_coding_genes.txt from make_gene_lists.py). The scEEMS priors exist only
for genes with scEEMS predictions, which are protein-coding, so a scEEMS prior is never fitted for a
lncRNA, while the uniform prior is. Left unrestricted, those genes would read as eGenes the scEEMS priors
lost, when the prior was never evaluated there; restricting makes every prior compare on the same genes.

Credible-set size, top PIP and purity are added from the step 11 credible-set tables
({aggregate_dir}/finemapping_comparison/{cell}_{prior}_cs.tsv) when they exist, and the credible-set counts are checked
against them.

usage:   python aggregate_crosscell.py
outputs (in aggregate_dir):
  crosscell_coloc_pairs.tsv   every row of the per-gene tables (coloc tests and note rows)
  crosscell_cs_status.tsv     one row per credible set: shared (0/1), max PP.H4, the cell types it is shared with
  crosscell_egene_status.tsv  one row per (prior, cell type, gene): credible sets, shared credible sets, cell types
  crosscell_summary.tsv       per (prior, cell type): % credible sets shared, % eGenes shared
"""
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

XCD = path("crosscell_dir")
AGG = path("aggregate_dir")
CMP = os.path.join(AGG, "finemapping_comparison")
H4 = 0.8
KEY = ["prior", "cell", "gene_id", "idx"]
os.makedirs(AGG, exist_ok=True)

genes = pd.read_csv(f"{XCD}/genes.tsv", sep="\t", header=None, names=["gene_id", "chr"])
expected = set(genes.gene_id + "." + genes.chr + ".tsv")
have = set(os.listdir(f"{XCD}/per_gene"))
missing = sorted(expected - have)
print(f"{len(have):,} TSVs present of {len(expected):,} expected; {len(missing):,} missing (failed tasks)")
if missing:
    print("  e.g.", missing[:5])

t = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(f"{XCD}/per_gene/*.tsv"))],
              ignore_index=True)
pc = set(open(f"{XCD}/protein_coding_genes.txt").read().split())
n_all = t.gene_id.nunique()
t = t[t.gene_id.isin(pc)].copy()
print(f"protein-coding restriction: {t.gene_id.nunique():,} of {n_all:,} genes kept")
t.to_csv(f"{AGG}/crosscell_coloc_pairs.tsv", sep="\t", index=False)
print(f"{len(t):,} rows; note rows:")
print(t.note.fillna("coloc test").value_counts().to_string())

# denominator: every focal credible set
cs = t[t.note == "focal credible set"][KEY].drop_duplicates().copy()
tested = t[t.note.isna()]
best = tested.groupby(KEY).agg(n_other_tested=("other_cell", "nunique"), max_H4=("PP.H4.abf", "max")).reset_index()
part = (tested[tested["PP.H4.abf"] > H4].groupby(KEY).other_cell
        .apply(lambda s: ",".join(sorted(set(s)))).rename("shared_with").reset_index())
cs = cs.merge(best, on=KEY, how="left").merge(part, on=KEY, how="left")
cs["n_other_tested"] = cs.n_other_tested.fillna(0).astype(int)
cs["shared_with"] = cs.shared_with.fillna("")
cs["shared"] = (cs.max_H4 > H4).fillna(False).astype(int)
# credible-set size / top PIP / purity from the step 11 credible-set tables, for stratifying
meta = []
for f in glob.glob(f"{CMP}/*_cs.tsv"):
    m = pd.read_csv(f, sep="\t", usecols=["gene_id", "cell_type", "prior_label", "cs_idx", "cs_size", "top_pip", "min_abs_corr"])
    meta.append(m.rename(columns={"cell_type": "cell", "prior_label": "prior", "cs_idx": "idx"}))
if meta:
    cs = cs.merge(pd.concat(meta, ignore_index=True), on=KEY, how="left")
else:
    print(f"no credible-set tables in {CMP}: cs_size, top_pip and min_abs_corr not added")
cs.to_csv(f"{AGG}/crosscell_cs_status.tsv", sep="\t", index=False)

eg = (cs.groupby(["prior", "cell", "gene_id"])
        .agg(n_cs=("shared", "size"), n_cs_shared=("shared", "sum"),
             shared_with=("shared_with", lambda s: ",".join(sorted(set(",".join(s).split(",")) - {""}))))
        .reset_index())
eg["shared"] = (eg.n_cs_shared > 0).astype(int)
eg.to_csv(f"{AGG}/crosscell_egene_status.tsv", sep="\t", index=False)

fit = t[t.note.isin(["focal credible set", "no credible set"])].groupby(["prior", "cell"]).gene_id.nunique().rename("n_egenes_fit")
summ = (cs.groupby(["prior", "cell"])
          .agg(n_cs=("shared", "size"), n_cs_shared=("shared", "sum"),
               n_cs_with_partner=("n_other_tested", lambda s: int((s > 0).sum())))
          .join(fit).join(eg.groupby(["prior", "cell"]).agg(n_egenes_with_cs=("shared", "size"), n_egenes_shared=("shared", "sum")))
          .reset_index())
summ["pct_cs_shared"] = (100 * summ.n_cs_shared / summ.n_cs).round(1)
summ["pct_cs_shared_with_partner"] = (100 * summ.n_cs_shared / summ.n_cs_with_partner).round(1)
summ["pct_egenes_shared"] = (100 * summ.n_egenes_shared / summ.n_egenes_with_cs).round(1)
summ.to_csv(f"{AGG}/crosscell_summary.tsv", sep="\t", index=False)
print("\n" + summ.to_string(index=False))

# check: credible-set counts must equal the step 11 credible-set tables
if meta:
    # the step 11 tables cover all genes, so restrict them to protein-coding genes too
    ref = pd.concat(meta)
    ref = ref[ref.gene_id.isin(pc)].groupby(["prior", "cell"]).size().rename("n_cs_ref")
    chk = summ.set_index(["prior", "cell"]).n_cs.to_frame().join(ref)
    bad = chk[chk.n_cs != chk.n_cs_ref]
    print("\ncredible-set counts match the step 11 credible-set tables" if bad.empty else f"\nCOUNT MISMATCH:\n{bad}")

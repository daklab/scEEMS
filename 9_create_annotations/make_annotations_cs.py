#!/usr/bin/env python
"""
Credible-set fine-mapped eQTL annotations for one chromosome and cell type: 31 binary columns in one
annotation file, so LD scores are computed once and each column is then scored on its own
(ldscore_regression.py).

Fine-mapped eQTLs are defined by credible-set membership rather than by a bare PIP cut, so the definition
means the same under every fine-mapping prior: a prior that spreads posterior mass thinly puts many
variants above a PIP cut that are in no credible set.

The two column families cover different genes, deliberately:

  {cell}_cs                      1 column, eQTL + other genes: the fine-mapped comparator of the
                                 predicted-eQTL annotation (pred_prob > tau*), which covers the same genes.
                                 1 if the variant is a member of any gene's 95% credible set (FunGen-xQTL
                                 fine-mapping, step 1 exports PIP_top_parquet and PIP_top_other_parquet),
                                 with no PIP cut on top: a PIP cut would remove whole large credible sets
                                 (in a 50-variant set the mean PIP is 0.02) and keep only small, sharp ones.

  {cell}_{prior}_cs_pip{NN}      30 columns, five fine-mapping priors (step 10) x six PIP thresholds within
                                 the credible sets: PIP > 0 (pip00, the whole credible set), 0.10, 0.20,
                                 0.30, 0.40 and 0.50. Restricted to the eQTL genes fine-mapped under all
                                 five priors (gene_consensus.tsv, in_all_priors), so a gene missing under
                                 one prior cannot read as a heritability difference.

Variants are joined to the baseline on (BP, A1 = alternative allele, A2 = reference allele), exactly as the
predicted-eQTL annotations are, so both arms of a comparison are built the same way.

usage:   python make_annotations_cs.py CHR COHORT
output:  {sldsc_dir}/cs/{cohort}/MLxQTL_chr{CHR}.annot.gz and MLxQTL_chr{CHR}.l2.M
"""
import os
import sys

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import cohort_name, path

PRIORS = ["uniform", "EMS", "scEEMS_Weighted_Full", "scEEMS_Weighted_Restricted", "scEEMS_Unweighted_Full"]
THRESHOLDS = [0.10, 0.20, 0.30, 0.40, 0.50]

chrom = sys.argv[1]
cohort = cohort_name(sys.argv[2])
cell = cohort.replace("_mega_eQTL", "")
SUSIE = path("susie_pips_dir", cohort=cohort)
FM = f"{path('aggregate_dir')}/finemapping_comparison"
out_dir = f"{path('sldsc_dir')}/cs/{cohort}"
os.makedirs(out_dir, exist_ok=True)

baseline = pd.read_csv(f"{path('baseline_annot_dir')}/baseline_chr{chrom}.annot.gz", sep="\t",
                       usecols=["CHR", "SNP", "BP", "A1", "A2"])
annot = baseline.copy()


def add_column(name, keys):
    """keys: DataFrame with BP, A1, A2. Left-join onto the baseline, missing -> 0."""
    global annot
    if len(keys):
        k = keys.drop_duplicates().assign(**{name: 1})
        m = baseline.merge(k, on=["BP", "A1", "A2"], how="left")[["CHR", "SNP", "BP", "A1", "A2", name]]
    else:
        m = baseline.assign(**{name: np.nan})
    annot = annot.merge(m, on=["CHR", "SNP", "BP", "A1", "A2"], how="left")
    annot[name] = annot[name].fillna(0).astype(int)


# ---- column 1: eQTL + other genes, every member of a 95% credible set ----------------------------
# variant_id is chr:pos:REF:ALT, so the chromosome filter is a prefix test
cs_parts = []
for tag in ["PIP_top_parquet", "PIP_top_other_parquet"]:
    src = f"{SUSIE}/{tag}"
    if not os.path.isdir(src):
        continue
    t = ds.dataset(src, format="parquet").to_table(
        columns=["variant_id", "pos", "ref", "alt", "cs_coverage_0.95"]).to_pandas()
    t.columns = ["variant_id", "BP", "A2", "A1", "cs"]
    t = t[(t["cs"] > 0) & t["variant_id"].str.startswith(f"chr{chrom}:")]
    cs_parts.append(t[["BP", "A1", "A2"]])
union = pd.concat(cs_parts, ignore_index=True) if cs_parts else pd.DataFrame(columns=["BP", "A1", "A2"])
add_column(f"{cell}_cs", union)

# ---- columns 2-31: eQTL genes fine-mapped under all priors, credible-set member with PIP > t ------
consensus = pd.read_csv(f"{FM}/gene_consensus.tsv", sep="\t")
consensus = set(consensus.loc[consensus["in_all_priors"] & (consensus["cell_type"] == cell), "gene_id"])

for prior in PRIORS:
    d = pd.read_csv(f"{FM}/{cell}_{prior}_cs_variants.tsv", sep="\t")
    d = d[(d["chr"] == f"chr{chrom}") & d["gene_id"].isin(consensus)]
    parts = d["variant_id"].str.split(":", expand=True)
    d = d.assign(BP=parts[1].astype(np.int64), A2=parts[2], A1=parts[3])
    # PIP > 0: every credible-set member, named pip00 so it sorts with the rest of the sweep
    add_column(f"{cell}_{prior}_cs_pip00", d[["BP", "A1", "A2"]])
    for t in THRESHOLDS:
        add_column(f"{cell}_{prior}_cs_pip{int(round(t * 100)):02d}", d.loc[d["pip"] > t, ["BP", "A1", "A2"]])

annot.to_csv(f"{out_dir}/MLxQTL_chr{chrom}.annot.gz", sep="\t", index=False, compression="gzip")
sums = annot.drop(columns=["CHR", "SNP", "BP", "A1", "A2"]).to_numpy().sum(axis=0).reshape(1, -1)
np.savetxt(f"{out_dir}/MLxQTL_chr{chrom}.l2.M", sums, fmt="%f", delimiter=" ")
cols = [c for c in annot.columns if c not in ("CHR", "SNP", "BP", "A1", "A2")]
print(f"[{cohort} chr{chrom}] {len(annot):,} baseline variants, {len(cols)} columns", flush=True)
for c, s in zip(cols, sums[0]):
    print(f"    {c:<48} {int(s):>8,}", flush=True)

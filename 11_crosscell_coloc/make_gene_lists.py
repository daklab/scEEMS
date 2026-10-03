#!/usr/bin/env python
"""
The two gene lists of the cross-cell-type colocalization:

  genes.tsv                 every gene with an eQTL fine-mapping fit (step 10) in any cell type under any of
                            the four priors compared here (gene_id, chr; no header). One row per SLURM array
                            task of run_coloc_crosscell.sh.
  protein_coding_genes.txt  GENCODE protein-coding gene IDs (without version), from the GENCODE v45 basic
                            annotation GTF (gencode_gtf_file in config.yaml). aggregate_crosscell.py
                            restricts the comparison to these genes.

Prints the number of genes, which is the SLURM array size.

usage:   python make_gene_lists.py
outputs: {crosscell_dir}/genes.tsv, {crosscell_dir}/protein_coding_genes.txt
"""
import gzip
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

CELLS = ["Ast", "Exc", "Inh", "Mic", "Oli", "OPC"]
PRIORS = ["uniform", "scEEMS_Weighted_Full", "scEEMS_Unweighted_Full", "scEEMS_Weighted_Restricted"]
FIT_SUFFIX = ".univariate_bvsr.rds"

finemap_dir = path("finemap_dir")
out_dir = path("crosscell_dir")
os.makedirs(out_dir, exist_ok=True)

# genes with a fit: {gene_id}.{chr}.univariate_bvsr.rds in any fine_mapping_{prior}/ of any cell type
genes = set()
for cell in CELLS:
    for prior in PRIORS:
        d = f"{finemap_dir}/{cell}/fine_mapping_{prior}"
        if not os.path.isdir(d):
            print(f"  no fits for {cell} {prior} ({d})")
            continue
        for fn in os.listdir(d):
            if fn.endswith(FIT_SUFFIX):
                gene_id, chrom = fn[:-len(FIT_SUFFIX)].split(".")[:2]
                genes.add(f"{gene_id}\t{chrom}")
with open(f"{out_dir}/genes.tsv", "w") as fh:
    fh.writelines(f"{g}\n" for g in sorted(genes))

# protein-coding genes: "gene" features with gene_type "protein_coding"; the version suffix (and the _PAR_Y
# tag of pseudoautosomal copies) is dropped from gene_id
pattern = re.compile(r'gene_id "([^".]*)[^;]*; gene_type "protein_coding"')
coding = set()
with gzip.open(path("gencode_gtf_file"), "rt") as fh:
    for line in fh:
        fields = line.split("\t")
        if len(fields) > 8 and fields[2] == "gene":
            m = pattern.search(fields[8])
            if m:
                coding.add(m.group(1))
with open(f"{out_dir}/protein_coding_genes.txt", "w") as fh:
    fh.writelines(f"{g}\n" for g in sorted(coding))

print(f"{len(genes):,} fitted genes -> {out_dir}/genes.tsv (run_coloc_crosscell.sh: --array=1-{len(genes)}%200)")
print(f"{len(coding):,} protein-coding genes -> {out_dir}/protein_coding_genes.txt")

#!/usr/bin/env python
"""
Per-SNP functional priors (SNPVAR) for the PolyFun GWAS prior of precompute_gwas.R: PolyFun
(--compute-h2-L2) fits per-SNP heritability of the AD GWAS from 83 annotations, selected by S-LDSC:
  baseline_filtered (79 baseline-LF annotations)
  microglia_enhancer_promoter_union_atac_500                                     (brain_all)
  E051-H3K27ac.imputed.narrowPeak                                                (roadmap)
  chrombpnet_microglia-enhancer_promoter_union_atac_500-microglia_combined-zscore_2_GPN_pct_90
  chrombpnet_neuron-enhancer_promoter_union_atac_500-neuron_combined-zscore_2_GPN_pct_98
                                                                                 (chrombpnet_intersect_GPN)
The --anno order must match the LD-score columns, so each file's annotation names are read in file order
and filtered to the selection.

usage:   python compute_snpvar.py
inputs:  {ldsc_dir} (PolyFun installation), {sumstats_file} (LDSC-munged GWAS, parquet),
         {ldsc_annotation_dir}/{baseline_filtered/baseline,brain_all/brain_all,roadmap/roadmap,
         annotations_robust/chrombpnet_intersect_GPN/chrombpnet_intersect_GPN,weights/weights}_chr*
         annotation and LD-score files
output:  {snpvar_dir}/bellenguez_sldsc83.{1-22}.snpvar_ridge[_constrained].gz; then run snpvar_to_tabix.py
"""
import os
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

META = ["CHR", "SNP", "BP", "A1", "A2"]
A = path("ldsc_annotation_dir")
bl = f"{A}/baseline_filtered/baseline_chr"
brain = f"{A}/brain_all/brain_all_chr"
roadmap = f"{A}/roadmap/roadmap_chr"
gpn = f"{A}/annotations_robust/chrombpnet_intersect_GPN/chrombpnet_intersect_GPN_chr"
weight = f"{A}/weights/weights_chr"
out_dir = path("snpvar_dir")
os.makedirs(out_dir, exist_ok=True)

brain_selected = {"microglia_enhancer_promoter_union_atac_500"}
roadmap_selected = {"E051-H3K27ac.imputed.narrowPeak"}
gpn_selected = {"chrombpnet_microglia-enhancer_promoter_union_atac_500-microglia_combined-zscore_2_GPN_pct_90",
                "chrombpnet_neuron-enhancer_promoter_union_atac_500-neuron_combined-zscore_2_GPN_pct_98"}


def cols(prefix):
    return [c for c in pd.read_csv(f"{prefix}22.annot.gz", sep="\t", nrows=0).columns if c not in META]


bl_cols = cols(bl)
brain_cols = [c for c in cols(brain) if c in brain_selected]
road_cols = [c for c in cols(roadmap) if c in roadmap_selected]
gpn_cols = [c for c in cols(gpn) if c in gpn_selected]
for want, got, name in [(brain_selected, brain_cols, "brain_all"), (roadmap_selected, road_cols, "roadmap"),
                        (gpn_selected, gpn_cols, "chrombpnet_intersect_GPN")]:
    if len(got) != len(want):
        sys.exit(f"{name}: expected {len(want)} annotation(s), found {len(got)}")

anno = ",".join(bl_cols + brain_cols + road_cols + gpn_cols)
cmd = (f"python -u polyfun.py --compute-h2-L2 --no-partitions "
       f"--output-prefix {out_dir}/bellenguez_sldsc83 --sumstats {path('sumstats_file')} "
       f"--ref-ld-chr {bl},{brain},{roadmap},{gpn} --w-ld-chr {weight} "
       f"--allow-missing --nnls-exact --anno {anno} --q 100000")
print(f"{len(bl_cols)} baseline + {len(brain_cols) + len(road_cols) + len(gpn_cols)} selected annotations\n{cmd}\n",
      flush=True)
sys.exit(subprocess.call(cmd, shell=True, cwd=path("ldsc_dir")))

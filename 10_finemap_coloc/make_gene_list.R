#!/usr/bin/env Rscript
# List the genes whose GWAS window is fine-mapped: every autosomal gene in any cell type's region list
# that has a TADB cis window (the LD reference panel covers chromosomes 1-22). LD and GWAS fits depend only on the gene, so they are computed once per gene for all
# cell types. Prints the number of genes, the SLURM array size of run_precompute_ld.sh and
# run_precompute_gwas.sh.
#
# usage:   Rscript make_gene_list.R
# output:  {finemap_dir}/genes.tsv  (gene_id, chr; no header)

suppressMessages({ library(data.table) })

args_all   <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(file.path(script_dir, "common.R"))

regions <- rbindlist(lapply(CELLS, function(cell) {
  fread(eqtl_files(cell)$region_list, select = c("#chr", "ID"))
}))
setnames(regions, c("chr", "gene_id"))
tadb <- fread(tadb_file(), select = "gene_id")
genes <- unique(regions[gene_id %in% tadb$gene_id & chr %in% paste0("chr", 1:22), .(gene_id, chr)],
                by = "gene_id")
genes[, chr_int := as.integer(sub("chr", "", chr))]
setorder(genes, chr_int, gene_id)

out <- file.path(cfg_path("finemap_dir"), "genes.tsv")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
fwrite(genes[, .(gene_id, chr)], out, sep = "\t", col.names = FALSE)
cat(sprintf("%d genes -> %s\n", nrow(genes), out))

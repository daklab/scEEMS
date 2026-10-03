#!/usr/bin/env Rscript
# Colocalization posterior probability (CLPP; eCAVIAR, and the EMS paper) for every gene of one
# chromosome: for each cell type, eQTL prior and GWAS prior, the products PIP_eQTL(v) * PIP_GWAS(v) over
# the variants the two fits share. Unlike coloc.susie, which uses only the Bayes factors and so cannot see
# the EMS re-weighting of PIPs, CLPP multiplies PIPs and does respond to it. CLPP is not on the PP4 scale.
#   clpp          sum of the products (eCAVIAR)
#   max_product   largest product; its maximum over cell types is the EMS paper's gene-level CLPP
#   *_cs          the same restricted to variants in an eQTL credible set
#
# usage:   Rscript compute_clpp.R CHR
# outputs: {aggregate_dir}/clpp/clpp_bellenguez_chr{N}.tsv            one row per (GWAS prior, cell, prior, gene)
#          {aggregate_dir}/clpp/clppvar_bellenguez_chr{N}.tsv         variants with a product > 0.1
suppressMessages(library(data.table))

args_all   <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(file.path(script_dir, "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript compute_clpp.R CHR")
chrom <- args[1]
if (!grepl("^chr", chrom)) chrom <- paste0("chr", chrom)

FM    <- cfg_path("finemap_dir")
CACHE <- cfg_path("gwas_cache_dir")
OUT   <- file.path(cfg_path("aggregate_dir"), "clpp")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CUT <- 0.1

genes <- sub("\\.gwas\\.rds$", "", list.files(CACHE, pattern = paste0("\\.", chrom, "\\.gwas\\.rds$")))
have <- list()
for (cl in CELLS) for (p in PRIORS)
  have[[paste(cl, p)]] <- sub("\\.cs\\.tsv$", "", list.files(file.path(FM, cl, paste0("fine_mapping_", p)),
                                                             pattern = paste0("\\.", chrom, "\\.cs\\.tsv$")))
res <- list(); vres <- list()
for (g in genes) {
  o <- tryCatch(readRDS(file.path(CACHE, paste0(g, ".gwas.rds"))), error = function(e) NULL)
  if (is.null(o)) next
  for (arm in GWAS_PRIORS) {
    fit <- o[[paste0("gwas_", arm)]]
    if (is.null(fit) || is.null(fit$pip)) next
    nm <- names(fit$pip); if (is.null(nm)) nm <- o$joined$canon
    if (length(nm) != length(fit$pip)) next
    pg <- as.numeric(fit$pip); names(pg) <- nm
    for (cl in CELLS) for (p in PRIORS) {
      if (!(g %in% have[[paste(cl, p)]])) next
      cs <- tryCatch(fread(file.path(FM, cl, paste0("fine_mapping_", p), paste0(g, ".cs.tsv")),
                           showProgress = FALSE), error = function(e) NULL)
      if (is.null(cs) || !nrow(cs)) next
      shared <- intersect(cs$canon_key, nm)
      if (!length(shared)) next
      pe <- cs$pip[match(shared, cs$canon_key)]
      pv <- pe * pg[shared]
      in_cs <- !is.na(cs$cs_index[match(shared, cs$canon_key)]) & cs$cs_index[match(shared, cs$canon_key)] != ""
      res[[length(res) + 1]] <- data.table(
        gwas_id = GWAS_ID, gwas_prior = arm, cell_type = cl, eqtl_prior = p,
        gene_id = sub("\\..*$", "", g), chr = chrom,
        n_shared = length(shared), clpp = sum(pv), clpp_cs = sum(pv[in_cs]), max_product = max(pv),
        max_product_cs = if (any(in_cs)) max(pv[in_cs]) else NA_real_, n_cs_variants = sum(in_cs),
        top_variant = shared[which.max(pv)],
        eqtl_pip_at_top = pe[which.max(pv)], gwas_pip_at_top = pg[shared][which.max(pv)])
      hit <- which(pv > CUT)
      if (length(hit)) vres[[length(vres) + 1]] <- data.table(
        gwas_id = GWAS_ID, gwas_prior = arm, cell_type = cl, eqtl_prior = p,
        gene_id = sub("\\..*$", "", g), chr = chrom, variant = shared[hit],
        eqtl_pip = pe[hit], gwas_pip = pg[shared][hit], clpp_variant = pv[hit], in_cs = in_cs[hit])
    }
  }
}
fwrite(rbindlist(res), file.path(OUT, sprintf("clpp_%s_%s.tsv", GWAS_ID, chrom)), sep = "\t")
fwrite(rbindlist(vres), file.path(OUT, sprintf("clppvar_%s_%s.tsv", GWAS_ID, chrom)), sep = "\t")
cat(sprintf("[%s] %d genes with a GWAS fit, %d CLPP rows\n", chrom, length(genes), length(res)))

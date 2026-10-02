#!/usr/bin/env Rscript
# Collect the per-gene coloc.susie results of one cell type, eQTL prior and GWAS prior into one table,
# one row per tested (eQTL credible set, GWAS credible set) pair, joined to the eQTL credible-set metrics
# of aggregate_finemap.R on (gene, idx1 = cs_idx). Genes without a testable pair (a one-row note) are left
# out. hit1/hit2 are the variants coloc.susie used for each side, which need not be the credible set's
# top-PIP variant (cs_top_variant_id).
#
# usage:   Rscript aggregate_coloc.R CELL EQTL_PRIOR [GWAS_PRIOR=uniform]
# output:  {aggregate_dir}/finemapping_comparison/{CELL}_{EQTL_PRIOR}_bellenguez_{GWAS_PRIOR}_coloc.tsv
suppressMessages({ library(data.table) })

args_all   <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(file.path(script_dir, "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript aggregate_coloc.R CELL EQTL_PRIOR [GWAS_PRIOR]")
CELL       <- args[1]
EQTL_PRIOR <- args[2]
GWAS_PRIOR <- if (length(args) >= 3) args[3] else "uniform"
stopifnot(CELL %in% CELLS, EQTL_PRIOR %in% PRIORS, GWAS_PRIOR %in% GWAS_PRIORS)
OUT <- file.path(cfg_path("aggregate_dir"), "finemapping_comparison")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

in_dir <- file.path(cfg_path("finemap_dir"), CELL, sprintf("coloc_%s%s", EQTL_PRIOR, gwas_tag(GWAS_PRIOR)))
if (!dir.exists(in_dir)) stop("missing coloc dir: ", in_dir)
stem <- file.path(OUT, sprintf("%s_%s_%s_%s_coloc", CELL, EQTL_PRIOR, GWAS_ID, GWAS_PRIOR))

# result files start with "nsnps"; the no-pair notes start with "gene_id"
real <- system(sprintf("find -L %s -name '*.coloc.tsv' -print0 | xargs -0 -r grep -l '^nsnps' 2>/dev/null",
                       shQuote(in_dir)), intern = TRUE)
if (!length(real)) { cat("nothing to aggregate in", in_dir, "\n"); quit(save = "no", status = 0) }

dt <- rbindlist(lapply(real, fread), use.names = TRUE, fill = TRUE)
setnames(dt, gsub("^PP\\.(H[0-4])\\.abf$", "PP_\\1", names(dt)))
dt[, c("hit1_chr", "hit1_pos", "hit1_a1", "hit1_a2") := tstrsplit(hit1, ":", fixed = TRUE)]
dt[, c("hit2_chr", "hit2_pos", "hit2_a1", "hit2_a2") := tstrsplit(hit2, ":", fixed = TRUE)]
dt[, `:=`(hit1_chr = as.integer(sub("^chr", "", hit1_chr)), hit1_pos = as.integer(hit1_pos),
          hit2_chr = as.integer(sub("^chr", "", hit2_chr)), hit2_pos = as.integer(hit2_pos))]
stopifnot(!anyNA(dt$hit1_pos), !anyNA(dt$hit2_pos))
dt[, same_hit := hit1 == hit2]

cs_f <- file.path(OUT, sprintf("%s_%s_cs.tsv", CELL, EQTL_PRIOR))
if (file.exists(cs_f)) {
  cs <- fread(cs_f)[, .(gene_id, idx1 = cs_idx, cs_size, cs_rank, cs_coverage = coverage,
                        cs_min_abs_corr = min_abs_corr, cs_lbf = lbf, cs_top_pip = top_pip,
                        cs_top_variant_id = top_variant_id)]
  n_before <- nrow(dt)
  dt <- cs[dt, on = .(gene_id, idx1)]          # left join from the coloc pairs
  stopifnot(nrow(dt) == n_before)
  if (dt[is.na(cs_size), .N] > 0) warning("coloc pairs without a matching eQTL credible set", call. = FALSE)
} else {
  cat("NOTE: no", basename(cs_f), "-- run aggregate_finemap.R first to add the credible-set metrics\n")
}

dt[, `:=`(cell_type = CELL, eqtl_prior = EQTL_PRIOR, gwas_id = GWAS_ID, gwas_prior = GWAS_PRIOR)]
setcolorder(dt, c("cell_type", "eqtl_prior", "gwas_id", "gwas_prior", "gene_id", "chr", "idx1", "idx2"))
fwrite(dt, paste0(stem, ".tsv"), sep = "\t")
cat(sprintf("[%s/%s/%s] %d coloc pairs over %d genes, %d with PP.H4 > 0.8 -> %s.tsv\n", CELL, EQTL_PRIOR,
            GWAS_PRIOR, nrow(dt), uniqueN(dt$gene_id), dt[PP_H4 > 0.8, .N], stem))

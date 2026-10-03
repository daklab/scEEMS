#!/usr/bin/env Rscript
# Marginal (one variant at a time) eQTL association statistics for the colocalizing genes of one cell
# type, for the eQTL p-values of credset_all.tsv: SuSiE works on individual-level data and its saved fit
# keeps no summary statistics. The data are loaded exactly as in finemap_and_coloc.R (same window, QC and
# covariate residualization), so the variants line up with the credible sets. The statistics do not
# depend on the prior, so one file per (cell type, gene) serves all five priors. The standard errors use
# n - 2 residual degrees of freedom on covariate-residualized data, so the p-values are slightly
# anti-conservative; they annotate figures and tables and are not used for inference.
#
# usage:   Rscript compute_eqtl_marginal.R CELL       (genes: the eQTL side of credset_all.tsv)
# output:  {finemap_dir}/eqtl_marginal/{CELL}/{gene}.{chr}.marginal.tsv
#          variant_id, canon_key, betahat, sebetahat, z, neglog10p
suppressMessages({ library(pecotmr); library(susieR); library(data.table); library(dplyr); library(readr) })

args_all   <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(file.path(script_dir, "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript compute_eqtl_marginal.R CELL")
CELL <- args[1]
stopifnot(CELL %in% CELLS)
if (!nzchar(Sys.which("tabix"))) stop("tabix is not on PATH (pecotmr needs it to read the phenotypes)")

targets <- file.path(cfg_path("aggregate_dir"), "finemapping_comparison", "credset_all.tsv")
tg <- unique(fread(targets)[side == "eqtl" & cell_type == CELL, .(gene_id, chr)])
if (!nrow(tg)) { cat("no colocalizing genes for", CELL, "\n"); quit(save = "no", status = 0) }
out_dir <- file.path(cfg_path("finemap_dir"), "eqtl_marginal", CELL)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files     <- eqtl_files(CELL)
region_df <- read_delim(files$region_list, delim = "\t", escape_double = FALSE, trim_ws = TRUE,
                        show_col_types = FALSE) %>% rename(chr = `#chr`)
tadb      <- read_delim(files$tadb, delim = "\t", escape_double = FALSE, trim_ws = TRUE,
                        show_col_types = FALSE) %>% rename(chr = `#chr`)

n_ok <- 0; n_skip <- 0; n_fail <- 0
for (i in seq_len(nrow(tg))) {
  gid <- tg$gene_id[i]; chr_str <- tg$chr[i]
  out_f <- file.path(out_dir, sprintf("%s.%s.marginal.tsv", gid, chr_str))
  if (file.exists(out_f)) { n_skip <- n_skip + 1; next }
  region   <- region_df[region_df$ID == gid, ]
  tadb_row <- tadb %>% filter(gene_id == !!gid)
  if (nrow(region) == 0 || nrow(tadb_row) == 0) { n_fail <- n_fail + 1; next }

  fdat <- tryCatch(load_eqtl_data(CELL, chr_str,
                                  sprintf("%s:%d-%d", chr_str, as.integer(region$start[1]), as.integer(region$end[1])),
                                  sprintf("%s:%d-%d", chr_str, as.integer(tadb_row$start[1]), as.integer(tadb_row$end[1]))),
                   error = function(e) { message("  load failed ", gid, ": ", conditionMessage(e)); NULL })
  if (is.null(fdat)) { n_fail <- n_fail + 1; next }
  X  <- fdat$residual_X[[1]]
  Ym <- tryCatch(single_phenotype(fdat$residual_Y[[1]], gid), error = function(e) NULL)
  if (is.null(Ym)) { n_fail <- n_fail + 1; next }

  ss <- univariate_regression(X, as.vector(Ym))
  z  <- ss$betahat / ss$sebetahat
  fwrite(data.table(variant_id = colnames(X), canon_key = pred_id_to_canon(colnames(X)),
                    betahat = ss$betahat, sebetahat = ss$sebetahat, z = z,
                    neglog10p = -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)),
         out_f, sep = "\t")
  n_ok <- n_ok + 1
}
cat(sprintf("[%s] %d genes: %d computed, %d already present, %d failed\n", CELL, nrow(tg), n_ok, n_skip, n_fail))

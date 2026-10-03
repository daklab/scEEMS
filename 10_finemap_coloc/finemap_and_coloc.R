#!/usr/bin/env Rscript
# eQTL fine-mapping of one gene in one cell type under one prior, then colocalization with the AD GWAS.
#
#   1. Load the gene's eQTL data (genotypes, pseudobulk expression, covariates) over its TADB cis window.
#   2. Fit SuSiE with the chosen eQTL prior:
#        uniform                  susie(X, y)
#        scEEMS_Weighted_Full     susie(X, y, prior_weights = recalibrated scEEMS predictions); likewise
#        scEEMS_Unweighted_Full   for the other two models of step 5. Predictions are recalibrated from
#        scEEMS_Weighted_Restricted  the balanced training prior (0.5) to 1.8 / 14,716 per variant by a
#                                 logit shift; variants without a prediction get 1.8 / 14,716.
#        EMS                      the uniform fit, with alpha re-weighted within each pure credible set
#                                 by the gene's EMS scores (prepare_ems_priors.py), as in the EMS paper.
#   3. Save the fit and a per-variant credible-set table. The eQTL fit does not depend on the GWAS, so it
#      is computed once per (cell, prior, gene) and reused by both GWAS priors.
#   4. Load the gene's GWAS fit (precompute_gwas.R; uniform or PolyFun prior) and run coloc.susie on the
#      variants the two fits share, matched by canonical key (chr:pos:alleles sorted).
#
# usage:  Rscript finemap_and_coloc.R IDX CELL PRIOR [GWAS_PRIOR]
#   IDX         1-based row of the cell type's region list (one gene)
#   CELL        Mic | Ast | Exc | Inh | Oli | OPC
#   PRIOR       uniform | EMS | scEEMS_Weighted_Full | scEEMS_Unweighted_Full | scEEMS_Weighted_Restricted
#   GWAS_PRIOR  uniform (default) | polyfun
# outputs (finemap_dir = config setting finemap_dir):
#   {finemap_dir}/{CELL}/fine_mapping_{PRIOR}/{gene}.{chr}.univariate_bvsr.rds and .cs.tsv
#   {finemap_dir}/{CELL}/gwas_susie_{PRIOR}{tag}/{gene}.{chr}.gwas.rds
#   {finemap_dir}/{CELL}/coloc_{PRIOR}{tag}/{gene}.{chr}.coloc.tsv     tag "" (uniform) or ".polyfun"
# A task whose coloc table exists exits immediately, so resubmitting an array only reruns missing genes.

suppressMessages({
  library(pecotmr)
  library(susieR)
  library(coloc)
  library(data.table)
  library(dplyr)
  library(readr)
  library(Matrix)
})

args_all   <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(file.path(script_dir, "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: Rscript finemap_and_coloc.R IDX CELL PRIOR [GWAS_PRIOR]")
idx         <- as.integer(args[1])
cell_type   <- args[2]
prior_label <- args[3]
gwas_prior  <- if (length(args) >= 4) args[4] else "uniform"
stopifnot(cell_type %in% CELLS, prior_label %in% PRIORS, gwas_prior %in% GWAS_PRIORS)

FINEMAP_DIR    <- cfg_path("finemap_dir")
LD_CACHE_DIR   <- cfg_path("ld_cache_dir")
GWAS_CACHE_DIR <- cfg_path("gwas_cache_dir")
files          <- eqtl_files(cell_type)

# ----- helpers -----

# Logit-shift recalibration of pred_prob from the balanced training prior (0.5) to prior_probability.
recalibrate_predictions <- function(pred_path, eqtl_canon, prior_probability = 1.8 / 14716) {
  if (!file.exists(pred_path)) return(NULL)
  pred <- fread(pred_path, select = c("variant_id", "pred_prob"))
  logit_intercept <- log((prior_probability / (1 - prior_probability)) / (0.5 / 0.5))
  pred[, logit_p := log(pred_prob / (1 - pred_prob))]
  pred[, recal   := 1 / (1 + exp(-(logit_p + logit_intercept)))]
  pred[, canon   := pred_id_to_canon(variant_id)]
  lk <- setNames(pred$recal, pred$canon)
  vec <- unname(lk[eqtl_canon])
  vec[is.na(vec)] <- prior_probability          # variants without a prediction
  vec
}

# EMS: re-weight alpha within each pure credible set (min |r| >= 0.5) by the variant's EMS score f, keeping
# the effect's total probability mass; variants without a score get p_random. lbf_variable, V, sets and
# lbf are left as fitted, so credible sets are those of the uniform fit.
ems_alpha_reweight <- function(susie_fit, eqtl_canon, cell_type, gene_id) {
  ems_dir      <- file.path(FINEMAP_DIR, cell_type, "ems_vector")
  ems_tsv      <- file.path(ems_dir, sprintf("%s.tsv", gene_id))
  p_random_txt <- file.path(ems_dir, "p_random.txt")
  if (!file.exists(p_random_txt)) stop(sprintf("EMS p_random.txt missing for cell %s", cell_type))
  p_random <- as.numeric(readLines(p_random_txt)[1])

  if (file.exists(ems_tsv)) {
    ems_df <- fread(ems_tsv)
    ems_df[, canon := pred_id_to_canon(variant_id)]
    ems_lk <- setNames(ems_df$ems, ems_df$canon)
    f <- unname(ems_lk[eqtl_canon])
    f[is.na(f)] <- p_random
  } else {
    f <- rep(p_random, length(eqtl_canon))
  }

  L <- nrow(susie_fit$alpha)
  is_pure_l <- rep(FALSE, L)
  if (!is.null(susie_fit$sets$cs_index) && length(susie_fit$sets$cs_index) > 0) {
    pure_mask <- susie_fit$sets$purity$min.abs.corr >= 0.5
    is_pure_l[susie_fit$sets$cs_index[pure_mask]] <- TRUE
  }
  alpha_new <- susie_fit$alpha
  for (l in which(is_pure_l)) {
    w <- susie_fit$alpha[l, ] * f
    row_sum <- sum(susie_fit$alpha[l, ])
    if (sum(w) > 0) alpha_new[l, ] <- (w / sum(w)) * row_sum
  }
  susie_fit$alpha <- alpha_new
  susie_fit$pip   <- 1 - apply(1 - alpha_new, 2, prod)
  susie_fit
}

# Per-variant table: PIP and the effect index of the credible set the variant belongs to (NA if none).
cs_table <- function(susie_fit, eqtl_canon, eqtl_orig_id) {
  cs_id <- rep(NA_integer_, length(eqtl_canon))
  if (!is.null(susie_fit$sets$cs) && length(susie_fit$sets$cs) > 0) {
    cs_index <- susie_fit$sets$cs_index
    for (i in seq_along(susie_fit$sets$cs)) cs_id[susie_fit$sets$cs[[i]]] <- cs_index[i]
  }
  data.table(variant_id = eqtl_orig_id, canon_key = eqtl_canon, pip = susie_fit$pip, cs_index = cs_id)
}

# Restrict a susie fit to a subset of variants for coloc.susie, re-indexing the credible sets.
subset_susie_variants <- function(fit, keep_idx) {
  out <- list()
  out$alpha        <- fit$alpha[, keep_idx, drop = FALSE]
  out$lbf_variable <- fit$lbf_variable[, keep_idx, drop = FALSE]
  out$pip          <- fit$pip[keep_idx]
  out$V            <- fit$V
  out$lbf          <- fit$lbf
  out$sigma2       <- if (!is.null(fit$sigma2)) fit$sigma2 else 1
  if (!is.null(fit$sets) && length(fit$sets$cs) > 0) {
    pos_lookup <- match(seq_along(fit$pip), keep_idx)      # old index -> new index (NA if dropped)
    new_cs <- list(); new_idx <- integer(0)
    for (i in seq_along(fit$sets$cs)) {
      mapped <- pos_lookup[fit$sets$cs[[i]]]
      mapped <- mapped[!is.na(mapped)]
      if (length(mapped) > 0) {
        new_cs[[length(new_cs) + 1]] <- mapped
        new_idx <- c(new_idx, fit$sets$cs_index[i])
      }
    }
    out$sets <- list(cs = new_cs, cs_index = new_idx,
                     coverage = fit$sets$coverage[match(new_idx, fit$sets$cs_index)],
                     purity   = fit$sets$purity[match(new_idx, fit$sets$cs_index), , drop = FALSE])
  } else {
    out$sets <- list(cs = list(), cs_index = integer(0))
  }
  class(out) <- "susie"
  out
}

# ----- 1. region -----
region_df <- read_delim(files$region_list, delim = "\t", escape_double = FALSE, trim_ws = TRUE,
                        show_col_types = FALSE) %>% rename(chr = `#chr`)
region  <- region_df[idx, ]
gene_id <- region$ID
chr_str <- region$chr
chr_int <- as.integer(sub("chr", "", chr_str))

tadb <- read_delim(files$tadb, delim = "\t", escape_double = FALSE, trim_ws = TRUE,
                   show_col_types = FALSE) %>% rename(chr = `#chr`)
tadb_row <- tadb %>% filter(gene_id == !!gene_id)
if (nrow(tadb_row) == 0) {
  message("gene_id ", gene_id, " not in TADB; skipping")
  quit(save = "no", status = 0)
}
start_bp <- as.integer(tadb_row$start[1])
end_bp   <- as.integer(tadb_row$end[1])
association_id <- sprintf("%s:%d-%d", chr_str, start_bp, end_bp)
region_id      <- sprintf("%s:%d-%d", chr_str, as.integer(region$start), as.integer(region$end))
message(sprintf("[%s] %s %s  window=%s  GWAS prior=%s", cell_type, prior_label, gene_id,
                association_id, gwas_prior))

# ----- 2. output paths -----
tag       <- gwas_tag(gwas_prior)
fm_dir    <- file.path(FINEMAP_DIR, cell_type, sprintf("fine_mapping_%s", prior_label))
gwas_dir  <- file.path(FINEMAP_DIR, cell_type, sprintf("gwas_susie_%s%s", prior_label, tag))
coloc_dir <- file.path(FINEMAP_DIR, cell_type, sprintf("coloc_%s%s", prior_label, tag))
for (d in c(fm_dir, gwas_dir, coloc_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
eqtl_rds_path  <- file.path(fm_dir,    sprintf("%s.%s.univariate_bvsr.rds", gene_id, chr_str))
cs_tsv_path    <- file.path(fm_dir,    sprintf("%s.%s.cs.tsv", gene_id, chr_str))
gwas_rds_path  <- file.path(gwas_dir,  sprintf("%s.%s.gwas.rds", gene_id, chr_str))
coloc_tsv_path <- file.path(coloc_dir, sprintf("%s.%s.coloc.tsv", gene_id, chr_str))

if (file.exists(coloc_tsv_path)) {
  message("coloc table exists, skipping: ", coloc_tsv_path)
  quit(save = "no", status = 0)
}

# ----- 3. eQTL fit: reuse a saved fit (from the other GWAS prior), else compute -----
reuse_eqtl <- file.exists(eqtl_rds_path)
if (reuse_eqtl) {
  e_src <- tryCatch(readRDS(eqtl_rds_path), error = function(e) NULL)
  if (is.null(e_src)) {
    reuse_eqtl <- FALSE                       # unreadable (partly written) fit: recompute
  } else {
    message("reusing eQTL fit ", eqtl_rds_path)
    eqtl_fit     <- e_src$susie_fitted
    eqtl_canon   <- e_src$canon_keys
    eqtl_orig_id <- e_src$variant_names
  }
}

if (!reuse_eqtl) {
  fdat <- tryCatch(load_eqtl_data(cell_type, chr_str, region_id, association_id),
                   error = function(e) { message("eQTL load failed: ", conditionMessage(e)); NULL })
  if (is.null(fdat)) quit(save = "no", status = 0)

  X_eqtl       <- fdat$residual_X[[1]]
  Y_eqtl       <- as.vector(single_phenotype(fdat$residual_Y[[1]], gene_id))
  eqtl_orig_id <- colnames(X_eqtl)
  eqtl_canon   <- pred_id_to_canon(eqtl_orig_id)
  message(sprintf("eQTL data: %d samples x %d variants", nrow(X_eqtl), ncol(X_eqtl)))

  # prior weights (NULL = uniform); susie() normalizes them to sum to 1
  prior_probability <- 1.8 / 14716
  prior_weights <- if (prior_label %in% c("uniform", "EMS")) {
    NULL
  } else {
    pred_path <- file.path(cfg_path("predictions_dir", cohort = cohort_name(cell_type)),
                           SCEEMS_MODEL[[prior_label]], sprintf("%s_predictions.tsv", gene_id))
    pw <- recalibrate_predictions(pred_path, eqtl_canon, prior_probability)
    if (is.null(pw)) {
      message("no predictions for ", gene_id, " (", pred_path, "); skipping")
      quit(save = "no", status = 0)
    }
    pw
  }

  # SuSiE with adaptive L: start at 5 and add 5 (up to 30) while every effect forms a credible set
  set.seed(123)
  susie_wrapper <- function(X, y, prior_weights, init_L = 5, max_L = 30, l_step = 5, ...) {
    L <- init_L
    repeat {
      res <- susie(X, y, L = L, prior_weights = prior_weights, ...)
      if (!is.null(res$sets$cs) && length(res$sets$cs) >= L && L < max_L) {
        L <- L + l_step
      } else break
    }
    res
  }
  eqtl_fit <- tryCatch(
    susie_wrapper(X_eqtl, Y_eqtl, prior_weights = prior_weights,
                  scaled_prior_variance = 0.2, coverage = 0.95, min_abs_corr = 0.5,
                  estimate_residual_variance = TRUE, estimate_prior_variance = TRUE,
                  max_iter = 1000, verbose = FALSE),
    error = function(e) { message("eQTL susie failed: ", conditionMessage(e)); NULL })
  if (is.null(eqtl_fit)) quit(save = "no", status = 0)

  if (prior_label == "EMS") eqtl_fit <- ems_alpha_reweight(eqtl_fit, eqtl_canon, cell_type, gene_id)

  saveRDS(list(susie_fitted = eqtl_fit,
               variant_names = eqtl_orig_id,
               canon_keys    = eqtl_canon,
               region_info   = list(gene_id = gene_id, chr = chr_str, start = start_bp, end = end_bp,
                                    n_variants = length(eqtl_orig_id), cell_type = cell_type,
                                    prior_label = prior_label)),
          eqtl_rds_path, compress = "xz")
  fwrite(cs_table(eqtl_fit, eqtl_canon, eqtl_orig_id), cs_tsv_path, sep = "\t")
  message("eQTL fit saved: ", eqtl_rds_path)
}

# ----- 4. GWAS fit (precompute_gwas.R) -----
gwas_cache_path <- file.path(GWAS_CACHE_DIR, sprintf("%s.%s.gwas.rds", gene_id, chr_str))
ld_cache_path   <- file.path(LD_CACHE_DIR,   sprintf("%s.%s.ld.rds",  gene_id, chr_str))
if (!file.exists(gwas_cache_path)) {
  message("no GWAS fit for ", gene_id, ": ", gwas_cache_path, " (run precompute_gwas.R first)")
  quit(save = "no", status = 0)
}
gwas_cache <- readRDS(gwas_cache_path)
# the variants the GWAS fit indexes (saved with the fit by precompute_gwas.R, else read from the LD cache)
joined <- if (!is.null(gwas_cache$joined)) gwas_cache$joined else {
  if (!file.exists(ld_cache_path)) {
    message("LD cache missing for ", gene_id, ": ", ld_cache_path)
    quit(save = "no", status = 0)
  }
  readRDS(ld_cache_path)$joined
}
if (is.null(joined) || nrow(joined) < 10) {
  message("GWAS fit has too few variants; skipping coloc")
  quit(save = "no", status = 0)
}
gwas_fit <- gwas_cache[[paste0("gwas_", gwas_prior)]]
if (is.null(gwas_fit)) {
  message("no ", gwas_prior, " GWAS fit for ", gene_id, "; skipping")
  quit(save = "no", status = 0)
}
n_for_susie <- gwas_cache$n_for_susie

colnames(eqtl_fit$alpha)        <- eqtl_canon
colnames(eqtl_fit$lbf_variable) <- eqtl_canon
names(eqtl_fit$pip)             <- eqtl_canon

saveRDS(list(susie_fitted = gwas_fit,
             variant_names = joined$canon,
             region_info   = list(gene_id = gene_id, chr = chr_str, start = start_bp, end = end_bp,
                                  n_variants = nrow(joined), n_for_susie = n_for_susie,
                                  cell_type = cell_type, prior_label = prior_label)),
        gwas_rds_path, compress = "xz")

# ----- 5. coloc.susie on the shared variants -----
shared <- intersect(eqtl_canon, joined$canon)
if (length(shared) < 10) {
  message("too few shared variants for coloc; skipping")
  quit(save = "no", status = 0)
}
eqtl_for_coloc <- subset_susie_variants(eqtl_fit, match(shared, eqtl_canon))
gwas_for_coloc <- subset_susie_variants(gwas_fit, match(shared, joined$canon))
for (nm in c("alpha", "lbf_variable")) {
  colnames(eqtl_for_coloc[[nm]]) <- shared
  colnames(gwas_for_coloc[[nm]]) <- shared
}
names(eqtl_for_coloc$pip) <- shared
names(gwas_for_coloc$pip) <- shared

coloc_res <- tryCatch(coloc::coloc.susie(eqtl_for_coloc, gwas_for_coloc),
                      error = function(e) { message("coloc.susie failed: ", conditionMessage(e)); NULL })
if (is.null(coloc_res) || is.null(coloc_res$summary)) {
  # a gene without a testable pair of credible sets gets a one-row note instead of results
  fwrite(data.table(gene_id = gene_id, chr = chr_str, cell_type = cell_type, prior_label = prior_label,
                    n_shared = length(shared), note = "coloc returned NULL summary"),
         coloc_tsv_path, sep = "\t")
  quit(save = "no", status = 0)
}
coloc_dt <- as.data.table(coloc_res$summary)
coloc_dt[, `:=`(gene_id = gene_id, chr = chr_str, cell_type = cell_type, prior_label = prior_label,
                n_shared = length(shared), n_for_susie = n_for_susie)]
fwrite(coloc_dt, coloc_tsv_path, sep = "\t")
message("coloc saved: ", coloc_tsv_path)

#!/usr/bin/env Rscript
# Compare the fine-mapping priors against the uniform prior, from the aggregate_finemap.R tables of all
# cell types and priors.
#   gene_level_comparison.tsv   one row per (cell type, gene, prior): credible-set count, purity, size,
#                               log Bayes factor and top PIP summarized per gene; genes without a
#                               credible set are kept.
#   gene_consensus.tsv          which priors fine-mapped each (cell type, gene). The scEEMS priors need
#                               predictions for the gene, so they cover fewer genes than uniform and EMS;
#                               compare priors on genes with in_all_priors = TRUE.
#   cs_matched_comparison.tsv   credible sets matched between uniform and each other prior on (gene,
#                               effect index), with the Jaccard overlap of their variants and whether they
#                               share the top variant. MATCH_MODE=overlap pairs credible sets by
#                               mutual-best Jaccard overlap instead (a sensitivity analysis).
#   cs_match_summary.tsv        per (cell type, prior): matched, uniform-only and prior-only credible sets.
#
# usage:   Rscript build_comparison.R
# output:  {aggregate_dir}/finemapping_comparison/
suppressMessages({ library(data.table) })

args_all   <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(file.path(script_dir, "common.R"))

DIR <- file.path(cfg_path("aggregate_dir"), "finemapping_comparison")
REF <- "uniform"
COMPARATORS <- setdiff(PRIORS, REF)
MATCH_MODE <- Sys.getenv("MATCH_MODE", "cs_idx")
stopifnot(MATCH_MODE %in% c("cs_idx", "overlap"))

rd <- function(cell, prior, kind) {
  f <- file.path(DIR, sprintf("%s_%s_%s.tsv", cell, prior, kind))
  if (!file.exists(f)) { message("  missing: ", basename(f)); return(NULL) }
  fread(f)
}

# ---------------------------------------------------------------- gene level
gene_rows <- list()
for (cell in CELLS) for (pr in PRIORS) {
  g  <- rd(cell, pr, "gene"); if (is.null(g)) next
  cs <- rd(cell, pr, "cs")
  per_gene <- if (is.null(cs) || !nrow(cs)) NULL else cs[, .(
    mean_purity = mean(min_abs_corr, na.rm = TRUE),
    min_purity  = min(min_abs_corr,  na.rm = TRUE),
    mean_cs_size = mean(cs_size), median_cs_size = as.numeric(median(cs_size)),
    mean_lbf = mean(lbf, na.rm = TRUE), max_lbf = max(lbf, na.rm = TRUE),
    mean_top_pip = mean(top_pip), best_top_pip = max(top_pip)), by = gene_id]
  if (!is.null(per_gene)) g <- per_gene[g, on = "gene_id"]
  gene_rows[[paste(cell, pr)]] <- g
}
gene_dt <- rbindlist(gene_rows, fill = TRUE)
setcolorder(gene_dt, c("cell_type", "prior_label", "gene_id", "chr", "n_cs"))
fwrite(gene_dt, file.path(DIR, "gene_level_comparison.tsv"), sep = "\t")

# ---------------------------------------------------------------- gene coverage per prior
present <- unique(gene_dt[, .(cell_type, gene_id, prior_label)])
present[, is_present := TRUE]
prior_cols <- intersect(PRIORS, unique(present$prior_label))
if (length(prior_cols) < length(PRIORS))
  warning("no fits for prior(s) ", paste(setdiff(PRIORS, prior_cols), collapse = ", "),
          "; in_all_priors does not cover them", call. = FALSE)
cov_dt <- dcast(present, cell_type + gene_id ~ prior_label, value.var = "is_present", fill = FALSE)
for (p in prior_cols) set(cov_dt, j = p, value = as.logical(cov_dt[[p]]))
cov_dt[, n_priors      := rowSums(as.matrix(.SD)), .SDcols = prior_cols]
cov_dt[, in_all_priors := n_priors == length(prior_cols)]
setcolorder(cov_dt, c("cell_type", "gene_id", "n_priors", "in_all_priors"))
fwrite(cov_dt, file.path(DIR, "gene_consensus.tsv"), sep = "\t")

# ---------------------------------------------------------------- credible-set matching
REF_COLS  <- function(dt) dt[, .(gene_id, chr, ref_cs = cs_idx, ref_size = cs_size,
                                 ref_purity = min_abs_corr, ref_lbf = lbf, ref_top_pip = top_pip,
                                 ref_coverage = coverage, ref_top_variant = top_variant_id)]
COMP_COLS <- function(dt) dt[, .(gene_id, comp_cs = cs_idx, comp_size = cs_size,
                                 comp_purity = min_abs_corr, comp_lbf = lbf, comp_top_pip = top_pip,
                                 comp_coverage = coverage, comp_top_variant = top_variant_id)]

match_cs <- function(cell, comp) {
  rv <- rd(cell, REF,  "cs_variants"); cv <- rd(cell, comp, "cs_variants")
  rc <- rd(cell, REF,  "cs");          cc <- rd(cell, comp, "cs")
  if (is.null(rv) || is.null(cv) || is.null(rc) || is.null(cc)) return(NULL)

  if (MATCH_MODE == "cs_idx") {
    out <- merge(REF_COLS(rc), COMP_COLS(cc), by = c("gene_id"), allow.cartesian = TRUE)
    out <- out[ref_cs == comp_cs]
    inter <- merge(rv[, .(gene_id, cs_idx, canon_key)], cv[, .(gene_id, cs_idx, canon_key)],
                   by = c("gene_id", "cs_idx", "canon_key"))[, .(n_inter = .N), by = .(gene_id, cs_idx)]
    out[inter, n_inter := i.n_inter, on = c("gene_id", "ref_cs" = "cs_idx")]
    out[is.na(n_inter), n_inter := 0L]
  } else {
    inter <- merge(rv[, .(gene_id, canon_key, ref_cs = cs_idx)],
                   cv[, .(gene_id, canon_key, comp_cs = cs_idx)],
                   by = c("gene_id", "canon_key"), allow.cartesian = TRUE)
    if (!nrow(inter)) return(list(matched = NULL, ref = rc, comp = cc))
    pair <- inter[, .(n_inter = .N), by = .(gene_id, ref_cs, comp_cs)]
    pair[rc[, .(gene_id, ref_cs = cs_idx, s = cs_size)],  ref_size  := i.s, on = c("gene_id", "ref_cs")]
    pair[cc[, .(gene_id, comp_cs = cs_idx, s = cs_size)], comp_size := i.s, on = c("gene_id", "comp_cs")]
    pair[, j := n_inter / (ref_size + comp_size - n_inter)]
    pair[, best_for_ref  := j == max(j), by = .(gene_id, ref_cs)]
    pair[, best_for_comp := j == max(j), by = .(gene_id, comp_cs)]
    m <- pair[best_for_ref & best_for_comp & j > 0]
    setorder(m, gene_id, ref_cs, -j, comp_cs)
    m <- m[, .SD[1], by = .(gene_id, ref_cs)][, .SD[1], by = .(gene_id, comp_cs)]
    out <- merge(m[, .(gene_id, ref_cs, comp_cs, n_inter)], REF_COLS(rc), by = c("gene_id", "ref_cs"))
    out <- merge(out, COMP_COLS(cc), by = c("gene_id", "comp_cs"))
  }
  if (!nrow(out)) return(list(matched = NULL, ref = rc, comp = cc))
  out[, jaccard := n_inter / (ref_size + comp_size - n_inter)]
  out[, `:=`(cell_type = cell, ref_prior = REF, comp_prior = comp, match_mode = MATCH_MODE,
             same_top_variant = ref_top_variant == comp_top_variant)]
  list(matched = out, ref = rc, comp = cc)
}

matched_all <- list(); summary_rows <- list()
for (cell in CELLS) for (comp in COMPARATORS) {
  r <- match_cs(cell, comp); if (is.null(r)) next
  m <- r$matched
  n_ref <- nrow(r$ref); n_comp <- nrow(r$comp); n_m <- if (is.null(m)) 0L else nrow(m)
  matched_all[[paste(cell, comp)]] <- m
  summary_rows[[paste(cell, comp)]] <- data.table(
    cell_type = cell, ref_prior = REF, comp_prior = comp,
    n_cs_ref = n_ref, n_cs_comp = n_comp, n_matched = n_m,
    ref_only = n_ref - n_m, comp_only = n_comp - n_m,
    pct_ref_matched = 100 * n_m / max(1, n_ref),
    median_jaccard = if (n_m) median(m$jaccard) else NA_real_,
    pct_same_top_variant = if (n_m) 100 * mean(m$same_top_variant) else NA_real_)
}
fwrite(rbindlist(Filter(Negate(is.null), matched_all), fill = TRUE),
       file.path(DIR, "cs_matched_comparison.tsv"), sep = "\t")
fwrite(rbindlist(summary_rows), file.path(DIR, "cs_match_summary.tsv"), sep = "\t")
cat(sprintf("wrote gene_level_comparison.tsv (%d rows), gene_consensus.tsv (%d (cell, gene) pairs, %d in all priors), cs_matched_comparison.tsv, cs_match_summary.tsv\n",
            nrow(gene_dt), nrow(cov_dt), sum(cov_dt$in_all_priors)))

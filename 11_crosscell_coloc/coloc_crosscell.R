#!/usr/bin/env Rscript
# Cross-cell-type colocalization of eQTL credible sets for ONE gene, all priors.
#
# For every prior, each cell type with at least one credible set for the gene is compared, with
# coloc.susie on the stored SuSiE fits of step 10 (no refitting), against every other cell type that also
# has a credible set for the gene under the same prior. Variants are aligned on the canonical keys stored
# with each fit, as in the eQTL-GWAS colocalization of step 10. One row per (prior, cell, idx, other_cell,
# other_idx) with PP.H0-PP.H4. Every case where nothing can be tested is written as a note row instead of
# an error, so the aggregation can count it:
#   "focal credible set"                one row per credible set of (prior, cell): the denominator
#   "no fit" / "no credible set"        this cell type has no fit / no credible set under this prior
#   "... in other cell"                 the other cell type has no fit / no credible set
#   "fewer than 10 shared variants", "credible set has no shared variants", "coloc returned no result"
# A TSV is always written, so a missing TSV means a failed task.
#
# EMS is not compared: it re-weights PIPs inside the uniform fit but leaves its credible sets and Bayes
# factors unchanged, so its results would be identical to uniform's.
#
# usage:   Rscript coloc_crosscell.R GENE_INDEX       (row of {crosscell_dir}/genes.tsv: gene_id, chr)
# output:  {crosscell_dir}/per_gene/{gene_id}.{chr}.tsv
suppressMessages({library(data.table); library(coloc)})
args_all   <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(file.path(script_dir, "..", "shared", "config.R"))

idx    <- as.integer(commandArgs(trailingOnly = TRUE)[1])
CELLS  <- c("Ast", "Exc", "Inh", "Mic", "Oli", "OPC")
PRIORS <- c("uniform", "scEEMS_Weighted_Full", "scEEMS_Unweighted_Full", "scEEMS_Weighted_Restricted")
FMD    <- cfg_path("finemap_dir")
XCD    <- cfg_path("crosscell_dir")

genes <- fread(file.path(XCD, "genes.tsv"), header = FALSE, col.names = c("gene_id", "chr"))
stopifnot(idx >= 1, idx <= nrow(genes))
gene_id <- genes$gene_id[idx]; chr_str <- genes$chr[idx]
out_dir <- file.path(XCD, "per_gene")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_tsv <- file.path(out_dir, sprintf("%s.%s.tsv", gene_id, chr_str))
if (file.exists(out_tsv)) { message("exists, skipping: ", out_tsv); quit(save = "no", status = 0) }

# subset a susie object to variant indices, reindexing $sets (same function as step 10's
# finemap_and_coloc.R)
subset_susie_variants <- function(fit, keep_idx) {
  out <- list()
  out$alpha        <- fit$alpha[, keep_idx, drop = FALSE]
  out$lbf_variable <- fit$lbf_variable[, keep_idx, drop = FALSE]
  out$pip          <- fit$pip[keep_idx]
  out$V            <- fit$V
  out$lbf          <- fit$lbf
  out$sigma2       <- if (!is.null(fit$sigma2)) fit$sigma2 else 1
  if (!is.null(fit$sets) && length(fit$sets$cs) > 0) {
    pos_lookup <- match(seq_along(fit$pip), keep_idx)
    new_cs <- list(); new_idx <- integer(0)
    for (i in seq_along(fit$sets$cs)) {
      mapped <- pos_lookup[fit$sets$cs[[i]]]; mapped <- mapped[!is.na(mapped)]
      if (length(mapped) > 0) { new_cs[[length(new_cs) + 1]] <- mapped
                                new_idx <- c(new_idx, fit$sets$cs_index[i]) }
    }
    out$sets <- list(cs = new_cs, cs_index = new_idx,
                     coverage = fit$sets$coverage[match(new_idx, fit$sets$cs_index)],
                     purity   = fit$sets$purity[match(new_idx, fit$sets$cs_index), , drop = FALSE])
  } else out$sets <- list(cs = list(), cs_index = integer(0))
  class(out) <- "susie"; out
}

# one cell type's fit for this gene under one prior, with a status instead of an error
load_fit <- function(cell, prior) {
  f <- file.path(FMD, cell, paste0("fine_mapping_", prior),
                 sprintf("%s.%s.univariate_bvsr.rds", gene_id, chr_str))
  if (!file.exists(f)) return(list(status = "no fit"))
  o <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(o) || is.null(o$susie_fitted)) return(list(status = "unreadable fit"))
  s <- o$susie_fitted
  if (is.null(s$sets$cs) || length(s$sets$cs) == 0) return(list(status = "no credible set"))
  list(status = "ok", fit = s, canon = o$canon_keys)
}

note_row <- function(prior, cell, other = NA_character_, idx = NA_integer_, note)
  data.table(gene_id = gene_id, chr = chr_str, prior = prior, cell = cell, idx = idx,
             other_cell = other, other_idx = NA_integer_, note = note)

rows <- list(); add <- function(r) rows[[length(rows) + 1]] <<- r
for (prior in PRIORS) {
  fits <- setNames(lapply(CELLS, load_fit, prior = prior), CELLS)
  for (cell in CELLS) {
    a <- fits[[cell]]
    if (a$status != "ok") { add(note_row(prior, cell, note = a$status)); next }
    for (k in a$fit$sets$cs_index) add(note_row(prior, cell, idx = as.integer(k), note = "focal credible set"))
    for (other in setdiff(CELLS, cell)) {
      b <- fits[[other]]
      if (b$status != "ok") { add(note_row(prior, cell, other, note = paste(b$status, "in other cell"))); next }
      shared <- intersect(a$canon, b$canon)
      if (length(shared) < 10) { add(note_row(prior, cell, other, note = "fewer than 10 shared variants")); next }
      fa <- subset_susie_variants(a$fit, match(shared, a$canon))
      fb <- subset_susie_variants(b$fit, match(shared, b$canon))
      if (length(fa$sets$cs) == 0 || length(fb$sets$cs) == 0) {
        add(note_row(prior, cell, other, note = "credible set has no shared variants")); next }
      for (nm in c("alpha", "lbf_variable")) { colnames(fa[[nm]]) <- shared; colnames(fb[[nm]]) <- shared }
      names(fa$pip) <- shared; names(fb$pip) <- shared
      res <- tryCatch(suppressMessages(coloc::coloc.susie(fa, fb)),
                      error = function(e) { message(prior, " ", cell, "-", other, ": ", conditionMessage(e)); NULL })
      if (is.null(res) || is.null(res$summary)) { add(note_row(prior, cell, other, note = "coloc returned no result")); next }
      s <- as.data.table(res$summary)
      add(data.table(gene_id = gene_id, chr = chr_str, prior = prior, cell = cell, idx = as.integer(s$idx1),
                     other_cell = other, other_idx = as.integer(s$idx2), note = NA_character_,
                     nsnps = s$nsnps, PP.H0.abf = s$PP.H0.abf, PP.H1.abf = s$PP.H1.abf,
                     PP.H2.abf = s$PP.H2.abf, PP.H3.abf = s$PP.H3.abf, PP.H4.abf = s$PP.H4.abf))
    }
  }
}
out <- rbindlist(rows, fill = TRUE)
fwrite(out, out_tsv, sep = "\t")
message(sprintf("%s: %d rows (%d coloc rows) -> %s", gene_id, nrow(out), sum(is.na(out$note)), out_tsv))

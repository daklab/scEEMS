#!/usr/bin/env Rscript
# Summarize the per-gene eQTL fits of one cell type and prior into tables (read from the SuSiE objects,
# which carry purity, Bayes factors and coverage that the per-gene .cs.tsv does not):
#   {CELL}_{PRIOR}_gene.tsv          one row per fine-mapped gene, including genes without a credible set
#   {CELL}_{PRIOR}_cs.tsv            one row per credible set
#   {CELL}_{PRIOR}_cs_variants.tsv   one row per (credible set, variant) for variants in a credible set
#   {CELL}_{PRIOR}_pip010.tsv        every variant with PIP > 0.10, in a credible set or not
# A fit has two indexings: sets$cs, sets$coverage and sets$purity run over credible sets (1..nCS), while
# lbf, alpha, lbf_variable and V run over effects (1..L); sets$cs_index maps one to the other. cs_idx in
# the tables is the effect index, the same as coloc's idx1 and the release's CS_INDEX.
#
# usage:   Rscript aggregate_finemap.R CELL PRIOR [PIP_THRESHOLD=0.10]
# output:  {aggregate_dir}/finemapping_comparison/
suppressMessages({ library(data.table); library(parallel) })

args_all   <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(file.path(script_dir, "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript aggregate_finemap.R CELL PRIOR [PIP_THRESHOLD]")
CELL  <- args[1]
PRIOR <- args[2]
PIP_T <- if (length(args) >= 3) as.numeric(args[3]) else 0.10
stopifnot(CELL %in% CELLS, PRIOR %in% PRIORS)
OUT <- file.path(cfg_path("aggregate_dir"), "finemapping_comparison")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

in_dir <- file.path(cfg_path("finemap_dir"), CELL, paste0("fine_mapping_", PRIOR))
files <- list.files(in_dir, pattern = "univariate_bvsr\\.rds$", full.names = TRUE)
if (!length(files)) stop("no fits in ", in_dir)
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

one_gene <- function(path) {
  x <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(x) || is.null(x$susie_fitted)) return(NULL)
  s  <- x$susie_fitted
  ri <- x$region_info
  gene_id <- ri$gene_id; chr <- ri$chr
  vid <- x$variant_names; canon <- x$canon_keys
  pip <- s$pip
  cs  <- s$sets$cs

  gene_row <- data.table(
    gene_id = gene_id, chr = chr, cell_type = CELL, prior_label = PRIOR,
    n_variants_region = as.integer(ri$n_variants),
    n_cs = if (is.null(cs)) 0L else length(cs),
    max_pip = if (length(pip)) max(pip, na.rm = TRUE) else NA_real_,
    converged = isTRUE(s$converged), niter = if (is.null(s$niter)) NA_integer_ else as.integer(s$niter))

  hi <- which(pip > PIP_T)
  pip_dt <- if (!length(hi)) NULL else data.table(
    gene_id = gene_id, chr = chr, cell_type = CELL, prior_label = PRIOR,
    variant_id = vid[hi], canon_key = canon[hi], pip = pip[hi],
    in_cs = if (is.null(cs) || !length(cs)) FALSE else hi %in% unlist(cs))

  if (is.null(cs) || !length(cs)) return(list(gene = gene_row, cs = NULL, var = NULL, pip = pip_dt))

  Ls <- as.integer(s$sets$cs_index)          # effect index of each credible set
  cs_list <- vector("list", length(cs)); var_list <- vector("list", length(cs))
  for (j in seq_along(cs)) {
    idx <- cs[[j]]
    L   <- Ls[j]
    p_j <- pip[idx]
    top <- which.max(p_j)
    cs_list[[j]] <- data.table(
      gene_id = gene_id, chr = chr, cell_type = CELL, prior_label = PRIOR,
      cs_idx = L, cs_size = length(idx),
      coverage        = s$sets$coverage[j],
      min_abs_corr    = s$sets$purity$min.abs.corr[j],
      mean_abs_corr   = s$sets$purity$mean.abs.corr[j],
      median_abs_corr = s$sets$purity$median.abs.corr[j],
      lbf             = s$lbf[L],
      prior_variance  = if (is.null(s$V)) NA_real_ else s$V[L],
      top_pip         = p_j[top],
      sum_pip         = sum(p_j),
      top_variant_id  = vid[idx][top],
      top_canon_key   = canon[idx][top])
    var_list[[j]] <- data.table(
      gene_id = gene_id, chr = chr, cell_type = CELL, prior_label = PRIOR,
      cs_idx = L, variant_id = vid[idx], canon_key = canon[idx], pip = p_j,
      alpha        = if (is.null(s$alpha))        NA_real_ else s$alpha[L, idx],
      lbf_variable = if (is.null(s$lbf_variable)) NA_real_ else s$lbf_variable[L, idx],
      prior_weight = if (is.null(s$pi))           NA_real_ else s$pi[idx])
  }
  cs_dt <- rbindlist(cs_list)
  setorder(cs_dt, -top_pip)
  cs_dt[, cs_rank := seq_len(.N)]            # 1 = the gene's credible set with the highest top PIP
  var_dt <- rbindlist(var_list)
  var_dt[cs_dt, cs_rank := i.cs_rank, on = "cs_idx"]
  list(gene = gene_row, cs = cs_dt, var = var_dt, pip = pip_dt)
}

res <- mclapply(files, one_gene, mc.cores = ncores, mc.preschedule = TRUE)
bad <- vapply(res, is.null, logical(1))
res <- res[!bad]
gene_dt <- rbindlist(lapply(res, `[[`, "gene"))
cs_dt   <- rbindlist(Filter(Negate(is.null), lapply(res, `[[`, "cs")))
var_dt  <- rbindlist(Filter(Negate(is.null), lapply(res, `[[`, "var")))
pip_dt  <- rbindlist(Filter(Negate(is.null), lapply(res, `[[`, "pip")))

stem <- file.path(OUT, sprintf("%s_%s", CELL, PRIOR))
fwrite(gene_dt, paste0(stem, "_gene.tsv"),        sep = "\t")
fwrite(cs_dt,   paste0(stem, "_cs.tsv"),          sep = "\t")
fwrite(var_dt,  paste0(stem, "_cs_variants.tsv"), sep = "\t")
fwrite(pip_dt,  sprintf("%s_pip%03d.tsv", stem, round(PIP_T * 100)), sep = "\t")
cat(sprintf("[%s/%s] %d genes (%d unreadable fits), %d with a credible set; %d credible sets, %d variants in them\n",
            CELL, PRIOR, nrow(gene_dt), sum(bad), sum(gene_dt$n_cs > 0), nrow(cs_dt), nrow(var_dt)))

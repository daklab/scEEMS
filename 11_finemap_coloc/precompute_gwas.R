#!/usr/bin/env Rscript
# AD GWAS fine-mapping (susie_rss) for one gene's window, computed once per gene and shared by every cell
# type and eQTL prior. Two fits:
#   gwas_uniform   uniform prior
#   gwas_polyfun   functionally informed prior: per-SNP heritability (SNPVAR) estimated by PolyFun from
#                  the GWAS and 83 annotations (compute_snpvar.py). Variants without SNPVAR get the region
#                  minimum; a window without any SNPVAR falls back to the uniform fit.
# Both: susie_rss with L = 10, 95% coverage, refine = TRUE, lambda = 1e-3, residual variance estimated, n =
# median GWAS sample size, LD from precompute_ld.R.
#
# When a fit fails (susie_rss reports a negative residual variance or an unreasonably large prior
# variance), the window usually holds a few variants whose z-scores contradict the reference LD
# (segmental duplications such as MS4A or 16p11.2). Those variants (|z_std_diff| > 3 in susieR's
# kriging_rss) are dropped and both fits are redone on the remaining variants, if the window has at most
# RESCUE_MAX_VAR variants (default 25,000; kriging_rss is O(n^3)). Genes left without a fit can be rerun
# with RESCUE_MAX_VAR=100000 and more memory.
#
# usage:   Rscript precompute_gwas.R IDX            IDX = 1-based row of {finemap_dir}/genes.tsv
# output:  {gwas_cache_dir}/{gene}.{chr}.gwas.rds (both fits + the variant table they index); skipped if
#          it exists

suppressMessages({ library(susieR); library(Matrix); library(data.table) })

args_all   <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(file.path(script_dir, "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript precompute_gwas.R IDX")
idx <- as.integer(args[1])

LD_CACHE_DIR   <- cfg_path("ld_cache_dir")
GWAS_CACHE_DIR <- cfg_path("gwas_cache_dir")
SNPVAR_TABIX   <- file.path(cfg_path("snpvar_dir"), "bellenguez_sldsc83.snpvar.hg38.tsv.gz")
stopifnot(file.exists(SNPVAR_TABIX))     # without it the PolyFun fit would silently equal the uniform one
dir.create(GWAS_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
RESCUE_MAX_VAR <- as.integer(Sys.getenv("RESCUE_MAX_VAR", "25000"))

genes <- fread(file.path(cfg_path("finemap_dir"), "genes.tsv"), header = FALSE,
               col.names = c("gene_id", "chr_str"))
if (idx > nrow(genes)) stop("idx ", idx, " > number of genes ", nrow(genes))
gene_id <- genes$gene_id[idx]
chr_str <- genes$chr_str[idx]
chr_int <- as.integer(sub("chr", "", chr_str))

out_path <- file.path(GWAS_CACHE_DIR, sprintf("%s.%s.gwas.rds", gene_id, chr_str))
if (file.exists(out_path)) {
  cat(sprintf("[%s] GWAS fit exists, skipping\n", gene_id))
  quit(save = "no", status = 0)
}

# ----- LD cache -----
ld_path <- file.path(LD_CACHE_DIR, sprintf("%s.%s.ld.rds", gene_id, chr_str))
if (!file.exists(ld_path)) {
  message(sprintf("[%s] LD cache missing; run precompute_ld.R first: %s", gene_id, ld_path))
  quit(save = "no", status = 0)
}
ld     <- readRDS(ld_path)
joined <- ld$joined
R_gwas <- as.matrix(Matrix::forceSymmetric(ld$R_sparse, uplo = "U"))
R_gwas[is.na(R_gwas)] <- 0
diag(R_gwas) <- 1
if (nrow(joined) < 10) {
  message(sprintf("[%s] fewer than 10 variants; skipping", gene_id))
  quit(save = "no", status = 0)
}
n_for_susie <- as.integer(round(median(joined$N, na.rm = TRUE)))

fit_gwas <- function(prior_weights = NULL) {
  tryCatch(
    susie_rss(z = joined$z_oriented, R = R_gwas, n = n_for_susie,
              L = 10, coverage = 0.95, refine = TRUE, max_iter = 1000,
              estimate_residual_variance = TRUE, lambda = 1e-3, prior_weights = prior_weights),
    error = function(e) { message("susie_rss failed: ", conditionMessage(e)); NULL })
}
attach_names <- function(fit) {
  if (is.null(fit)) return(fit)
  colnames(fit$alpha)        <- joined$canon
  colnames(fit$lbf_variable) <- joined$canon
  names(fit$pip)             <- joined$canon
  fit
}

# ----- PolyFun prior: SNPVAR over the window, matched by canonical key -----
sv <- tryCatch(
  fread(cmd = sprintf("tabix %s %d:%d-%d", SNPVAR_TABIX, chr_int, ld$start_bp, ld$end_bp), header = FALSE,
        col.names = c("CHR", "POS", "A1", "A2", "SNPVAR")),
  error = function(e) { message("SNPVAR tabix failed: ", conditionMessage(e)); NULL })
prior_weights <- NULL
n_match <- 0L
if (!is.null(sv) && nrow(sv) > 0) {
  sv[, canon := mapply(canonical_key, CHR, POS, A1, A2)]
  sv <- sv[, .(SNPVAR = SNPVAR[1]), by = canon]          # one value per canonical key
  pw <- unname(setNames(sv$SNPVAR, sv$canon)[joined$canon])
  n_match <- sum(!is.na(pw))
  if (n_match > 0) {
    pw[is.na(pw)] <- min(pw, na.rm = TRUE)               # missing -> region minimum
    prior_weights <- pw / sum(pw)
  }
}

# ----- fits -----
gwas_uniform <- attach_names(fit_gwas(NULL))
if (is.null(gwas_uniform)) message(sprintf("[%s] uniform fit failed; continuing with the PolyFun prior", gene_id))
gwas_polyfun <- if (is.null(prior_weights)) gwas_uniform else attach_names(fit_gwas(prior_weights))

# ----- rescue: drop LD-inconsistent variants and refit -----
n_dropped <- 0L
if (is.null(gwas_uniform) || is.null(gwas_polyfun)) {
  if (nrow(joined) > RESCUE_MAX_VAR) {
    message(sprintf("[%s] fit failed; %d variants > RESCUE_MAX_VAR = %d, not rescued",
                    gene_id, nrow(joined), RESCUE_MAX_VAR))
  } else {
    kr <- tryCatch(susieR::kriging_rss(z = joined$z_oriented, R = R_gwas, n = n_for_susie),
                   error = function(e) { message("kriging_rss failed: ", conditionMessage(e)); NULL })
    if (!is.null(kr)) {
      d    <- abs(as.data.table(kr$conditional_dist)$z_std_diff)
      keep <- which(!is.na(d) & d <= 3)
      if (length(keep) >= 10 && length(keep) < nrow(joined)) {
        n_dropped <- nrow(joined) - length(keep)
        message(sprintf("[%s] dropping %d LD-inconsistent variants (max |z_std_diff| = %.1f)",
                        gene_id, n_dropped, max(d, na.rm = TRUE)))
        # subset the variant table, LD and prior together, then refit both priors on the same variants
        joined <- joined[keep, , drop = FALSE]
        R_gwas <- R_gwas[keep, keep, drop = FALSE]
        if (!is.null(prior_weights)) prior_weights <- prior_weights[keep] / sum(prior_weights[keep])
        gwas_uniform <- attach_names(fit_gwas(NULL))
        gwas_polyfun <- if (is.null(prior_weights)) gwas_uniform else attach_names(fit_gwas(prior_weights))
      }
    }
  }
}
if (is.null(gwas_polyfun)) gwas_polyfun <- gwas_uniform
if (is.null(gwas_polyfun)) {
  message(sprintf("[%s] both fits failed; skipping", gene_id))
  quit(save = "no", status = 0)
}

# ----- save (written to a temporary file and renamed, so a concurrent duplicate cannot corrupt it) -----
tmp_out <- paste0(out_path, ".tmp", Sys.getpid())
saveRDS(list(gene_id = gene_id, chr_str = chr_str, n_variants = nrow(joined), n_for_susie = n_for_susie,
             n_snpvar_matched = n_match, snpvar_used = !is.null(prior_weights),
             n_dropped = n_dropped, rescued = n_dropped > 0L,
             # the variants the fits index; coloc matches against these
             joined = joined[, .(canon, physical.pos, allele1, allele2, z_oriented, N)],
             gwas_uniform = gwas_uniform, gwas_polyfun = gwas_polyfun),
        tmp_out, compress = "xz")
invisible(file.rename(tmp_out, out_path))
cat(sprintf("[%s] %d variants, %d with SNPVAR%s -> %s\n", gene_id, nrow(joined), n_match,
            if (n_dropped > 0L) sprintf(", %d dropped by the rescue", n_dropped) else "", out_path))

#!/usr/bin/env Rscript
# Expand every colocalization (PP.H4 > 0.8 in coloc_all.tsv) into the full credible sets on both sides,
# one row per (colocalization, side, variant), for the mirrored PIP plots and the supplementary tables.
#   eQTL side: {CELL}_{EQTL_PRIOR}_cs_variants.tsv (aggregate_finemap.R)
#   GWAS side: the GWAS fit of precompute_gwas.R ($gwas_uniform / $gwas_polyfun); idx2 is an effect
#              index, so the credible set is sets$cs[[match(idx2, sets$cs_index)]].
# The GWAS fit depends only on the gene, so its rows repeat for every cell type colocalizing there.
# Added when available: the GWAS -log10(p) of every variant (both sides; tabix on the summary
# statistics), the marginal eQTL statistics of compute_eqtl_marginal.R and gene names from
# gene_mapping_file (gene_id, gene_name, gene_TSS).
#
# usage:   Rscript build_credset_table.R        (ALL_PAIRS=1: every tested pair, not only colocalizations)
# output:  {aggregate_dir}/finemapping_comparison/credset_all.tsv
suppressMessages({ library(data.table) })

args_all   <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(file.path(script_dir, "common.R"))

DIR       <- file.path(cfg_path("aggregate_dir"), "finemapping_comparison")
GWAS_DIR  <- cfg_path("gwas_cache_dir")
SUMSTATS  <- cfg_path("gwas_sumstats_file")
EQTL_P    <- file.path(cfg_path("finemap_dir"), "eqtl_marginal")
GENE_MAP  <- tryCatch(cfg_path("gene_mapping_file"), error = function(e) "")
SCRATCH   <- cfg_path("scratch_dir")
PP4_CUT   <- 0.8
ALL_PAIRS <- Sys.getenv("ALL_PAIRS", "0") == "1"

co  <- fread(file.path(DIR, "coloc_all.tsv"))
sel <- if (ALL_PAIRS) co else co[PP_H4 > PP4_CUT]
if (!nrow(sel)) stop("no coloc pairs pass PP.H4 > ", PP4_CUT, call. = FALSE)
# one colocalization = one (cell, priors, gene, eQTL effect, GWAS effect): a gene can colocalize through
# several pairs, and one GWAS credible set can pair with two eQTL credible sets
sel[, coloc_id := sprintf("%s|%s|%s|%s|%s|cs%d-%d", cell_type, eqtl_prior, gwas_id, gwas_prior,
                          gene_id, idx1, idx2)]
stopifnot(uniqueN(sel$coloc_id) == nrow(sel))

# ---- eQTL side ----
combos <- unique(sel[, .(cell_type, eqtl_prior)])
eq <- rbindlist(lapply(seq_len(nrow(combos)), function(i) {
  f <- file.path(DIR, sprintf("%s_%s_cs_variants.tsv", combos$cell_type[i], combos$eqtl_prior[i]))
  if (!file.exists(f)) { warning("missing ", basename(f), call. = FALSE); return(NULL) }
  want <- unique(sel[cell_type == combos$cell_type[i] & eqtl_prior == combos$eqtl_prior[i],
                     .(gene_id, cs_idx = idx1)])
  v <- fread(f, select = c("gene_id", "chr", "cell_type", "prior_label", "cs_idx", "variant_id",
                           "canon_key", "pip"))
  v[want, on = .(gene_id, cs_idx), nomatch = 0L]
}), use.names = TRUE)
setnames(eq, c("prior_label", "cs_idx"), c("eqtl_prior", "cs_idx_side"))

# ---- GWAS side: one fit per gene, both GWAS priors read from it ----
genes <- unique(sel[, .(gene_id, chr, gwas_id)])
gw <- rbindlist(lapply(seq_len(nrow(genes)), function(i) {
  gid <- genes$gene_id[i]; ch <- genes$chr[i]; arm <- genes$gwas_id[i]
  f <- file.path(GWAS_DIR, sprintf("%s.%s.gwas.rds", gid, ch))
  if (!file.exists(f)) { warning("no GWAS fit for ", gid, call. = FALSE); return(NULL) }
  cache <- readRDS(f)
  want <- unique(sel[gene_id == gid & gwas_id == arm, .(gwas_prior, cs_idx_side = idx2)])
  rbindlist(lapply(seq_len(nrow(want)), function(k) {
    fit <- cache[[paste0("gwas_", want$gwas_prior[k])]]
    if (is.null(fit) || is.null(fit$sets$cs)) return(NULL)
    j <- match(want$cs_idx_side[k], as.integer(fit$sets$cs_index))
    if (is.na(j)) { warning(sprintf("%s: idx2=%d not a credible set", gid, want$cs_idx_side[k]), call. = FALSE)
                    return(NULL) }
    mem <- fit$sets$cs[[j]]
    nm <- names(fit$pip)
    if (is.null(nm)) nm <- cache$joined$canon
    if (length(nm) != length(fit$pip)) return(NULL)
    data.table(gene_id = gid, chr = ch, gwas_id = arm, gwas_prior = want$gwas_prior[k],
               cs_idx_side = want$cs_idx_side[k], canon_key = nm[mem], pip = as.numeric(fit$pip[mem]))
  }), use.names = TRUE)
}), use.names = TRUE)

# ---- attach both sides to the colocalizations ----
key <- sel[, .(coloc_id, cell_type, eqtl_prior, gwas_id, gwas_prior, gene_id, idx1, idx2, PP_H4,
               cs_size_eqtl = cs_size, hit1, hit2)]
key[, `:=`(join_idx1 = idx1, join_idx2 = idx2)]
e <- merge(eq, key, by.x = c("cell_type", "eqtl_prior", "gene_id", "cs_idx_side"),
           by.y = c("cell_type", "eqtl_prior", "gene_id", "join_idx1"), allow.cartesian = TRUE)
e[, `:=`(side = "eqtl", is_coloc_lead = canon_key == hit1)]
g <- merge(gw, key, by.x = c("gene_id", "gwas_id", "gwas_prior", "cs_idx_side"),
           by.y = c("gene_id", "gwas_id", "gwas_prior", "join_idx2"), allow.cartesian = TRUE)
g[, `:=`(side = "gwas", is_coloc_lead = canon_key == hit2, variant_id = canon_key)]
setnames(e, "cs_idx_side", "cs_idx"); setnames(g, "cs_idx_side", "cs_idx")
cols <- c("coloc_id", "cell_type", "eqtl_prior", "gwas_id", "gwas_prior", "gene_id", "chr", "side", "cs_idx",
          "idx1", "idx2", "variant_id", "canon_key", "pip", "is_coloc_lead", "PP_H4", "cs_size_eqtl")
dt <- rbindlist(list(e[, ..cols], g[, ..cols]), use.names = TRUE)
sides <- dt[, .(n_side = uniqueN(side)), by = coloc_id]
if (any(sides$n_side != 2)) warning(sprintf("%d colocalizations have only one side", sides[n_side != 2, .N]),
                                    call. = FALSE)

dt[, c("k_chr", "k_pos", "a1", "a2") := tstrsplit(canon_key, ":", fixed = TRUE)]
dt[, `:=`(pos = as.integer(k_pos), k_chr = NULL, k_pos = NULL)]
stopifnot(!anyNA(dt$pos))
dt[, y := fifelse(side == "gwas", pip, -pip)]                 # mirror: GWAS up, eQTL down
dt[, cs_size_side := .N, by = .(coloc_id, side)]
setorder(dt, coloc_id, side, -pip)
dt[, pip_rank := seq_len(.N), by = .(coloc_id, side)]
dt[, is_top_pip := pip_rank == 1L]

# ---- gene names ----
if (nzchar(GENE_MAP) && file.exists(GENE_MAP)) {
  gm <- unique(fread(GENE_MAP, select = c("gene_id", "gene_name", "gene_TSS")))
  gm <- gm[!duplicated(gene_id)]
  dt[gm, `:=`(gene_name = i.gene_name, gene_TSS = i.gene_TSS), on = "gene_id"]
  dt[is.na(gene_name) | !nzchar(gene_name), gene_name := gene_id]
} else {
  dt[, gene_name := gene_id]
}

# ---- GWAS -log10(p) of every variant ----
if (nzchar(Sys.which("tabix"))) {
  dir.create(SCRATCH, recursive = TRUE, showWarnings = FALSE)
  p <- unique(dt[, .(c = sub("^chr", "", sub(":.*", "", canon_key)), pos)])
  bed <- file.path(SCRATCH, "credset_hits.bed")
  fwrite(p[order(c, pos), .(c, start = pos - 1L, end = pos)], bed, sep = "\t", col.names = FALSE)
  x <- fread(cmd = sprintf("tabix -R %s %s", shQuote(bed), shQuote(SUMSTATS)), header = FALSE,
             col.names = c("CHR", "POS", "REF", "ALT", "Z", "N", "MAF", "SNP"))
  if (nrow(x)) {
    x[, key := sprintf("chr%s:%d:%s:%s", CHR, POS, pmin(toupper(REF), toupper(ALT)), pmax(toupper(REF), toupper(ALT)))]
    x[, neglog10p := -(pnorm(-abs(Z), log.p = TRUE) + log(2)) / log(10)]
    s <- x[, .(neglog10p = max(neglog10p), rsid = SNP[1]), by = key][, gwas_id := GWAS_ID][]
    dt[s, `:=`(gwas_neglog10p = i.neglog10p, rsid = i.rsid), on = c("gwas_id", canon_key = "key")]
  }
} else {
  warning("tabix is not on PATH; GWAS p-values omitted", call. = FALSE)
}

# ---- marginal eQTL statistics ----
if (dir.exists(EQTL_P)) {
  need <- unique(dt[side == "eqtl", .(cell_type, gene_id, chr)])
  mg <- rbindlist(lapply(seq_len(nrow(need)), function(i) {
    f <- file.path(EQTL_P, need$cell_type[i], sprintf("%s.%s.marginal.tsv", need$gene_id[i], need$chr[i]))
    if (!file.exists(f)) return(NULL)
    cbind(fread(f, select = c("canon_key", "betahat", "sebetahat", "z", "neglog10p")),
          cell_type = need$cell_type[i], gene_id = need$gene_id[i])
  }), use.names = TRUE)
  if (nrow(mg)) {
    mg <- unique(mg, by = c("cell_type", "gene_id", "canon_key"))
    dt[mg, `:=`(eqtl_beta = i.betahat, eqtl_se = i.sebetahat, eqtl_z = i.z, eqtl_neglog10p = i.neglog10p),
       on = .(cell_type, gene_id, canon_key)]
  }
}

fwrite(dt, file.path(DIR, "credset_all.tsv"), sep = "\t")
cat(sprintf("wrote credset_all.tsv: %d rows for %d colocalizations\n", nrow(dt), uniqueN(dt$coloc_id)))

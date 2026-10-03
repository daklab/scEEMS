#!/usr/bin/env Rscript
# LD for one gene's GWAS fine-mapping, computed once per gene and shared by every cell type and prior.
#
#   1. Look up the gene's TADB cis window.
#   2. Read the AD GWAS summary statistics over the window (tabix; prep_gwas_sumstats.py).
#   3. Attach the ADSP European reference panel for the chromosome (bigsnpr backing; build_backing.R).
#   4. Match GWAS and panel variants by canonical key, orient the GWAS z-score to the panel's allele1, and
#      keep one row per canonical key (preferring a panel variant whose alleles match without flipping).
#   5. Compute the sparse LD matrix with bigsnpr::snp_cor (window 3 Mb, r2 >= 0.001); NaN entries
#      (monomorphic edge cases) are set to 0.
#
# usage:   Rscript precompute_ld.R IDX            IDX = 1-based row of {finemap_dir}/genes.tsv (make_gene_list.R)
# output:  {ld_cache_dir}/{gene}.{chr}.ld.rds  (joined variant table + R_sparse); skipped if it exists

suppressMessages({ library(bigsnpr); library(Matrix); library(data.table); library(readr) })

args_all   <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(file.path(script_dir, "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript precompute_ld.R IDX")
idx <- as.integer(args[1])

GWAS_TSVGZ   <- cfg_path("gwas_sumstats_file")
BACKING_DIR  <- cfg_path("adsp_backing_dir")
LD_CACHE_DIR <- cfg_path("ld_cache_dir")
dir.create(LD_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

genes <- fread(file.path(cfg_path("finemap_dir"), "genes.tsv"), header = FALSE,
               col.names = c("gene_id", "chr_str"))
if (idx > nrow(genes)) stop("idx ", idx, " > number of genes ", nrow(genes))
gene_id <- genes$gene_id[idx]
chr_str <- genes$chr_str[idx]
chr_int <- as.integer(sub("chr", "", chr_str))

out_path <- file.path(LD_CACHE_DIR, sprintf("%s.%s.ld.rds", gene_id, chr_str))
if (file.exists(out_path)) {
  cat(sprintf("[%s] LD cache exists, skipping\n", gene_id))
  quit(save = "no", status = 0)
}

# ----- TADB window -----
tadb <- read_delim(tadb_file(), delim = "\t", escape_double = FALSE, trim_ws = TRUE,
                   show_col_types = FALSE)
names(tadb)[1] <- "chr"
tadb_row <- tadb[tadb$gene_id == gene_id, ]
if (nrow(tadb_row) == 0) {
  message(sprintf("[%s] gene not in TADB; skipping", gene_id))
  quit(save = "no", status = 0)
}
start_bp <- as.integer(tadb_row$start[1])
end_bp   <- as.integer(tadb_row$end[1])

# ----- GWAS summary statistics over the window -----
gwas <- tryCatch(
  fread(cmd = sprintf("tabix %s %d:%d-%d", GWAS_TSVGZ, chr_int, start_bp, end_bp),
        col.names = c("CHR", "POS", "REF", "ALT", "Z", "N", "MAF", "SNP")),
  error = function(e) { message("GWAS tabix failed: ", conditionMessage(e)); NULL })
if (is.null(gwas) || nrow(gwas) < 10) {
  message(sprintf("[%s] too few GWAS variants; skipping", gene_id))
  quit(save = "no", status = 0)
}
gwas[, canon := mapply(function(p, a, r) canonical_key(chr_int, p, a, r), POS, ALT, REF)]

# ----- reference panel -----
bk_rds <- file.path(BACKING_DIR, sprintf("bigsnpr_ADSP_chr%d.rds", chr_int))
if (!file.exists(bk_rds)) stop(sprintf("[%s] backing file missing: %s (run build_backing.R)", gene_id, bk_rds))
bigSNP <- snp_attach(bk_rds)
G <- bigSNP$genotypes

map <- as.data.table(bigSNP$map)
map[, idx := .I]
map_region <- map[physical.pos >= start_bp & physical.pos <= end_bp]
if (nrow(map_region) < 10) {
  message(sprintf("[%s] too few panel variants in window; skipping", gene_id))
  quit(save = "no", status = 0)
}
map_region[, canon := mapply(function(p, a1, a2) canonical_key(chr_int, p, a1, a2),
                             physical.pos, allele1, allele2)]

# ----- match, orient, deduplicate -----
joined <- merge(map_region, gwas, by = "canon", suffixes = c("_adsp", "_bell"))
joined[, orient := fcase(allele1 == ALT & allele2 == REF, "keep",
                         allele1 == REF & allele2 == ALT, "flip",
                         default = "drop")]
joined[, z_oriented := fcase(orient == "keep", Z, orient == "flip", -Z, default = NA_real_)]
joined <- joined[!is.na(z_oriented)]
# multi-allelic panel sites can match one GWAS key twice; keep one row, preferring orient == "keep"
joined[, prefer := orient != "keep"]
setorder(joined, canon, prefer, physical.pos)
joined <- joined[, .SD[1], by = canon]
joined[, prefer := NULL]
setorder(joined, physical.pos)
stopifnot(!is.unsorted(joined$physical.pos))
if (nrow(joined) < 10) {
  message(sprintf("[%s] too few matched variants; skipping", gene_id))
  quit(save = "no", status = 0)
}

# ----- LD -----
options(bigstatsr.check.parallel.blas = FALSE)
slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK")
ncores_use <- if (nzchar(slurm_cpus)) as.integer(slurm_cpus) else nb_cores()
R_sparse <- snp_cor(G, ind.col = joined$idx, size = 3000, alpha = 1, thr_r2 = 0.001,
                    fill.diag = TRUE, infos.pos = joined$physical.pos, ncores = ncores_use)
R_sparse@x[is.na(R_sparse@x)] <- 0

saveRDS(list(gene_id = gene_id, chr_str = chr_str, chr_int = chr_int, start_bp = start_bp, end_bp = end_bp,
             n_variants = nrow(joined), joined = joined, R_sparse = R_sparse),
        out_path, compress = "xz")
cat(sprintf("[%s] %d variants -> %s\n", gene_id, nrow(joined), out_path))

#!/usr/bin/env Rscript
# Combine the per-(cell type, eQTL prior, GWAS prior) coloc tables of aggregate_coloc.R into
# coloc_all.tsv and add the GWAS association of both lead variants (hit1 = eQTL side, hit2 = GWAS side),
# looked up in the GWAS summary statistics by canonical key. -log10(p) is computed in log space, because
# 2 * pnorm(-|z|) underflows to 0 beyond |z| ~ 38. colocalized = PP.H4 > 0.8.
#
# usage:   Rscript build_coloc_table.R             (needs tabix on PATH)
# output:  {aggregate_dir}/finemapping_comparison/coloc_all.tsv
suppressMessages({ library(data.table) })

args_all   <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(file.path(script_dir, "common.R"))

DIR      <- file.path(cfg_path("aggregate_dir"), "finemapping_comparison")
SUMSTATS <- cfg_path("gwas_sumstats_file")
SCRATCH  <- cfg_path("scratch_dir")
PP4_CUT  <- 0.8
dir.create(SCRATCH, recursive = TRUE, showWarnings = FALSE)
if (!nzchar(Sys.which("tabix"))) stop("tabix is not on PATH", call. = FALSE)

files <- list.files(DIR, pattern = "_coloc\\.tsv$", full.names = TRUE)
files <- files[basename(files) != "coloc_all.tsv"]
dt <- rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)
cat(sprintf("%d coloc pairs from %d tables\n", nrow(dt), length(files)))

# one tabix -R pass over every hit position
pos <- unique(rbindlist(list(dt[, .(chr = hit1_chr, pos = hit1_pos)], dt[, .(chr = hit2_chr, pos = hit2_pos)])))
pos <- pos[!is.na(chr) & !is.na(pos)]
setorder(pos, chr, pos)
bed <- file.path(SCRATCH, "coloc_hits.bed")
fwrite(pos[, .(chr, start = pos - 1L, end = pos)], bed, sep = "\t", col.names = FALSE)
g <- fread(cmd = sprintf("tabix -R %s %s", shQuote(bed), shQuote(SUMSTATS)), header = FALSE,
           col.names = c("CHR", "POS", "REF", "ALT", "Z", "N", "MAF", "SNP"))
if (!nrow(g)) stop("tabix returned no rows from ", SUMSTATS, call. = FALSE)
g[, `:=`(a_lo = pmin(toupper(REF), toupper(ALT)), a_hi = pmax(toupper(REF), toupper(ALT)))]
g[, key := sprintf("chr%s:%d:%s:%s", CHR, POS, a_lo, a_hi)]
g[, neglog10p := -(pnorm(-abs(Z), log.p = TRUE) + log(2)) / log(10)]
g <- g[, .(neglog10p = max(neglog10p), Z = Z[which.max(abs(Z))], MAF = MAF[1], rsid = SNP[1]),
       by = key][, gwas_id := GWAS_ID][]

for (h in c("hit1", "hit2")) {
  dt[g, paste0(h, c("_neglog10p", "_z", "_maf", "_rsid")) := .(i.neglog10p, i.Z, i.MAF, i.rsid),
     on = setNames(c("gwas_id", "key"), c("gwas_id", h))]
}

# plot_* columns: the GWAS lead (hit2), the natural x-position when y is a GWAS p-value
dt[, `:=`(plot_variant   = hit2, plot_chr = hit2_chr, plot_pos = hit2_pos,
          plot_neglog10p = hit2_neglog10p, plot_rsid = hit2_rsid, plot_hit_side = "hit2",
          colocalized    = PP_H4 > PP4_CUT)]
fwrite(dt, file.path(DIR, "coloc_all.tsv"), sep = "\t")
cc <- dt[colocalized == TRUE]
cat(sprintf("wrote coloc_all.tsv: %d rows; PP.H4 > %.1f: %d pairs, %d (cell type, gene) pairs\n",
            nrow(dt), PP4_CUT, nrow(cc), uniqueN(cc[, .(cell_type, gene_id)])))

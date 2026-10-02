#!/usr/bin/env Rscript
# Convert one chromosome of the ADSP European reference panel (PLINK) into a bigsnpr file-backed matrix,
# which precompute_ld.R attaches to compute LD. Run once per chromosome; an existing pair is reused.
#
# usage:   Rscript build_backing.R CHR            (CHR = 1-22)
# input:   {adsp_plink_dir}/ADSP_EUR_chr{CHR}.{bed,bim,fam}
# output:  {adsp_backing_dir}/bigsnpr_ADSP_chr{CHR}.{bk,rds}  (the .bk file is tens of GB)

suppressMessages(library(bigsnpr))

args_all   <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(file.path(script_dir, "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript build_backing.R CHR")
chr_int <- as.integer(args[1])
stopifnot(!is.na(chr_int), chr_int >= 1, chr_int <= 22)

backing_dir <- cfg_path("adsp_backing_dir")
dir.create(backing_dir, recursive = TRUE, showWarnings = FALSE)
bed_file    <- file.path(cfg_path("adsp_plink_dir"), sprintf("ADSP_EUR_chr%d.bed", chr_int))
backingfile <- file.path(backing_dir, sprintf("bigsnpr_ADSP_chr%d", chr_int))
rds_file    <- paste0(backingfile, ".rds")
bk_file     <- paste0(backingfile, ".bk")
if (!file.exists(bed_file)) stop("PLINK file missing: ", bed_file)

if (file.exists(rds_file) && file.exists(bk_file)) {
  cat(sprintf("chr%d: backing files exist, skipping\n", chr_int))
  quit(save = "no", status = 0)
}
for (f in c(rds_file, bk_file)) if (file.exists(f)) file.remove(f)     # partial output of a failed run

options(bigstatsr.check.parallel.blas = FALSE)
snp_readBed2(bed_file, backingfile = backingfile, ncores = nb_cores())
bigSNP <- snp_attach(rds_file)
cat(sprintf("chr%d: %d samples x %d variants -> %s\n", chr_int, nrow(bigSNP$genotypes),
            ncol(bigSNP$genotypes), rds_file))

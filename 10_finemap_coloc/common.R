# Shared by the R scripts of step 10: configuration, fine-mapping prior labels, canonical variant keys and
# the locations of the eQTL input files.
#
#   args_all   <- commandArgs(trailingOnly = FALSE)
#   script_dir <- dirname(normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
#   source(file.path(script_dir, "common.R"))

source(file.path(script_dir, "..", "shared", "config.R"))

CELLS <- c("Mic", "Ast", "Exc", "Inh", "Oli", "OPC")

# eQTL fine-mapping priors. The three scEEMS priors are built from the step 6 predictions of one model.
PRIORS <- c("uniform", "EMS", "scEEMS_Weighted_Full", "scEEMS_Unweighted_Full", "scEEMS_Weighted_Restricted")
SCEEMS_MODEL <- c(scEEMS_Weighted_Full       = "weighted_full",
                  scEEMS_Unweighted_Full     = "unweighted_full",
                  scEEMS_Weighted_Restricted = "weighted_restricted")
GWAS_PRIORS <- c("uniform", "polyfun")
GWAS_ID <- "bellenguez"             # the AD GWAS (Bellenguez et al. 2022, stage 1)

# coloc output tag of a GWAS prior: "" for uniform, ".polyfun" for the PolyFun (SNPVAR) prior
gwas_tag <- function(gwas_prior) if (gwas_prior == "uniform") "" else paste0(".", gwas_prior)

# Canonical variant key chr<N>:<pos>:<allele>:<allele>, alleles sorted, so the eQTL data, the GWAS, the LD
# panel and the predictions join regardless of which allele each calls the reference.
canonical_key <- function(chr_int, pos, a1, a2) {
  pair <- sort(c(as.character(a1), as.character(a2)))
  sprintf("chr%d:%s:%s:%s", chr_int, as.character(pos), pair[1], pair[2])
}
pred_id_to_canon <- function(vid_vec) {          # "chrN:pos:REF:ALT" -> canonical key
  m <- stringr::str_match(vid_vec, "^chr([0-9XYMT]+):([0-9]+):([^:]+):([^:]+)$")
  vapply(seq_along(vid_vec), function(i) {
    if (is.na(m[i, 1])) return(NA_character_)
    canonical_key(as.integer(m[i, 2]), m[i, 3], m[i, 4], m[i, 5])
  }, character(1))
}

# eQTL inputs: FunGen-xQTL ROSMAP single-nucleus pseudobulk data, laid out as
#   {eqtl_data_dir}/genotype/ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.{N}.{bed,bim,fam}
#   {eqtl_data_dir}/{cell}/phenotype/snuc_pseudo_bulk.{cell}.mega.normalized.log2cpm.region_list.txt
#   {eqtl_data_dir}/{cell}/phenotype/phenotype_by_chrom/snuc_pseudo_bulk.{cell}.mega.normalized.log2cpm.bed.chr{N}.bed.gz
#   {eqtl_data_dir}/{cell}/covariate/snuc_pseudo_bulk.{cell}.mega.normalized.log2cpm.rosmap_cov.<...>.Marchenko_PC.gz
#   {eqtl_data_dir}/reference/TADB_enhanced_cis.bed       (cis window of each gene)
tadb_file <- function() file.path(cfg_path("eqtl_data_dir"), "reference", "TADB_enhanced_cis.bed")

eqtl_files <- function(cell, chr_str = NULL) {
  d <- cfg_path("eqtl_data_dir")
  out <- list(
    region_list = file.path(d, cell, "phenotype",
                            sprintf("snuc_pseudo_bulk.%s.mega.normalized.log2cpm.region_list.txt", cell)),
    tadb = tadb_file(),
    covariate = file.path(d, cell, "covariate", sprintf(paste0(
      "snuc_pseudo_bulk.%s.mega.normalized.log2cpm.rosmap_cov.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.",
      "plink_qc.snuc_pseudo_bulk_mega.related.plink_qc.extracted.pca.projected.Marchenko_PC.gz"), cell)))
  if (!is.null(chr_str)) {
    out$genotype <- file.path(d, "genotype", sprintf("ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.%d",
                                                     as.integer(sub("chr", "", chr_str))))
    out$phenotype <- file.path(d, cell, "phenotype", "phenotype_by_chrom",
                               sprintf("snuc_pseudo_bulk.%s.mega.normalized.log2cpm.bed.%s.bed.gz", cell, chr_str))
  }
  out
}

# eQTL data loading (pecotmr), shared by the fine-mapping and the marginal statistics so both see the same
# variants: residualized on the covariates, MAF >= 0.001, minor allele count >= 10, missingness <= 0.1.
load_eqtl_data <- function(cell, chr_str, region_id, association_id) {
  f <- eqtl_files(cell, chr_str)
  pecotmr::load_regional_univariate_data(
    genotype = f$genotype, phenotype = f$phenotype, covariate = f$covariate,
    region = region_id, association_window = association_id, conditions = c(cell),
    maf_cutoff = 0.001, mac_cutoff = 10, imiss_cutoff = 0.1,
    keep_indel = TRUE, scale_residuals = FALSE)
}

# Genes whose transcription start sites coincide are duplicate phenotype rows with identical expression;
# keep one column (stop if they ever differ).
single_phenotype <- function(Y_mat, gene_id) {
  if (ncol(Y_mat) > 1) {
    same <- vapply(seq_len(ncol(Y_mat)), function(j) isTRUE(all.equal(Y_mat[, j], Y_mat[, 1])), logical(1))
    if (!all(same)) stop(sprintf("co-located phenotypes differ for %s", gene_id))
    message(sprintf("[note] %d co-located duplicate phenotypes for %s; using column 1", ncol(Y_mat), gene_id))
    Y_mat <- Y_mat[, 1, drop = FALSE]
  }
  Y_mat
}

# Baseline SuSiE fine-mapping with uniform priors.
#
# Runs SuSiE fine-mapping without scEEMS-informed priors (uniform prior),
# serving as the baseline comparison for evaluating the effect of
# prediction-based priors.
#
# Usage:
#   Rscript template_uniform.R <gene_index> <cell_type>
#
# Arguments:
#   gene_index: 1-based index into the region list file
#   cell_type: Cell type abbreviation (e.g., Mic, Ast, Exc)

library(pecotmr)
library(tidyverse)
library(readr)
library(susieR)
library(yaml)

idx <- as.integer(commandArgs(trailingOnly = TRUE))[1]
cell_type <- commandArgs(trailingOnly = TRUE)[2]

# Load configuration
config <- yaml::read_yaml("../config.yaml")
data_dir <- config$paths$data_dir
fine_mapping_dir <- config$paths$fine_mapping_dir

predictions_dir <- file.path(data_dir, paste0("training_data/", cell_type, "_mega_eQTL/predictions_catboost"))
genotype_dir <- file.path(fine_mapping_dir, "genotype")
covariate_dir <- file.path(fine_mapping_dir, cell_type, "covariate")
phenotype_dir <- file.path(fine_mapping_dir, cell_type, "phenotype")
reference_dir <- file.path(fine_mapping_dir, "reference")

region_file <- file.path(phenotype_dir,
    paste0("snuc_pseudo_bulk.", cell_type, ".mega.normalized.log2cpm.region_list.txt"))
association_file <- file.path(reference_dir, "TADB_enhanced_cis.bed")

association_df <- read_delim(association_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    rename(chr = `#chr`)
region_df <- read_delim(region_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    rename(chr = `#chr`)

region_of_interest <- region_df[idx, ]
start_bp_region <- as.integer(region_of_interest$start)
end_bp_region <- as.integer(region_of_interest$end)
gene_ID <- region_of_interest$ID
chr <- region_of_interest$chr

association_of_interest <- association_df %>% filter(gene_id == gene_ID)
start_bp_association <- as.integer(association_of_interest$start[1])
end_bp_association <- as.integer(association_of_interest$end[1])

chr_int <- gsub("chr", "", chr)
region_id <- paste0(chr, ":", start_bp_region, "-", end_bp_region)
association_id <- paste0(chr, ":", start_bp_association, "-", end_bp_association)

genotype_file <- paste0(genotype_dir, "/ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.", chr_int)
phenotype_file <- file.path(phenotype_dir, "phenotype_by_chrom",
    paste0("snuc_pseudo_bulk.", cell_type, ".mega.normalized.log2cpm.bed.", chr, ".bed.gz"))
covariate_file <- file.path(covariate_dir, paste0("snuc_pseudo_bulk.", cell_type,
    ".mega.normalized.log2cpm.rosmap_cov.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.",
    "snuc_pseudo_bulk_mega.related.plink_qc.extracted.pca.projected.Marchenko_PC.gz"))

# Fine-mapping parameters
maf_cutoff <- 0.001
mac_cutoff <- 10
imiss_cutoff <- 0.1
max_L <- 30
coverage <- 0.95
seed <- 123
signal_cutoff <- 0.025

# Load data
tryCatch({
    fdat <- load_regional_univariate_data(
        genotype = genotype_file,
        phenotype = phenotype_file,
        covariate = covariate_file,
        region = region_id,
        association_window = association_id,
        conditions = c(cell_type),
        maf_cutoff = maf_cutoff,
        mac_cutoff = mac_cutoff,
        imiss_cutoff = imiss_cutoff,
        keep_indel = TRUE,
        scale_residuals = FALSE
    )
}, error = function(e) {
    message("Error loading data: ", e$message)
    quit(save = "no")
})

prior_probability <- 1 / 1000

X_colnames <- colnames(fdat$residual_X[[1]])
prior_vector <- rep(prior_probability, length(X_colnames))
names(prior_vector) <- X_colnames

X <- fdat$residual_X[[1]]
Y <- as.vector(fdat$residual_Y[[1]])
X_scalar <- fdat$X_scalar[[1]]
Y_scalar <- fdat$Y_scalar[[1]]
maf <- fdat$maf[[1]]

# SuSiE with uniform priors
set.seed(seed)

susie_wrapper <- function(X, y, init_L = 5, max_L = 30, l_step = 5, ...) {
  if (init_L == max_L) return(susie(X, y, L = init_L, ...))
  L <- init_L
  gst <- proc.time()
  while (TRUE) {
    st <- proc.time()
    res <- susie(X, y, L = L, ...)
    res$time_elapsed <- proc.time() - st
    if (!is.null(res$sets$cs)) {
      if (length(res$sets$cs) >= L && L <= max_L) {
        L <- L + l_step
      } else break
    } else break
  }
  message(paste("Total time elapsed for susie_wrapper:", (proc.time() - gst)[3]))
  return(res)
}

susie_fit <- susie_wrapper(
    X = X, y = Y, init_L = 5, max_L = max_L, l_step = 5,
    scaled_prior_variance = 0.2,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE,
    coverage = coverage,
    refine = TRUE,
    verbose = TRUE
)

# Organize output
region_info <- list(
    region_coord = parse_region(region_id),
    grange = parse_region(association_id),
    region_name = gene_ID
)

finemapping_output <- list(
    pip = setNames(susie_get_pip(susie_fit), X_colnames),
    cs = susie_get_cs(susie_fit, coverage = coverage, min_abs_corr = 0.5),
    variant_names = X_colnames,
    maf = fdat$maf[[1]],
    region_info = region_info,
    susie_fitted = susie_fit
)

susie_result_trimmed <- susie_post_processor(
    finemapping_output$susie_fitted, X, Y, X_scalar, Y_scalar, maf,
    signal_cutoff = signal_cutoff, min_abs_corr = 0.5, mode = "susie"
)

# Helper functions
extract_data <- function(i, susie_fitted, susie_result_trimmed) {
  variant_idx <- susie_fitted$sets$cs[[i]]
  num_variants <- length(variant_idx)
  tibble::tibble(
    variant_id = susie_result_trimmed$variant_names[variant_idx],
    betahat = susie_result_trimmed$sumstats$betahat[variant_idx],
    sebetahat = susie_result_trimmed$sumstats$sebetahat[variant_idx],
    z = susie_result_trimmed$sumstats$z_scores[variant_idx],
    pip = susie_fitted$pip[variant_idx],
    lbf_variable = susie_fitted$lbf_variable[i, variant_idx],
    lbf = rep(susie_fitted$lbf[i], num_variants),
    alpha = susie_fitted$alpha[i, variant_idx],
    coverage = rep(susie_fitted$sets$coverage[i], num_variants),
    min_abs_corr = rep(susie_fitted$sets$purity$min.abs.corr[i], num_variants),
    pi = susie_fitted$pi[variant_idx],
    CS_idx = paste0("CS_", i)
  )
}

retrieve_cs <- function(susie_fitted, susie_result_trimmed) {
  pip_sets <- susie_result_trimmed$susie_result_trimmed$sets$cs_index
  purrr::map_dfr(pip_sets, extract_data,
                 susie_fitted = susie_fitted,
                 susie_result_trimmed = susie_result_trimmed)
}

finemapping_output$susie_fitted_trimmed <- susie_result_trimmed

# Save full results
output_dir <- file.path(fine_mapping_dir, cell_type, "fine_mapping_uniform")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- paste0(output_dir, "/", gene_ID, ".", chr, ".univariate_bvsr.rds")
saveRDS(finemapping_output, output_file, compress = "xz")
message("Results saved to: ", output_file)

# Save top loci
top_output_dir <- file.path(fine_mapping_dir, cell_type, "fine_mapping_top_uniform")
top_pip_df <- retrieve_cs(susie_fitted = susie_fit, susie_result_trimmed = susie_result_trimmed)
top_pip_df$gene_id <- gene_ID
top_pip_df$chr <- chr
top_pip_df$job_idx <- idx

maf_df <- data.frame(variant_id = names(maf), maf = as.numeric(maf))
top_pip_df <- top_pip_df %>% left_join(maf_df, by = "variant_id")

dir.create(top_output_dir, recursive = TRUE, showWarnings = FALSE)
top_output_file <- paste0(top_output_dir, "/", gene_ID, ".", chr, ".univariate_bvsr_top_loci.tsv")
write_tsv(top_pip_df, top_output_file)

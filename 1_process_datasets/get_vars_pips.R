# Extract the variant PIPs and top loci of one gene from its SuSiE fine-mapping export.
#
# For each eQTL analysis in the gene's export (e.g. Mic_mega_eQTL, Mic_DeJager_eQTL, Mic_Kellis_eQTL;
# the DLPFC_Klein analyses are skipped), writes:
# - PIP_all: the PIP of every variant of the gene's cis window
# - PIP_top: the export's top loci table, with each variant's credible-set membership
#
# Usage:
#   Rscript get_vars_pips.R <gene_index>
#
# Arguments:
#   gene_index: 1-based line of {susie_dir}/genes_list.txt (create_job_list.py)
#
# Input:    {finemapping_rds_dir}/Fungen_xQTL.<gene_id>.cis_results_db.export{,_sumstats}.rds, which are not
#           distributed with this repository or the data release
# Outputs:  {susie_pips_dir}/PIP_all/<gene_id>_pips.csv and {susie_pips_dir}/PIP_top/<gene_id>_top_loci.csv,
#           one {susie_pips_dir} per eQTL analysis (by default {data_dir}/susie_vars_pips/<analysis>)
#
# Runs in the scEEMS_R environment (environment_r.yml).

library(tidyr)
library(dplyr)
library(tibble)

args_all <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(file.path(script_dir, "..", "shared", "config.R"))

args <- commandArgs(trailingOnly = TRUE)
idx <- args[1]

in_path <- cfg_path("finemapping_rds_dir")
out_path <- cfg_path("susie_dir")

if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

gene_list <- file.path(out_path, "genes_list.txt")
gene_ids <- read.table(gene_list, header = FALSE, stringsAsFactors = FALSE, col.names = "Ensembl_ID")

# Subset to the gene at the given index
gene_ids <- gene_ids[idx, ]

find_single_file <- function(base_dir, pattern) {
  hits <- list.files(base_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) {
    stop(paste("No file matched pattern:", pattern, "under", base_dir))
  }
  if (length(hits) > 1) {
    message(paste("Multiple matches found for", pattern, "- using first:", hits[1]))
  }
  return(hits[1])
}

file_PIP <- find_single_file(
  in_path,
  paste0("^Fungen_xQTL\\.", gene_ids, "\\.cis_results_db\\.export\\.rds$")
)
file_sumstats <- find_single_file(
  in_path,
  paste0("^Fungen_xQTL\\.", gene_ids, "\\.cis_results_db\\.export_sumstats\\.rds$")
)

list_PIP <- readRDS(file_PIP)
list_sumstats <- readRDS(file_sumstats)

cell_types <- names(list_PIP[[gene_ids]])

# Remove Klein DLPFC entries
cell_types <- cell_types[!grepl("DLPFC_Klein", cell_types)]

for (cell_type in cell_types) {
  print(cell_type)

  gene_top_loci <- list_PIP[[gene_ids]][[cell_type]][["top_loci"]]
  gene_pips <- list_PIP[[gene_ids]][[cell_type]][["pip"]]

  # Convert gene_pips to a dataframe
  if (!is.null(gene_pips)) {
    gene_pips <- as.data.frame(gene_pips)
    gene_pips <- gene_pips %>%
      rownames_to_column(var = "variant_id") %>%
      rename(pip = "gene_pips") %>%
      separate(variant_id, into = c("chr", "pos", "ref", "alt"), sep = ":", remove = FALSE) %>%
      mutate(gene_id = gene_ids)
  }

  if (!is.null(gene_top_loci)) {
    gene_top_loci <- gene_top_loci %>%
      separate(variant_id, into = c("chr", "pos", "ref", "alt"), sep = ":", remove = FALSE) %>%
      mutate(gene_id = gene_ids)
  }

  # Save results
  out_folder <- cfg_path("susie_pips_dir", cohort = cell_type)
  out_folder_PIP <- file.path(out_folder, "PIP_all")
  out_folder_top_PIP <- file.path(out_folder, "PIP_top")

  if (!dir.exists(out_folder_PIP)) {
    dir.create(out_folder_PIP, recursive = TRUE)
  }
  if (!dir.exists(out_folder_top_PIP)) {
    dir.create(out_folder_top_PIP, recursive = TRUE)
  }

  if (!is.null(gene_pips)) {
    write.table(gene_pips, file.path(out_folder_PIP, paste0(gene_ids, "_pips.csv")),
                row.names = FALSE, col.names = TRUE)
    if (!is.null(gene_top_loci)) {
      write.table(gene_top_loci, file.path(out_folder_top_PIP, paste0(gene_ids, "_top_loci.csv")),
                  row.names = FALSE, col.names = TRUE)
    }
  }
}

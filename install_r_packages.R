# R packages of the fine-mapping environment (environment_r.yml) that are not available from conda at the
# versions used for the manuscript: zlibbioc and snpStats from Bioconductor 3.21 (pecotmr reads PLINK files
# with snpStats), ebnm and flashier from the CRAN archive, and mr.mash.alpha, mvsusieR and pecotmr from
# GitHub at fixed commits. Their dependencies are in environment_r.yml.
#
#   conda env create -f environment_r.yml
#   conda activate scEEMS_R
#   R_LIBS_USER=NULL Rscript install_r_packages.R
#
# R_LIBS_USER=NULL keeps a personal R library out of the way, so the packages are built against, and
# installed into, the environment.

bioc <- c(zlibbioc = "1.54.0", snpStats = "1.58.0")
cran <- c(ebnm = "1.1-38", flashier = "1.0.7")
github <- c("stephenslab/mr.mash.alpha" = "7ec4b07c9cb297ec5e2a9df17a0b5021cea61020",
            "stephenslab/mvsusieR"      = "a152f4b3e848dcd98d6e31c5c1493bd2297a8cf9",
            "StatFunGen/pecotmr"        = "f7974d22fbfe20d5854b62a1de3b72bd4967e869")

for (p in names(bioc)) {
  remotes::install_url(sprintf("https://bioconductor.org/packages/3.21/bioc/src/contrib/%s_%s.tar.gz", p, bioc[[p]]),
                       dependencies = FALSE, upgrade = "never")
}
for (p in names(cran)) {
  remotes::install_version(p, version = cran[[p]], repos = "https://cloud.r-project.org",
                           dependencies = FALSE, upgrade = "never")
}
for (r in names(github)) {
  remotes::install_url(sprintf("https://github.com/%s/archive/%s.tar.gz", r, github[[r]]),
                       dependencies = FALSE, upgrade = "never")
}

for (p in c(names(bioc), names(cran), basename(names(github)))) {
  cat(sprintf("%-14s %s\n", p, as.character(packageVersion(p))))
}

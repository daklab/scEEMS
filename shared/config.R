# Settings for the R scripts of steps 10-11, read from config.yaml: the R counterpart of shared/config.py.
#
# The file is <repository>/config.yaml unless the SCEEMS_CONFIG environment variable names another one.
# Path settings are templates that may refer to {data_dir}, {release_dir}, {output_dir}, other path
# settings (e.g. {aggregate_dir}) and {cohort};
# settings config.yaml does not define take the defaults of shared/config.py (read from that file, so the
# two languages always agree).
#
#   source(file.path(repo_dir, "shared", "config.R"))
#   cfg_path("predictions_dir", cohort = "Mic_mega_eQTL")

suppressMessages(library(yaml))

.cfg_repo_dir <- function() {
  # the repository root = the parent of the directory holding this file (found via source()'s ofile)
  for (i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if (!is.null(f)) return(dirname(dirname(normalizePath(f))))
  }
  stop("source() shared/config.R with its path, e.g. source(file.path(script_dir, '..', 'shared', 'config.R'))")
}

SCEEMS_REPO_DIR <- .cfg_repo_dir()
SCEEMS_CONFIG_FILE <- Sys.getenv("SCEEMS_CONFIG", file.path(SCEEMS_REPO_DIR, "config.yaml"))

# DEFAULTS of shared/config.py, parsed from its dict literal: one '"key": "template",' entry per line
.cfg_defaults <- local({
  lines <- readLines(file.path(SCEEMS_REPO_DIR, "shared", "config.py"))
  m <- regmatches(lines, regexec('^\\s*"([A-Za-z0-9_]+)":\\s*"([^"]*)",?\\s*(#.*)?$', lines))
  m <- m[lengths(m) > 0]
  setNames(vapply(m, `[`, "", 3), vapply(m, `[`, "", 2))
})

.cfg <- NULL
cfg_load <- function() {
  if (is.null(.cfg)) {
    if (!file.exists(SCEEMS_CONFIG_FILE)) {
      stop(SCEEMS_CONFIG_FILE, " not found: copy config.yaml.example to config.yaml or set SCEEMS_CONFIG")
    }
    .cfg <<- yaml::read_yaml(SCEEMS_CONFIG_FILE)
  }
  .cfg
}

# Resolve a path setting: config.yaml paths[[key]] if set, else the shared/config.py default.
cfg_path <- function(key, cohort = NULL) {
  paths <- cfg_load()$paths
  template <- if (!is.null(paths[[key]])) paths[[key]] else unname(.cfg_defaults[key])
  if (is.null(template) || is.na(template)) {
    stop("path setting '", key, "' is not in ", SCEEMS_CONFIG_FILE, " and has no default")
  }
  for (k in unique(regmatches(template, gregexpr("\\{[A-Za-z0-9_]+\\}", template))[[1]])) {
    name <- substr(k, 2, nchar(k) - 1)
    value <- if (name == "cohort") {
      if (is.null(cohort)) stop("path setting '", key, "' needs a cohort")
      cohort
    } else if (name %in% c("data_dir", "release_dir", "output_dir")) {
      if (is.null(paths[[name]])) stop("path setting '", key, "' uses ", k, ", which is not set")
      paths[[name]]
    } else {
      cfg_path(name, cohort)
    }
    template <- gsub(k, value, template, fixed = TRUE)
  }
  path.expand(template)
}

cohort_name <- function(x) if (grepl("_mega_eQTL$", x)) x else paste0(x, "_mega_eQTL")

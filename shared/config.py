"""
Settings of every step of the pipeline, read from config.yaml.

The file is <repository>/config.yaml unless the SCEEMS_CONFIG environment variable names another one.
Every path setting is a template that may refer to {data_dir}, {release_dir}, {output_dir}, to other
path settings (e.g. {aggregate_dir}) and, for paths that depend on the cell type, {cohort} (for example
Mic_mega_eQTL). Settings that config.yaml does not define take the defaults below: steps 1-4 read and
write under data_dir, steps 5-11 read the scEEMS data release (release_dir, the folder downloaded from
Synapse) and write everything under output_dir.

    from config import path
    path("train_dir", cohort="Mic_mega_eQTL")   # -> <release_dir>/model_training/train/Mic_mega_eQTL
"""
import os
import re

import yaml

REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_FILE = os.environ.get("SCEEMS_CONFIG", os.path.join(REPO_DIR, "config.yaml"))

COHORTS = ["Ast_mega_eQTL", "Exc_mega_eQTL", "Inh_mega_eQTL", "Mic_mega_eQTL", "Oli_mega_eQTL", "OPC_mega_eQTL"]

DEFAULTS = {
    # ---- data release (model_training/, predictions/, fine_mapping/) ----
    "train_dir": "{release_dir}/model_training/train/{cohort}",
    "train_restricted_dir": "{release_dir}/model_training/train_restricted/{cohort}",
    "test_dir": "{release_dir}/model_training/test/{cohort}",
    "columns_dict_file": "{release_dir}/model_training/columns_dict/columns_dict.pkl",
    "gene_lof_file": "{release_dir}/model_training/gene_lof/41588_2024_1820_MOESM4_ESM.xlsx",
    "gnomad_maf_dir": "{release_dir}/model_training/gnomad_MAF",
    "gpn_star_file": "{release_dir}/model_training/gpn_star/gpn_star_scores_all.parquet",
    "feature_weights_file": "{release_dir}/model_training/feature_weights/best_configs_{cohort}.json",
    "published_models_dir": "{release_dir}/model_training/models/{cohort}",
    "release_predictions_dir": "{release_dir}/predictions/{cohort}",
    # ---- data release (featurization/): inputs for featurizing new variants (featurization/) ----
    "featurization_dir": "{release_dir}/featurization",
    "chrombpnet_models_dir": "{featurization_dir}/chrombpnet_models",
    "chrombpnet_peaks_file": "{featurization_dir}/chrombpnet_peaks.tsv.gz",
    "celltype_annotations_file": "{featurization_dir}/celltype_annotations.bed.gz",
    "abc_scores_file": "{featurization_dir}/abc_scores.tsv.gz",
    "genes_file": "{featurization_dir}/genes.tsv",
    "targets_file": "{featurization_dir}/targets_human.txt",
    # ---- steps 1-4: inputs and outputs of the featurization (not in the data release) ----
    "finemapping_rds_dir": "{data_dir}/release_04_2024",
    "susie_dir": "{data_dir}/susie_vars_pips",
    "variant_list_dir": "{susie_dir}/variant_list",
    "all_variants_dir": "{data_dir}/training_data/{cohort}/all_variants",
    "gene_list_dir": "{data_dir}/training_data/{cohort}",
    "training_sets_dir": "{data_dir}/training_data/{cohort}/training_data",
    # ---- outputs of steps 5-11 ----
    "feature_weight_search_dir": "{output_dir}/{cohort}/feature_weight_search",
    "model_dir": "{output_dir}/{cohort}/models",
    "test_predictions_dir": "{output_dir}/{cohort}/test_predictions",
    "published_models_eval_dir": "{output_dir}/{cohort}/published_models",
    "gene_list_file": "{output_dir}/{cohort}/list_genes_all.csv",
    "gpn_star_by_chr_dir": "{output_dir}/gpn_star_by_chr",
    "predictions_dir": "{output_dir}/{cohort}/predictions",
    "predictions_parquet_dir": "{output_dir}/{cohort}/predictions_parquet",
    "release_export_dir": "{output_dir}/release",
    "shap_dir": "{output_dir}/{cohort}/shap",
    "aggregate_dir": "{output_dir}/aggregate_results",
    "scratch_dir": "{output_dir}/scratch",
    # ---- steps 9-11 ----
    "susie_pips_dir": "{susie_dir}/{cohort}",
    "sldsc_dir": "{output_dir}/sldsc",
    "finemap_dir": "{output_dir}/fine_mapping",
    "ld_cache_dir": "{finemap_dir}/ld_cache",
    "gwas_cache_dir": "{finemap_dir}/gwas_cache",
    "adsp_backing_dir": "{finemap_dir}/adsp_backing",
    "snpvar_dir": "{finemap_dir}/snpvar",
    "gwas_sumstats_file": "{finemap_dir}/gwas/bellenguez_gcst90027158_hg38.tsv.gz",
    "crosscell_dir": "{output_dir}/crosscell_coloc",
}

_cfg = None


def load():
    """The parsed config.yaml (cached)."""
    global _cfg
    if _cfg is None:
        if not os.path.exists(CONFIG_FILE):
            raise FileNotFoundError(f"{CONFIG_FILE} not found: copy config.yaml.example to config.yaml "
                                    f"and fill in your paths, or set SCEEMS_CONFIG")
        with open(CONFIG_FILE) as fh:
            _cfg = yaml.safe_load(fh) or {}
    return _cfg


BASE_DIRS = ("data_dir", "release_dir", "output_dir")


def path(key, cohort=None):
    """Resolve a path setting: config.yaml paths[key] if set, else DEFAULTS[key], with its placeholders
    filled from the base directories, other path settings and `cohort`."""
    paths = load().get("paths", {}) or {}
    template = paths.get(key, DEFAULTS.get(key))
    if template is None:
        raise KeyError(f"path setting '{key}' is not in {CONFIG_FILE} and has no default")
    template = str(template)
    values = {}
    for name in set(re.findall(r"\{(\w+)\}", template)):
        if name == "cohort":
            if cohort is None:
                raise ValueError(f"path setting '{key}' needs a cohort")
            values[name] = cohort
        elif name in BASE_DIRS:
            if not paths.get(name):
                raise KeyError(f"path setting '{key}' uses {{{name}}}, which is not set in {CONFIG_FILE}")
            values[name] = paths[name]
        else:
            values[name] = path(name, cohort=cohort)
    return os.path.expanduser(template.format(**values))


def cohort_name(name):
    """'Mic' or 'Mic_mega_eQTL' -> 'Mic_mega_eQTL'."""
    return name if name.endswith("_mega_eQTL") else f"{name}_mega_eQTL"

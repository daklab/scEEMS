"""
Feature matrix shared by model training (step 5), inference (step 6) and SHAP (step 8), so the matrix a
model is trained on and the matrix it scores are built by the same code.

Feature contract (4,840 features, in this order):
  1. the columns_dict groups distance, ABC, celltype, baseline, chrombpnet_positive, diff and
     tf_positive (distance_TSS and abs_distance_TSS are dropped), then length_diff, is_SNP, is_indel,
     is_insertion, is_deletion, gene_lof and gnomad_MAF;
  2. absolute values of the diff (Enformer), tf_positive and chrombpnet_positive columns; infinite values
     and missing values set to 0;
  3. gpn_star_llr = |GPN-STAR log-likelihood ratio|, appended last. GPN-STAR scores only SNVs, so it is
     missing (NaN, handled natively by CatBoost) for insertions, deletions and other non-SNVs.

is_insertion is 1 when length_diff = len(REF) - len(ALT) > 0, and is_deletion when length_diff < 0, so the
two names are swapped relative to their meaning. The published models were trained with this definition;
both columns are in the variant_type category, so category weights and SHAP summaries are unaffected.

Typical use (training):
    aux = load_aux(); gpn = load_gpn_map()
    df = load_training_tables(cohort, "train", chroms, aux)
    X, FEATS, cols, abscols = build_X(df, aux["column_dict"], gpn)
    w = sample_weights(df)
    fw = feature_weights(weights, FEATS, aux["column_dict"])          # CatBoost feature_weights=
Typical use (inference), with the (cols, abscols) saved at training time so the column order matches:
    X, FEATS, _, _ = build_X(df_gene, aux["column_dict"], gpn, cols=cols, abscols=abscols)
"""
import os
import pickle
import re

import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.parquet as pq
from dask import dataframe as dd

from config import path

# ------------------------------------------------------------------ feature contract
SUBSET_KEYS = ["distance", "ABC", "celltype", "baseline", "chrombpnet_positive", "diff", "tf_positive"]
VARFEAT = ["length_diff", "is_SNP", "is_indel", "is_insertion", "is_deletion", "gene_lof", "gnomad_MAF"]
ABS_KEYS = ["diff", "tf_positive", "chrombpnet_positive"]
DROP_COLS = ["abs_distance_TSS", "distance_TSS"]
# columns that are never features (variant_id and gene_id stay in the frame as row keys)
META = ["variant_id", "pip", "CHR", "BP", "REF", "ALT", "SNP", "label", "weight", "_chrom", "gene_id"]
VARTYPE = ["is_SNP", "is_indel", "is_insertion", "is_deletion", "length_diff"]
CATS = ["CRE", "Enformer", "chromBPNet", "ABC", "TF", "abs_gpn", "gene_lof", "variant_type"]
GPN_COL = "gpn_star_llr"
SEED = 9448

# CatBoost hyperparameters, fixed for every model; only the feature-category weights are tuned
CONS = dict(depth=5, iterations=1000, learning_rate=0.03, l2_leaf_reg=5.0, min_data_in_leaf=10,
            bagging_temperature=1.0, leaf_estimation_method="Newton", leaf_estimation_iterations=10,
            loss_function="Logloss", verbose=False, random_seed=SEED)

# training data sets -> config path settings
LANES = {"train": "train_dir", "train_restricted": "train_restricted_dir", "test": "test_dir"}
ALL_CHR = [f"chr{i}" for i in range(1, 23)]


# ------------------------------------------------------------------ variant features
def make_variant_features(df):
    """_chrom, length_diff, is_SNP, is_indel, is_insertion, is_deletion from variant_id (chr:pos:ref:alt).

    is_insertion = length_diff > 0 and is_deletion = length_diff < 0: swapped relative to their meaning,
    kept because the models were trained this way (see the module docstring).
    """
    p = df["variant_id"].str.split(":", expand=True)
    ld = p[2].str.len() - p[3].str.len()          # len(ref) - len(alt)
    df["_chrom"] = p[0]
    df["length_diff"] = ld
    df["is_SNP"] = (ld == 0).astype(int)
    df["is_indel"] = (ld != 0).astype(int)
    df["is_insertion"] = (ld > 0).astype(int)
    df["is_deletion"] = (ld < 0).astype(int)
    return df


# ------------------------------------------------------------------ loaders
def load_aux(maf_chrom=None):
    """Gene constraint (gene_lof), gnomAD MAF and the columns dictionary. `maf_chrom` (e.g. 'chr7') loads
    only that chromosome's MAF table (per-gene inference); the default loads all 22 (training)."""
    glof = pd.read_excel(path("gene_lof_file"), "Supplementary Table 1")[["ensg", "post_mean"]]
    glof = glof.rename(columns={"ensg": "gene_id", "post_mean": "gene_lof"})
    glof["gene_lof"] = np.log2(glof["gene_lof"])
    maf_dir = path("gnomad_maf_dir")
    if maf_chrom is not None:
        maf = pd.read_csv(f"{maf_dir}/gnomad_MAF_{maf_chrom}.tsv", sep="\t")[["variant_id", "gnomad_MAF"]]
    else:
        maf = dd.read_csv(f"{maf_dir}/gnomad_MAF_chr*.tsv", sep="\t")[["variant_id", "gnomad_MAF"]].compute()
    with open(path("columns_dict_file"), "rb") as fh:
        column_dict = pickle.load(fh)
    return dict(glof=glof, maf=maf, column_dict=column_dict)


def load_gpn_map(chrom=None):
    """variant_id -> signed GPN-STAR LLR (SNVs only). `chrom` (e.g. 'chr7') loads only that chromosome:
    from the per-chromosome split written by 6_model_inference/split_gpn_by_chr.py if it exists, else
    by filtering the full file. The default loads all 15.3 million scores (training)."""
    if chrom is None:
        g = pd.read_parquet(path("gpn_star_file"), columns=["variant_id", GPN_COL])
    else:
        split = f"{path('gpn_star_by_chr_dir')}/{chrom}.parquet"
        if os.path.exists(split):
            g = pd.read_parquet(split, columns=["variant_id", GPN_COL])
        else:
            t = pq.read_table(path("gpn_star_file"), columns=["variant_id", GPN_COL])
            g = t.filter(pc.starts_with(t["variant_id"], f"{chrom}:")).to_pandas()
    return dict(zip(g["variant_id"].to_numpy(), g[GPN_COL].to_numpy()))


def _natural_key(p):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", p)]


def _list_parquet(cohort, lane, chroms):
    """Part files of annotated_data_{cohort}_{chr}.parquet for the given chromosomes, in natural order
    (chromosome, then part number). dask.read_parquet sorts the same way, so this order is also the
    row order the published models were trained on, whatever order the filesystem lists files in."""
    lane_dir = path(LANES[lane], cohort=cohort)
    files = []
    for ch in chroms:
        d = f"{lane_dir}/annotated_data_{cohort}_{ch}.parquet"
        if os.path.isdir(d):
            files += [os.path.join(dp, f) for dp, _, fns in os.walk(d) for f in fns if f.endswith(".parquet")]
        elif os.path.isfile(d):
            files.append(d)
    return sorted(files, key=_natural_key)


def _attach_aux(df, aux):
    for c in ("variant_id", "gene_id"):
        if c in df.columns:
            df[c] = df[c].astype(object)
    df = make_variant_features(df)
    df = df.merge(aux["glof"], on="gene_id", how="left").merge(aux["maf"], on="variant_id", how="left")
    df["gene_lof"] = df["gene_lof"].fillna(df["gene_lof"].median())
    df["gnomad_MAF"] = df["gnomad_MAF"].fillna(df["gnomad_MAF"].median())
    return df


def load_training_tables(cohort, lane, chroms, aux):
    """Read a training data set ('train', 'train_restricted' or 'test') for the given chromosomes and
    attach the variant features, gene_lof and gnomad_MAF. Missing chromosome files are skipped (microglia
    has no chr21 test or restricted file); returns an empty frame if there is nothing to read."""
    files = _list_parquet(cohort, lane, chroms)
    if not files:
        return pd.DataFrame()
    df = dd.read_parquet(files, engine="pyarrow").compute().reset_index(drop=True)
    return _attach_aux(df, aux)


def load_gene_variants(gene_path, gene_id, aux):
    """Read one gene's all_variants table (step 3 output) and attach the same features (inference)."""
    df = pd.read_parquet(gene_path).reset_index(drop=True)
    if "gene_id" not in df.columns:
        df["gene_id"] = gene_id
    return _attach_aux(df, aux)


# ------------------------------------------------------------------ feature matrix
def compute_cols(df, column_dict):
    """The feature columns present in df (order-stable) and the subset that is absolute-valued."""
    X = df.drop(columns=META, errors="ignore")
    sub = []
    for k in SUBSET_KEYS:
        sub += [c for c in column_dict.get(k, []) if c in X.columns]
    for vf in VARFEAT:
        if vf in X.columns and vf not in sub:
            sub.append(vf)
    sub = [c for c in sub if c not in DROP_COLS]
    abscols = [c for k in ABS_KEYS for c in column_dict.get(k, []) if c in X.columns and c in sub]
    return sub, abscols


def make_X(df, cols, abscols):
    """Reindex to cols (absent columns 0), categories -> numeric, inf -> 0, abs() on abscols, NaN -> 0."""
    Xs = df.reindex(columns=cols, fill_value=0).copy()
    for c in Xs.columns:
        if str(Xs[c].dtype) == "category":
            Xs[c] = pd.to_numeric(Xs[c].astype("object"), errors="coerce")
    Xs = Xs.replace([np.inf, -np.inf], 0)
    for c in abscols:
        if c in Xs.columns:
            Xs[c] = Xs[c].abs()
    return Xs.fillna(0)


def build_X(df, column_dict, gpn_map, cols=None, abscols=None):
    """Return (X[FEATS], FEATS, cols, abscols), with gpn_star_llr = |GPN-STAR LLR| appended last (NaN for
    variants without a score). Pass the (cols, abscols) saved at training time when scoring new data."""
    if cols is None:
        cols, abscols = compute_cols(df, column_dict)
    X = make_X(df, cols, abscols)
    gpn = pd.Series(np.abs(df["variant_id"].map(gpn_map).to_numpy()), index=X.index, name=GPN_COL)
    X = pd.concat([X, gpn], axis=1)
    FEATS = list(cols) + [GPN_COL]
    return X[FEATS], FEATS, cols, abscols


# ------------------------------------------------------------------ weights
def sample_weights(df):
    """Negatives 1; positive i gets pip_i * N_neg / sum(pip of positives), so both classes carry equal
    total weight and higher-PIP positives weigh more."""
    nneg = int((df["label"] == 0).sum())
    tot = df.loc[df["label"] == 1, "pip"].sum()
    pct = nneg / tot if tot > 0 else 1.0
    return np.where(df["label"] == 0, 1.0, df["pip"] * pct)


def category_map(FEATS, column_dict):
    """feature -> one of the 8 weighted categories (features outside them keep weight 1)."""
    fs = set(FEATS)
    catf = {
        "CRE": [c for c in column_dict.get("celltype", []) if c in fs],
        "Enformer": [c for c in column_dict.get("diff", []) if c in fs],
        "chromBPNet": [c for c in column_dict.get("chrombpnet_positive", []) if c in fs],
        "ABC": [c for c in column_dict.get("ABC", []) if c in fs],
        "TF": [c for c in column_dict.get("tf_positive", []) if c in fs],
        "abs_gpn": [GPN_COL],
        "gene_lof": [c for c in ["gene_lof"] if c in fs],
        "variant_type": [c for c in VARTYPE if c in fs],
    }
    return {f: cat for cat, fl in catf.items() for f in fl}


def feature_weights(weights, FEATS, column_dict):
    """{feature: weight of its category} for CatBoost's feature_weights; other features get 1.0.
    `weights` maps each of CATS to a weight (a "weights" entry of best_configs_{cohort}.json)."""
    f2c = category_map(FEATS, column_dict)
    return {f: float(weights.get(f2c.get(f), 1.0)) for f in FEATS}


NULL_WEIGHTS = {c: 1.0 for c in CATS}

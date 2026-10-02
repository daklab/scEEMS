#!/usr/bin/env python
"""
Feature-category weight search for one cell type on one chromosome parity (odd or even).

CatBoost feature weights multiply the score of every candidate split on a feature. The search tunes one
weight per feature category (featurize.CATS), each in [1, 100] on a log scale, with two objectives:
  * maximize AUPRC under leave-one-chromosome-out cross-validation within the parity: for each
    chromosome of the parity that has a test set, train on the other chromosomes of the same parity and
    score its test set; the AUPRC is computed on the pooled held-out predictions;
  * minimize complexity = sum over categories of (1 + log10 n_c) * |log10 w_c|, where n_c is the number of
    features in category c, so weighting few features by little is preferred.
Trials: 3 seed configurations (all weights 1; Enformer, chromBPNet, TF and GPN-STAR x10; the same with
GPN-STAR x50), 40 Sobol quasi-random trials, then 40 Gaussian-process (GPSampler) trials. The study is
kept in an SQLite database, so a rerun resumes where the last one stopped.

Weights selected on one parity are applied to the chromosomes of the other parity (train_loco.py), so
the search never sees the chromosome a model is evaluated on. Run both parities, then
select_feature_weights.py.

usage:   python search_feature_weights.py COHORT {odd|even} [N_SOBOL=40] [N_GP=40]
output:  {feature_weight_search_dir}/optuna_{cohort}_{parity}.db
"""
import math
import os
import sys
import time

import numpy as np
import optuna
from catboost import CatBoostClassifier
from optuna.samplers import GPSampler, QMCSampler, RandomSampler
from optuna.trial import TrialState
from sklearn.metrics import average_precision_score

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
import featurize as F
from config import cohort_name, path

optuna.logging.set_verbosity(optuna.logging.WARNING)
np.random.seed(F.SEED)

cohort = cohort_name(sys.argv[1])
parity = sys.argv[2].lower()
assert parity in ("odd", "even"), "parity must be odd or even"
N_SOBOL = int(sys.argv[3]) if len(sys.argv) > 3 else 40
N_GP = int(sys.argv[4]) if len(sys.argv) > 4 else 40

IN_ALL = [f"chr{i}" for i in range(1, 23) if i % 2 == (1 if parity == "odd" else 0)]
OUT_DIR = path("feature_weight_search_dir", cohort=cohort)
os.makedirs(OUT_DIR, exist_ok=True)

SEED_TRIALS = [
    {"CRE": 1., "Enformer": 1., "chromBPNet": 1., "ABC": 1., "TF": 1., "abs_gpn": 1., "gene_lof": 1., "variant_type": 1.},
    {"CRE": 1., "Enformer": 10., "chromBPNet": 10., "ABC": 1., "TF": 10., "abs_gpn": 10., "gene_lof": 1., "variant_type": 1.},
    {"CRE": 1., "Enformer": 10., "chromBPNet": 10., "ABC": 1., "TF": 10., "abs_gpn": 50., "gene_lof": 1., "variant_type": 1.},
]
N_SEED = len(SEED_TRIALS)

# ------------------------------------------------------------------ data: all 22 chromosomes, loaded once
t0 = time.time()
aux = F.load_aux()
gpn = F.load_gpn_map()
df = F.load_training_tables(cohort, "train", F.ALL_CHR, aux)
dft = F.load_training_tables(cohort, "test", F.ALL_CHR, aux)
X, FEATS, cols, abscols = F.build_X(df, aux["column_dict"], gpn)
Xte, _, _, _ = F.build_X(dft, aux["column_dict"], gpn, cols=cols, abscols=abscols)
y, w, chrom = df["label"].to_numpy(), F.sample_weights(df), df["_chrom"].to_numpy()
yte, chrom_te = dft["label"].to_numpy(), dft["_chrom"].to_numpy()
in_train = np.isin(chrom, IN_ALL)
IN_TEST = [c for c in IN_ALL if (chrom_te == c).sum() > 0]     # in-parity chromosomes with a test set
FOLDS = [(in_train & (chrom != j), chrom_te == j) for j in IN_TEST]
print(f"[{cohort} {parity}] train n={len(df):,} features={len(FEATS)} folds={IN_TEST} "
      f"(loaded in {time.time() - t0:.0f}s)", flush=True)

feat2cat = F.category_map(FEATS, aux["column_dict"])
NFEAT = {c: max(sum(1 for f in FEATS if feat2cat.get(f) == c), 1) for c in F.CATS}
GCOST = {c: 1.0 + math.log10(NFEAT[c]) for c in F.CATS}


def complexity(wts):
    return float(sum(GCOST[c] * abs(math.log10(wts[c])) for c in F.CATS))


def pooled_auprc(yy, pp):
    if len(yy) == 0 or yy.sum() == 0 or yy.sum() == len(yy):
        return float("nan")
    return float(average_precision_score(yy, pp))


def objective(trial):
    wts = {cat: trial.suggest_float(cat, 1.0, 100.0, log=True) for cat in F.CATS}
    fw = {f: wts.get(feat2cat.get(f), 1.0) for f in FEATS}
    ys, ps = [], []
    for tr, te in FOLDS:
        clf = CatBoostClassifier(**F.CONS, feature_weights=fw)
        clf.fit(X.loc[tr, FEATS], y[tr], sample_weight=w[tr])
        ps.append(clf.predict_proba(Xte.loc[te, FEATS])[:, 1])
        ys.append(yte[te])
    ap, cx = pooled_auprc(np.concatenate(ys), np.concatenate(ps)), complexity(wts)
    trial.set_user_attr("auprc_all", ap)
    trial.set_user_attr("complexity", cx)
    return ap, cx


study = optuna.create_study(study_name=f"{cohort}_{parity}_feature_weights",
                            storage=f"sqlite:///{OUT_DIR}/optuna_{cohort}_{parity}.db",
                            load_if_exists=True, directions=["maximize", "minimize"],
                            sampler=QMCSampler(qmc_type="sobol", scramble=True, seed=F.SEED))


def n_done():
    return len([t for t in study.trials if t.state == TrialState.COMPLETE])


def report(st, tr):
    phase = "seed" if tr.number < N_SEED else ("sobol" if tr.number < N_SEED + N_SOBOL else "gp")
    a = tr.user_attrs
    print(f"  trial {tr.number:3d} [{phase:5}] AUPRC={a.get('auprc_all', float('nan')):.5f} "
          f"complexity={a.get('complexity', float('nan')):.2f} | Pareto front {len(st.best_trials)}", flush=True)


done = n_done()
print(f"[{cohort} {parity}] {done} trials complete; target {N_SEED}+{N_SOBOL}+{N_GP}", flush=True)
if done < N_SEED:
    study.sampler = RandomSampler(seed=F.SEED)
    for s in SEED_TRIALS:
        study.enqueue_trial(s, skip_if_exists=True)
    study.optimize(objective, n_trials=N_SEED - done, callbacks=[report])
    done = n_done()
if done < N_SEED + N_SOBOL:
    study.sampler = QMCSampler(qmc_type="sobol", scramble=True, seed=F.SEED)
    study.optimize(objective, n_trials=(N_SEED + N_SOBOL) - done, callbacks=[report])
    done = n_done()
if done < N_SEED + N_SOBOL + N_GP:
    study.sampler = GPSampler(seed=F.SEED, n_startup_trials=0, deterministic_objective=False)
    study.optimize(objective, n_trials=(N_SEED + N_SOBOL + N_GP) - done, callbacks=[report])
print(f"[{cohort} {parity}] done: {n_done()} trials ({time.time() - t0:.0f}s)", flush=True)

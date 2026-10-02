#!/usr/bin/env python
"""
Select the feature-category weights of each cell type from its odd and even searches
(search_feature_weights.py).

For each parity:
  * "best" (used by train_loco.py) = the highest-AUPRC trial from the Gaussian-process phase (trial number
    >= 43), the model-guided part of the search, where a single lucky quasi-random point is less likely to
    win. If a search has not reached that phase, all completed trials are used.
  * "band_knee" (reference only) = the lowest-complexity Pareto-optimal trial within 0.005 AUPRC of "best".

usage:   python select_feature_weights.py [COHORT ...]            (default: all six cell types)
output:  {feature_weight_search_dir}/best_configs_{cohort}.json, in the format of the data release's
         model_training/feature_weights/ files. To train with your own selection, set
         paths.feature_weights_file in config.yaml to "{output_dir}/{cohort}/feature_weight_search/best_configs_{cohort}.json".
"""
import json
import os
import sys

import optuna
from optuna.trial import TrialState

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import COHORTS, cohort_name, path

optuna.logging.set_verbosity(optuna.logging.ERROR)
CATS = ["CRE", "Enformer", "chromBPNet", "ABC", "TF", "abs_gpn", "gene_lof", "variant_type"]
BAND = 0.005
GP_START = 3 + 40               # 3 seed + 40 Sobol trials; later trials come from the GPSampler

for coh in [cohort_name(c) for c in sys.argv[1:]] or COHORTS:
    wd = path("feature_weight_search_dir", cohort=coh)
    out = {"cohort": coh}
    for parity in ("odd", "even"):
        db = f"{wd}/optuna_{coh}_{parity}.db"
        if not os.path.exists(db):
            print(f"{coh} {parity}: no search database ({db})")
            continue
        storage = f"sqlite:///{db}"
        study = optuna.load_study(study_name=optuna.get_all_study_names(storage)[0], storage=storage)
        pts = []                # (complexity, AUPRC, weights, trial number)
        for t in study.get_trials(deepcopy=False, states=(TrialState.COMPLETE,)):
            if t.values and t.values[0] is not None and all(c in t.params for c in CATS):
                pts.append((float(t.values[1]), float(t.values[0]),
                            {c: float(t.params[c]) for c in CATS}, t.number))
        if not pts:
            continue
        pool = [p for p in pts if p[3] >= GP_START] or pts
        best = max(pool, key=lambda p: p[1])
        front = [p for p in pool if not any(q[1] >= p[1] and q[0] <= p[0] and q is not p for q in pool)]
        knee = min([p for p in front if p[1] >= best[1] - BAND], key=lambda p: p[0])
        out[parity] = {"n_trials": len(pts),
                       "best": {"auprc_in_parity": best[1], "complexity": best[0], "weights": best[2]},
                       "band_knee": {"auprc_in_parity": knee[1], "complexity": knee[0], "weights": knee[2]}}
    with open(f"{wd}/best_configs_{coh}.json", "w") as fh:
        json.dump(out, fh, indent=2)
    have = [p for p in ("odd", "even") if p in out]
    print(f"{coh}: {'both parities' if len(have) == 2 else f'only {have} -- train_loco.py needs both'}  "
          + "  ".join(f"{p} best AUPRC={out[p]['best']['auprc_in_parity']:.4f}" for p in have))

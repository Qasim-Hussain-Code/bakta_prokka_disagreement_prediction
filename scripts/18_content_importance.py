#!/usr/bin/env python3
"""Feature importance: permutation on held-out genomes, not impurity alone.

Reads:  results/metrics/16_content_cv_<model>.json
        data/interim/content_features.tsv.gz
Writes: results/metrics/18_content_importance_<model>.json

Impurity importance is reported because it is what a tree model hands you, but
it is not the ranking to read. It is computed on the training data and it
favours continuous, high-cardinality features over binary ones regardless of
whether the model relies on them. Several features here are correlated by
construction -- twenty amino acid fractions that sum to one, protein length
against interval length -- and impurity splits credit between them in a way
that misleads.

Permutation importance is measured on the held-out genomes, in MCC, so it
answers a question worth asking: how much test-set performance is lost when
this column is shuffled. It has its own failure mode with correlated features
-- shuffling one lets the model fall back on its partners, so the measured
drop understates every member of a correlated group. Correlations above 0.7
are listed next to each feature so a small number can be read for what it is.

Neither ranking is evidence about biology. The label is tool output.
"""

import json
import sys
from pathlib import Path

import numpy as np
from sklearn.inspection import permutation_importance
from sklearn.metrics import make_scorer, matthews_corrcoef

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_content

ROOT = Path(__file__).resolve().parent.parent
METRICS = ROOT / "results" / "metrics"

MODELS = ("decision_tree", "random_forest", "gradient_boosting")
N_REPEATS = 10
CORRELATION_THRESHOLD = 0.7


def correlations(X, names):
    """feature -> [(partner, r)] for |r| above the threshold."""
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = np.corrcoef(X, rowvar=False)
    corr = np.nan_to_num(corr)
    out = {}
    for i, name in enumerate(names):
        partners = [
            {"feature": names[j], "pearson_r": round(float(corr[i, j]), 3)}
            for j in range(len(names))
            if j != i and abs(corr[i, j]) >= CORRELATION_THRESHOLD
        ]
        out[name] = sorted(partners, key=lambda p: -abs(p["pearson_r"]))[:5]
    return out


def main():
    data = lib_content.prepare()
    names = data["names"]
    X_tr, y_tr = data["X_train"], data["y_train"]
    X_te, y_te = data["X_test"], data["y_test"]
    manifest = data["manifest"]["manifest"]

    corr = correlations(X_te, names)
    scorer = make_scorer(matthews_corrcoef)

    for model in MODELS:
        cv_path = METRICS / f"16_content_cv_{model}.json"
        if not cv_path.exists():
            print(f"  {cv_path.name} absent; skipping {model}")
            continue
        params = json.loads(cv_path.read_text())["selection"]["chosen_params"]
        fitted = lib_content.build(model, params).fit(X_tr, y_tr)

        impurity = getattr(fitted, "feature_importances_", None)
        perm = permutation_importance(
            fitted, X_te, y_te, scoring=scorer, n_repeats=N_REPEATS,
            random_state=lib_content.RANDOM_STATE, n_jobs=-1)

        ranked = []
        for i, name in enumerate(names):
            ranked.append({
                "feature": name,
                "group": manifest[name]["group"],
                "caller_derived": manifest[name]["caller_derived"],
                "db_derived": manifest[name]["db_derived"],
                "impurity_importance": (round(float(impurity[i]), 6)
                                        if impurity is not None else None),
                "permutation_importance_mean": round(float(perm.importances_mean[i]), 6),
                "permutation_importance_std": round(float(perm.importances_std[i]), 6),
                "correlated_with": corr[name],
            })
        ranked.sort(key=lambda r: -r["permutation_importance_mean"])

        by_flag = {}
        for flag in ("caller_derived", "db_derived"):
            sel = [r for r in ranked if r[flag]]
            by_flag[flag] = {
                "n_features": len(sel),
                "summed_permutation_importance": round(
                    sum(r["permutation_importance_mean"] for r in sel), 6),
                "top_3": [r["feature"] for r in sel[:3]],
            }
        neither = [r for r in ranked
                   if not r["caller_derived"] and not r["db_derived"]]
        by_flag["neither"] = {
            "n_features": len(neither),
            "summed_permutation_importance": round(
                sum(r["permutation_importance_mean"] for r in neither), 6),
            "top_3": [r["feature"] for r in neither[:3]],
        }

        payload = {
            "step": f"18_content_importance_{model}",
            "model": model,
            "hyperparameters": params,
            "scoring": "matthews_corrcoef",
            "n_repeats": N_REPEATS,
            "measured_on": "held-out genomes, not the training set",
            "correlation_threshold": CORRELATION_THRESHOLD,
            "caveats": {
                "read_permutation_not_impurity": (
                    "Impurity importance is computed on the training data and "
                    "favours continuous, high-cardinality features over binary "
                    "ones regardless of usefulness. It is reported for "
                    "completeness; the permutation column is the one to read."),
                "correlated_features": (
                    "Twenty amino acid fractions sum to one, and protein "
                    "length tracks interval length. When one member of a "
                    "correlated group is permuted the model falls back on its "
                    "partners, so each measured drop understates the group. "
                    "See correlated_with."),
                "not_causal": (
                    "These ranks describe what this model used. They are not "
                    "evidence about biology, and the label is tool output."),
            },
            "summary_by_flag": by_flag,
            "ranked_by_permutation": ranked,
        }
        out = METRICS / f"18_content_importance_{model}.json"
        out.write_text(json.dumps(payload, indent=2) + "\n")

        print(f"{model}")
        for r in ranked[:6]:
            print(f"  {r['feature']:32s} {r['group']:9s} "
                  f"perm {r['permutation_importance_mean']:+.5f}")
        print(f"  -> {out.relative_to(ROOT)}\n")


if __name__ == "__main__":
    main()

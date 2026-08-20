#!/usr/bin/env python3
"""Importances plus a permutation check. Correlated features split credit.

Reads:  data/interim/forest.joblib
        data/interim/features.tsv.gz
Writes: results/metrics/09_feature_importance.json

Two importance measures, because the cheap one is misleading on its own.

Impurity importance is what the forest reports for free. It is biased toward
features with many distinct values -- length_bp and dist_to_contig_end are
continuous and take thousands of values, while start_atg takes two. A binary
feature can matter more and still score lower.

Permutation importance asks a better question: shuffle one column in the
held-out data and see how much average precision falls. It is measured on
genomes the model never saw, so it reflects what the feature is worth for
generalisation rather than for fitting.

Neither survives correlation intact. gc_content, gc3 and genome_gc move
together; when one is permuted the forest leans on the others and the drop
understates all three. Correlated pairs are therefore reported alongside, so
a small importance can be read as "redundant" rather than "useless". Ranks
here are evidence about the model, not about biology.
"""

import json
from pathlib import Path

import joblib
import numpy as np
from sklearn.inspection import permutation_importance

import lib_model

ROOT = Path(__file__).resolve().parent.parent
MODEL_IN = ROOT / "data" / "interim" / "forest.joblib"
METRICS = ROOT / "results" / "metrics" / "09_feature_importance.json"

N_REPEATS = 10
CORRELATION_THRESHOLD = 0.7


def main():
    if not MODEL_IN.exists():
        raise SystemExit(f"{MODEL_IN} absent; run scripts/08_forest.py first")

    bundle = joblib.load(MODEL_IN)
    model, names = bundle["model"], bundle["names"]
    data = lib_model.prepare()
    if names != data["names"]:
        raise SystemExit(
            "feature set in the saved model does not match the current "
            "features table; re-run 08 before 09"
        )

    groups = {
        name: meta["group"]
        for name, meta in data["manifest"]["manifest"].items()
    }

    impurity = model.feature_importances_

    print(f"permutation importance, {N_REPEATS} repeats on held-out genomes")
    perm = permutation_importance(
        model, data["X_test"], data["y_test"],
        scoring="average_precision",
        n_repeats=N_REPEATS,
        random_state=lib_model.RANDOM_STATE,
        n_jobs=-1,
    )

    # Correlations among features, so a suppressed importance can be told
    # apart from an irrelevant feature.
    X = np.vstack([data["X_train"], data["X_test"]])
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = np.corrcoef(X, rowvar=False)
    corr = np.nan_to_num(corr)

    correlated = {}
    for i, name in enumerate(names):
        partners = [
            {"feature": names[j], "pearson_r": round(float(corr[i, j]), 3)}
            for j in range(len(names))
            if j != i and abs(corr[i, j]) >= CORRELATION_THRESHOLD
        ]
        if partners:
            correlated[name] = sorted(
                partners, key=lambda p: -abs(p["pearson_r"]))

    rows = []
    for i, name in enumerate(names):
        rows.append({
            "feature": name,
            "group": groups.get(name, "?"),
            "impurity_importance": round(float(impurity[i]), 6),
            "permutation_importance_mean": round(float(perm.importances_mean[i]), 6),
            "permutation_importance_std": round(float(perm.importances_std[i]), 6),
            "correlated_with": correlated.get(name, []),
        })
    by_perm = sorted(rows, key=lambda r: -r["permutation_importance_mean"])

    payload = {
        "step": "09_feature_importance",
        "scoring": "average_precision",
        "n_repeats": N_REPEATS,
        "measured_on": "held-out genomes, not the training set",
        "correlation_threshold": CORRELATION_THRESHOLD,
        "caveats": {
            "impurity_bias": (
                "Impurity importance favours continuous, high-cardinality "
                "features over binary ones regardless of usefulness."
            ),
            "correlated_features": (
                "Correlated features split credit. When one is permuted the "
                "forest falls back on its partners, so the measured drop "
                "understates each of them. See correlated_with."
            ),
            "not_causal": (
                "These ranks describe what this model used. They are not "
                "evidence about biology, and the label is tool output."
            ),
        },
        "ranked_by_permutation": by_perm,
    }
    METRICS.write_text(json.dumps(payload, indent=2) + "\n")

    print(f"\n{'feature':<22} {'group':<9} {'perm':>9} {'impurity':>9}")
    for r in by_perm[:12]:
        print(f"  {r['feature']:<20} {r['group']:<9} "
              f"{r['permutation_importance_mean']:>9.5f} "
              f"{r['impurity_importance']:>9.5f}")
    print(f"wrote {METRICS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Random forest, grouped split by genome, depth swept before test set is opened.

Reads:  data/interim/features.tsv.gz
Writes: data/interim/forest.joblib   (gitignored; step 09 reads it)
        results/metrics/08_forest.json

The order matters and is the reason this is its own step. max_depth is chosen
by cross-validation over the training genomes only. The held-out genomes are
scored once, afterwards, and that single number is the result. Sweeping depth
against the test set and then reporting the best one is how a model that has
learned nothing comes to look convincing.

This run uses every feature, including the gene-caller-derived group. That is
deliberate and the number is not the answer to the question. Step 10 removes
that group and reports both. Read them together.
"""

import json
from pathlib import Path

import joblib

import lib_model

ROOT = Path(__file__).resolve().parent.parent
MODEL_OUT = ROOT / "data" / "interim" / "forest.joblib"
METRICS = ROOT / "results" / "metrics" / "08_forest.json"


def main():
    data = lib_model.prepare()
    print(f"{len(data['names'])} features, "
          f"{len(data['train_genomes'])} train / {len(data['test_genomes'])} test genomes")
    print(f"train rows {len(data['y_train']):,} "
          f"(positive rate {data['y_train'].mean():.4f})")

    print("sweeping max_depth on training folds (test set untouched)")
    best_depth, sweep = lib_model.sweep_depth(
        data["X_train"], data["y_train"], data["fold_train"])
    for r in sweep:
        mark = " <-" if r["max_depth"] == best_depth else ""
        print(f"  depth {str(r['max_depth']):<5} "
              f"mean AP {r['mean_average_precision']}{mark}")

    model = lib_model.new_forest(best_depth).fit(data["X_train"], data["y_train"])
    MODEL_OUT.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump({"model": model, "names": data["names"]}, MODEL_OUT)

    overall, per_genome, _proba = lib_model.evaluate(
        model, data["X_test"], data["y_test"], data["genomes_test"])

    payload = {
        "step": "08_forest",
        "feature_set": "all groups, including gene-caller-derived",
        "feature_set_caveat": (
            "Includes the caller group, so this number does not answer "
            "'predictable from sequence alone'. See 10_circularity_audit."
        ),
        "n_features": len(data["names"]),
        "features": data["names"],
        "model": {
            "estimator": "RandomForestClassifier",
            "n_estimators": lib_model.N_ESTIMATORS,
            "class_weight": "balanced_subsample",
            "random_state": lib_model.RANDOM_STATE,
            "max_depth_grid": [str(d) for d in lib_model.DEPTH_GRID],
            "max_depth_selected": best_depth,
            "selection": "mean average precision over training folds only",
        },
        "split": data["split"],
        "depth_sweep": sweep,
        "test": overall,
        "test_per_genome": per_genome,
    }
    METRICS.write_text(json.dumps(payload, indent=2) + "\n")

    print(f"\nheld-out test: AP {overall['average_precision']} "
          f"(no-skill {overall['no_skill_average_precision']}), "
          f"ROC-AUC {overall['roc_auc']}")
    print(f"wrote {METRICS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

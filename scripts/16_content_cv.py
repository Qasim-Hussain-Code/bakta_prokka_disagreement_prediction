#!/usr/bin/env python3
"""Hyperparameter selection, on the training genomes only.

Reads:  data/interim/content_features.tsv.gz
        results/metrics/14_content_features.json
Writes: results/metrics/16_content_cv_<model>.json

Three model families, in increasing order of capacity: a single decision tree,
then a random forest, then gradient boosting. Each gets the same folds, the
same selection rule and the same seed, so the only thing that differs between
them is the model.

The grids are in lib_content.GRIDS and were written before the first fit. The
selection rule is mean MCC across five genome-grouped folds followed by the
one-standard-error rule -- among candidates within one standard error of the
best mean, take the simplest. Preferring the simplest model that is not
distinguishable from the best is what stops a sweep from chasing fold noise.

The test genomes are not read by this script at all.
"""

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_content

ROOT = Path(__file__).resolve().parent.parent
METRICS = ROOT / "results" / "metrics"

MODELS = ("decision_tree", "random_forest", "gradient_boosting")


def main():
    data = lib_content.prepare()
    X, y, folds = data["X_train"], data["y_train"], data["fold_train"]

    print(f"training rows {len(y):,}  positives {int(y.sum()):,} "
          f"({y.mean():.1%})  features {len(data['names'])}")
    print(f"cv folds: {lib_content.N_CV_FOLDS} genome groups, "
          f"test genomes untouched\n")

    for model in MODELS:
        t0 = time.time()
        chosen, results, selection = lib_content.sweep(model, X, y, folds)
        elapsed = round(time.time() - t0, 1)

        payload = {
            "step": f"16_content_cv_{model}",
            "model": model,
            "seed": lib_content.RANDOM_STATE,
            "selection": selection,
            "n_train_rows": int(len(y)),
            "n_train_positive": int(y.sum()),
            "n_features": len(data["names"]),
            "features": data["names"],
            "split": data["split_assertion"],
            "cv_folds": data["cv_folds"],
            "grid": lib_content.GRIDS[model],
            "results": sorted(results, key=lambda r: -r["mean_mcc"]),
            "elapsed_seconds": elapsed,
            "test_set_untouched": (
                "This script reads only the training genomes. The held-out "
                "genomes play no part in choosing hyperparameters."),
        }
        out = METRICS / f"16_content_cv_{model}.json"
        out.write_text(json.dumps(payload, indent=2) + "\n")

        print(f"{model}")
        print(f"  best      mean MCC {selection['best_mean_mcc']:.4f} "
              f"{selection['best_params']}")
        print(f"  1-SE      threshold {selection['one_se_threshold']:.4f}, "
              f"{selection['n_within_one_se']} of {selection['n_candidates']} within")
        print(f"  chosen    mean MCC {selection['chosen_mean_mcc']:.4f} "
              f"{selection['chosen_params']}"
              f"{'  (simpler than best)' if selection['chosen_is_simpler_than_best'] else ''}")
        print(f"  {elapsed}s -> {out.relative_to(ROOT)}\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Fit each selected model once and score it on the held-out genomes.

Reads:  results/metrics/16_content_cv_<model>.json   (chosen hyperparameters)
        data/interim/content_features.tsv.gz
Writes: results/metrics/17_content_final_<model>.json
        data/interim/content_model_<model>.joblib    (gitignored)
        data/interim/content_preds_<model>.npz       (gitignored)

One fit, one evaluation, on genomes that played no part in choosing anything.
Hyperparameters come from step 16's files rather than being re-selected here,
so the selection cannot drift towards the test set between the two steps.

Accuracy, precision, recall, F1, MCC and the confusion matrix in counts, for
every model. Counts and not only rates: on a near-balanced problem a rate can
be read as reassuring while the underlying cell is small.

For the forest, OOB is reported alongside the grouped CV mean. The two are not
competing estimates of the same thing and the file says so.
"""

import json
import sys
from pathlib import Path

import numpy as np
from joblib import dump

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_content

ROOT = Path(__file__).resolve().parent.parent
METRICS = ROOT / "results" / "metrics"
INTERIM = ROOT / "data" / "interim"

MODELS = ("decision_tree", "random_forest", "gradient_boosting")


def main():
    data = lib_content.prepare()
    X_tr, y_tr = data["X_train"], data["y_train"]
    X_te, y_te = data["X_test"], data["y_test"]

    baselines = json.loads((METRICS / "15_content_baselines.json").read_text())

    for model in MODELS:
        cv_path = METRICS / f"16_content_cv_{model}.json"
        if not cv_path.exists():
            raise SystemExit(f"{cv_path} absent; run scripts/16_content_cv.py first")
        cv = json.loads(cv_path.read_text())
        params = cv["selection"]["chosen_params"]

        fitted = lib_content.build(
            model, params, oob=(model == "random_forest")).fit(X_tr, y_tr)
        overall, per_genome, proba = lib_content.evaluate(
            fitted, X_te, y_te, data["genomes_test"])

        train_overall, _, _ = lib_content.evaluate(
            fitted, X_tr, y_tr, data["genomes_train"])

        payload = {
            "step": f"17_content_final_{model}",
            "model": model,
            "seed": lib_content.RANDOM_STATE,
            "hyperparameters": params,
            "hyperparameters_chosen_by": (
                f"step 16, mean MCC over {lib_content.N_CV_FOLDS} genome-grouped "
                "folds with the one-standard-error rule; test genomes untouched"),
            "grouped_cv_mean_mcc": cv["selection"]["chosen_mean_mcc"],
            "n_features": len(data["names"]),
            "split": data["split_assertion"],
            "test": overall,
            "per_genome_test": per_genome,
            "train_for_reference": {
                "note": ("in-sample, reported only so the gap to the test "
                         "numbers is visible; not a performance claim"),
                "accuracy": train_overall["accuracy"],
                "mcc": train_overall["mcc"],
                "f1": train_overall["f1"],
            },
            "baseline_comparison_same_test_rows": {
                name: {"mcc": b["test"]["mcc"], "accuracy": b["test"]["accuracy"],
                       "f1": b["test"]["f1"]}
                for name, b in baselines["baselines"].items()
                if b.get("computed", True)
            },
            "label_note": (
                "The label is disagreement between two software products about "
                "a name. A positive prediction is not a claim that either name "
                "is wrong, and not a claim about what the protein does."),
        }

        if model == "random_forest":
            oob = lib_content.oob_report(fitted, y_tr)
            payload["oob"] = oob
            payload["oob_vs_grouped_cv"] = {
                "oob_mcc": oob["oob_mcc"],
                "grouped_cv_mean_mcc": cv["selection"]["chosen_mean_mcc"],
                "gap": round(oob["oob_mcc"] - cv["selection"]["chosen_mean_mcc"], 5),
                "interpretation": (
                    "The gap is a measurement of genome-level leakage. OOB "
                    "resamples rows, so a held-out row is scored by trees that "
                    "saw its neighbours from the same genome and the same "
                    "annotation run. Grouped CV holds whole genomes out. The "
                    "difference is what that leakage is worth, not an error in "
                    "either estimate."),
            }
            dump(fitted, INTERIM / f"content_model_{model}.joblib")

        np.savez_compressed(
            INTERIM / f"content_preds_{model}.npz",
            y_true=y_te, proba=proba, pred=fitted.predict(X_te),
            genomes=data["genomes_test"])

        out = METRICS / f"17_content_final_{model}.json"
        out.write_text(json.dumps(payload, indent=2) + "\n")

        t = overall
        print(f"{model:20s} MCC {t['mcc']:>7.4f}  acc {t['accuracy']:.4f}  "
              f"P {t['precision']:.4f}  R {t['recall']:.4f}  F1 {t['f1']:.4f}")
        c = t["confusion_matrix_counts"]
        print(f"{'':20s} TN {c['true_negative']:>6,} FP {c['false_positive']:>6,} "
              f"FN {c['false_negative']:>6,} TP {c['true_positive']:>6,}")
        if model == "random_forest":
            print(f"{'':20s} OOB MCC {payload['oob']['oob_mcc']:.4f} vs grouped CV "
                  f"{cv['selection']['chosen_mean_mcc']:.4f} "
                  f"(gap {payload['oob_vs_grouped_cv']['gap']:+.4f})")
        print(f"{'':20s} -> {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

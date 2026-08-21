#!/usr/bin/env python3
"""The circularity audit: refit with the flagged features removed. Three arms.

Reads:  data/interim/content_features.tsv.gz
        results/metrics/14_content_features.json
Writes: results/metrics/19_content_circularity.json
        data/interim/content_preds_arm_<arm>.npz   (gitignored)

Part one's audit fired as designed and showed the model was reading the sORF
module's footprint rather than the sequence. The same audit is run here, with
one more arm, because this half has two ways of being circular rather than
one.

  full            every feature.

  no_caller       every feature except those encoding a gene-caller decision:
                  the chosen boundaries, the frame they imply, the strand, and
                  everything computed from the translated protein.

  sequence_only   caller-derived AND database-derived features both removed.
                  What is left is computed from DNA and from whole-genome
                  properties, and nothing else.

The db_derived arm is the one that matters most here. A model that predicts
name disagreement from Bakta's cross-references or from whether Prokka found a
similarity hit has learned which regions are thinly represented in db-light.
That is a fact about reference databases, not about DNA, and it would answer a
different question from the one on the tin while looking like an answer to
this one.

Which features go in which arm is read from the flags in
14_content_features.json. No feature list is hard-coded here.

Everything else is held identical between arms: same folds, same grid, same
selection rule, same seed. The only thing that varies is the feature set, so
the gaps are attributable to it.

Whichever way the numbers come out is what gets reported.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_content

ROOT = Path(__file__).resolve().parent.parent
METRICS = ROOT / "results" / "metrics" / "19_content_circularity.json"
INTERIM = ROOT / "data" / "interim"

MODEL = "random_forest"

ARMS = {
    "full":          {"drop_caller": False, "drop_db": False},
    "no_caller":     {"drop_caller": True,  "drop_db": False},
    "sequence_only": {"drop_caller": True,  "drop_db": True},
}


def main():
    baselines = json.loads(
        (ROOT / "results" / "metrics" / "15_content_baselines.json").read_text())

    arms = {}
    for arm, flags in ARMS.items():
        t0 = time.time()
        data = lib_content.prepare(**flags)
        names = data["names"]

        full_manifest = data["manifest"]["manifest"]
        dropped = sorted(set(full_manifest) - set(names))

        params, results, selection = lib_content.sweep(
            MODEL, data["X_train"], data["y_train"], data["fold_train"])
        fitted = lib_content.build(MODEL, params, oob=True).fit(
            data["X_train"], data["y_train"])
        overall, per_genome, proba = lib_content.evaluate(
            fitted, data["X_test"], data["y_test"], data["genomes_test"])

        np.savez_compressed(
            INTERIM / f"content_preds_arm_{arm}.npz",
            y_true=data["y_test"], proba=proba,
            pred=fitted.predict(data["X_test"]), genomes=data["genomes_test"])

        arms[arm] = {
            "n_features": len(names),
            "features": names,
            "n_dropped": len(dropped),
            "dropped_features": dropped,
            "hyperparameters": params,
            "grouped_cv_mean_mcc": selection["chosen_mean_mcc"],
            "selection": selection,
            "oob": lib_content.oob_report(fitted, data["y_train"]),
            "test": overall,
            "per_genome_test": per_genome,
            "elapsed_seconds": round(time.time() - t0, 1),
        }
        print(f"{arm:15s} {len(names):>3} features  "
              f"CV MCC {selection['chosen_mean_mcc']:.4f}  "
              f"test MCC {overall['mcc']:.4f}  acc {overall['accuracy']:.4f}",
              flush=True)

    full, nocall, seq = (arms["full"], arms["no_caller"], arms["sequence_only"])
    db_baseline = baselines["baselines"]["database_coverage_substitute"]["test"]["mcc"]
    length_baseline = baselines["baselines"]["protein_length_only"]["test"]["mcc"]

    payload = {
        "step": "19_content_circularity",
        "model": MODEL,
        "seed": lib_content.RANDOM_STATE,
        "arms_definition": {
            "full": "every feature",
            "no_caller": "caller_derived features removed",
            "sequence_only": "caller_derived and db_derived features both removed",
        },
        "feature_membership_read_from": (
            "the caller_derived and db_derived flags in "
            "14_content_features.json; no feature list is hard-coded in this "
            "script"),
        "held_identical_between_arms": (
            "folds, grid, selection rule, seed. Only the feature set varies, "
            "so the gaps are attributable to it."),
        "arms": arms,
        "side_by_side": {
            arm: {
                "n_features": a["n_features"],
                "grouped_cv_mean_mcc": a["grouped_cv_mean_mcc"],
                "test_mcc": a["test"]["mcc"],
                "test_accuracy": a["test"]["accuracy"],
                "test_f1": a["test"]["f1"],
                "test_precision": a["test"]["precision"],
                "test_recall": a["test"]["recall"],
                "test_roc_auc": a["test"]["roc_auc"],
                "confusion_matrix_counts": a["test"]["confusion_matrix_counts"],
            }
            for arm, a in arms.items()
        },
        "gaps_in_test_mcc": {
            "full_minus_no_caller": round(full["test"]["mcc"] - nocall["test"]["mcc"], 5),
            "no_caller_minus_sequence_only": round(
                nocall["test"]["mcc"] - seq["test"]["mcc"], 5),
            "full_minus_sequence_only": round(
                full["test"]["mcc"] - seq["test"]["mcc"], 5),
        },
        "reference_points_same_test_rows": {
            "majority_class_mcc": 0.0,
            "protein_length_only_mcc": length_baseline,
            "database_coverage_substitute_mcc": db_baseline,
        },
        "label_note": (
            "The label is disagreement between two software products about a "
            "name. None of these arms measures whether a name is correct."),
    }
    METRICS.write_text(json.dumps(payload, indent=2) + "\n")

    print(f"\ngaps in test MCC: "
          f"full-no_caller {payload['gaps_in_test_mcc']['full_minus_no_caller']:+.4f}, "
          f"no_caller-sequence_only "
          f"{payload['gaps_in_test_mcc']['no_caller_minus_sequence_only']:+.4f}")
    print(f"wrote {METRICS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

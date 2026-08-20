#!/usr/bin/env python3
"""Refit with every gene-caller-derived feature removed. Report both.

Reads:  data/interim/features.tsv.gz
        results/metrics/08_forest.json
Writes: results/metrics/10_circularity_audit.json

The question is whether disagreement is predictable from the sequence. Step
08 answers a weaker question, because it is allowed features that come from
Bakta's own decision rather than from DNA: the boundaries it chose, the frame
those imply, whether its database produced a name for the product.

is_hypothetical is the clearest case. A CDS Bakta could not name is exactly
the sort of call Prokka is also likely not to make. A model given it can
score well while having learned nothing about sequence, and the score would
still be real -- it would simply be answering "do the two tools' database
lookups agree", which is not the question on the tin.

This step refits with the whole caller group dropped, using the identical
procedure from lib_model: same folds, same depth grid, same seed, depth swept
on training folds only. The only difference is the feature set, so the gap
between the two numbers is attributable to the dropped group.

The sequence-only number is the one that answers the question. If it collapses
to the no-skill floor, the honest conclusion is that the disagreement is not
predictable from sequence with these features -- not that the experiment
failed.
"""

import json
from pathlib import Path

import numpy as np

import lib_model

ROOT = Path(__file__).resolve().parent.parent
FOREST_METRICS = ROOT / "results" / "metrics" / "08_forest.json"
METRICS = ROOT / "results" / "metrics" / "10_circularity_audit.json"
PRED_OUT = ROOT / "data" / "interim" / "preds_sequence_only.npz"

EXCLUDED_GROUP = "caller"


def main():
    data = lib_model.prepare(exclude_groups=(EXCLUDED_GROUP,))
    manifest = data["manifest"]
    dropped = manifest["features_by_group"].get(EXCLUDED_GROUP, [])

    print(f"refitting without the '{EXCLUDED_GROUP}' group "
          f"({len(dropped)} features dropped, {len(data['names'])} kept)")

    best_depth, sweep = lib_model.sweep_depth(
        data["X_train"], data["y_train"], data["fold_train"])
    for r in sweep:
        mark = " <-" if r["max_depth"] == best_depth else ""
        print(f"  depth {str(r['max_depth']):<5} "
              f"mean AP {r['mean_average_precision']}{mark}")

    model = lib_model.new_forest(best_depth).fit(data["X_train"], data["y_train"])
    overall, per_genome, proba = lib_model.evaluate(
        model, data["X_test"], data["y_test"], data["genomes_test"])
    np.savez(PRED_OUT, proba=proba, y=data["y_test"],
             genomes=data["genomes_test"])

    with_caller = json.loads(FOREST_METRICS.read_text())["test"] \
        if FOREST_METRICS.exists() else None

    comparison = None
    if with_caller:
        floor = overall["no_skill_average_precision"]
        full_ap = with_caller["average_precision"]
        seq_ap = overall["average_precision"]
        comparison = {
            "average_precision_all_features": full_ap,
            "average_precision_sequence_only": seq_ap,
            "no_skill_average_precision": floor,
            "absolute_drop": round(full_ap - seq_ap, 5),
            "lift_over_no_skill_all_features": round(full_ap - floor, 5),
            "lift_over_no_skill_sequence_only": round(seq_ap - floor, 5),
            "share_of_lift_lost": (
                round((full_ap - seq_ap) / (full_ap - floor), 4)
                if full_ap - floor > 0 else None
            ),
        }

    payload = {
        "step": "10_circularity_audit",
        "excluded_group": EXCLUDED_GROUP,
        "excluded_features": dropped,
        "kept_features": data["names"],
        "procedure_note": (
            "Identical to step 08 apart from the feature set: same folds, "
            "same depth grid, same seed, depth swept on training folds only. "
            "The difference between the two numbers is attributable to the "
            "dropped group."
        ),
        "model": {
            "max_depth_selected": best_depth,
            "max_depth_grid": [str(d) for d in lib_model.DEPTH_GRID],
            "n_estimators": lib_model.N_ESTIMATORS,
            "random_state": lib_model.RANDOM_STATE,
        },
        "split": data["split"],
        "depth_sweep": sweep,
        "test": overall,
        "test_per_genome": per_genome,
        "comparison": comparison,
        "reading_note": (
            "The sequence-only number is the one that answers the question. "
            "A collapse to the no-skill floor means the disagreement is not "
            "predictable from sequence with these features, which is a "
            "result, not a failure."
        ),
    }
    METRICS.write_text(json.dumps(payload, indent=2) + "\n")

    print(f"\nsequence-only test AP {overall['average_precision']} "
          f"(no-skill {overall['no_skill_average_precision']})")
    if comparison:
        print(f"all features      AP {comparison['average_precision_all_features']}")
        print(f"lift over no-skill: all {comparison['lift_over_no_skill_all_features']}, "
              f"sequence-only {comparison['lift_over_no_skill_sequence_only']}")
    print(f"wrote {METRICS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

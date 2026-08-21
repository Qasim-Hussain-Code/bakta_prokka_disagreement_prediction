#!/usr/bin/env python3
"""Baselines, scored on the same held-out genomes as every model.

Reads:  data/interim/content_features.tsv.gz
        results/metrics/14_content_features.json
Writes: results/metrics/15_content_baselines.json

Scored on the test genomes, not on the whole dataset. A baseline measured on a
different denominator from the models it is meant to calibrate is not a
baseline, and the comparison built on it would not be like for like.

Three baselines were specified:

  majority class            the floor any model has to clear.

  protein length only       one feature, one threshold, direction and cut
                            chosen on the training genomes only.

  is either product         SPECIFIED BUT DEGENERATE HERE, and reported as
  hypothetical              such rather than quietly replaced. The primary
                            analysis set is defined as both tools having
                            assigned a non-placeholder product, so this
                            predictor is constant 0 across every row it would
                            be scored on. It cannot be computed, and saying so
                            is the result.

The intent behind the third one -- if database coverage predicts naming
disagreement as well as a forest does, the answer is that this is a coverage
story rather than a sequence story -- still needs testing. The nearest
computable predictor that carries it is whether Prokka found a UniProtKB
similarity hit at all. That is added as a fourth baseline, labelled as a
substitute and not as the thing that was asked for.
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_content

ROOT = Path(__file__).resolve().parent.parent
METRICS = ROOT / "results" / "metrics" / "15_content_baselines.json"


class Constant:
    def __init__(self, value):
        self.value = int(value)

    def predict(self, X):
        return np.full(len(X), self.value, dtype=int)

    def predict_proba(self, X):
        p = np.full((len(X), 2), 0.0)
        p[:, self.value] = 1.0
        return p


class Threshold:
    """Predict 1 when feature <op> cut. One feature, one number."""

    def __init__(self, col, cut, above):
        self.col, self.cut, self.above = col, cut, bool(above)

    def predict(self, X):
        v = X[:, self.col]
        return ((v > self.cut) if self.above else (v <= self.cut)).astype(int)

    def predict_proba(self, X):
        pred = self.predict(X)
        p = np.zeros((len(X), 2))
        p[np.arange(len(X)), pred] = 1.0
        return p


def fit_threshold(X, y, col):
    """Choose cut and direction on the TRAINING rows only."""
    from sklearn.metrics import matthews_corrcoef
    values = np.unique(np.percentile(X[:, col], np.arange(1, 100)))
    best = (-2.0, None, None)
    for cut in values:
        for above in (True, False):
            mcc = matthews_corrcoef(y, Threshold(col, cut, above).predict(X))
            if mcc > best[0]:
                best = (float(mcc), float(cut), above)
    return Threshold(col, best[1], best[2]), best


def main():
    data = lib_content.prepare()
    names = data["names"]
    X_tr, y_tr = data["X_train"], data["y_train"]
    X_te, y_te = data["X_test"], data["y_test"]
    genomes_te = data["genomes_test"]

    results = {}

    # --- 1. majority class -------------------------------------------------
    majority = int(y_tr.mean() >= 0.5)
    overall, per_genome, _ = lib_content.evaluate(
        Constant(majority), X_te, y_te, genomes_te)
    results["majority_class"] = {
        "description": "always predict the majority class of the training genomes",
        "predicted_class": majority,
        "train_positive_rate": round(float(y_tr.mean()), 5),
        "test": overall,
        "per_genome": per_genome,
        "note": ("MCC is 0 by construction for a constant predictor. It is "
                 "reported so the floor is visible next to every other row."),
    }

    # --- 2. protein length, single threshold ------------------------------
    col = names.index("protein_length_aa")
    model, (train_mcc, cut, above) = fit_threshold(X_tr, y_tr, col)
    overall, per_genome, _ = lib_content.evaluate(model, X_te, y_te, genomes_te)
    results["protein_length_only"] = {
        "description": "one feature, one threshold",
        "feature": "protein_length_aa",
        "rule": f"predict disagreement when protein_length_aa "
                f"{'>' if above else '<='} {cut:g}",
        "threshold": cut,
        "direction": "above" if above else "at or below",
        "chosen_on": "training genomes only",
        "train_mcc": round(train_mcc, 5),
        "test": overall,
        "per_genome": per_genome,
    }

    # --- 3. is either product hypothetical: specified, degenerate ----------
    results["is_either_product_hypothetical"] = {
        "description": ("specified baseline: predict disagreement when either "
                        "tool's product is a placeholder"),
        "computed": False,
        "reason": (
            "Degenerate by construction. The primary analysis set is defined "
            "as both tools having assigned a non-placeholder product, so this "
            "predictor takes the value 0 on every row it would be scored on. "
            "It is identical to the majority-class baseline with the minority "
            "class predicted, and carries no information about anything."),
        "not_silently_replaced": (
            "Reported as uncomputable rather than swapped for something that "
            "would produce a number. The substitute below is labelled as a "
            "substitute."),
        "where_this_predictor_does_live": (
            "It is the cohort definition itself. Its content is the "
            "cross-tabulation in 12_content_cohort.json: 8,653 pairs where "
            "both tools produced a placeholder, 26,699 where only Bakta named "
            "the region, 2,297 where only Prokka did."),
    }

    # --- 4. database coverage: substitute carrying the same intent --------
    col = names.index("prokka_has_similarity_hit")
    model, (train_mcc, cut, above) = fit_threshold(X_tr, y_tr, col)
    overall, per_genome, _ = lib_content.evaluate(model, X_te, y_te, genomes_te)
    results["database_coverage_substitute"] = {
        "description": ("substitute for the degenerate baseline above: predict "
                        "from whether Prokka found a UniProtKB similarity hit"),
        "is_a_substitute": True,
        "what_it_tests": (
            "If this scores close to the forest, naming disagreement is mostly "
            "a database-coverage story and not a sequence story. That was the "
            "point of the baseline it replaces."),
        "feature": "prokka_has_similarity_hit",
        "feature_is_db_derived": True,
        "rule": f"predict disagreement when prokka_has_similarity_hit "
                f"{'>' if above else '<='} {cut:g}",
        "chosen_on": "training genomes only",
        "train_mcc": round(train_mcc, 5),
        "test": overall,
        "per_genome": per_genome,
    }

    payload = {
        "step": "15_content_baselines",
        "scored_on": "held-out test genomes only, the same rows every model is scored on",
        "denominator_note": (
            "Baselines and models share one denominator. Scoring baselines on "
            "the whole dataset and models on a held-out split would make every "
            "comparison downstream meaningless."),
        "n_test_rows": int(len(y_te)),
        "n_test_positive": int(y_te.sum()),
        "test_positive_rate": round(float(y_te.mean()), 5),
        "test_genomes": data["test_genomes"],
        "split": data["split_assertion"],
        "baselines": results,
    }
    METRICS.write_text(json.dumps(payload, indent=2) + "\n")

    print(f"test set {len(y_te):,} rows, {int(y_te.sum()):,} positive "
          f"({y_te.mean():.1%})")
    for name, r in results.items():
        if not r.get("computed", True):
            print(f"  {name:34s} not computable (degenerate by construction)")
            continue
        t = r["test"]
        print(f"  {name:34s} MCC {t['mcc']:>7.4f}  acc {t['accuracy']:.4f}  "
              f"F1 {t['f1']:.4f}")
    print(f"wrote {METRICS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Length-only, GC-only, majority class. Nothing beats these by default.

Reads:  data/interim/features.tsv.gz
Writes: results/metrics/07_baselines.json

These exist so the forest in step 08 has something to beat. A model that
cannot beat "always predict the majority class" has learned nothing, and on
an imbalanced label that failure is easy to hide behind a high accuracy.

Accuracy is therefore not reported as a headline number. With a positive rate
in the low single digits, always answering "no disagreement" scores extremely
well on accuracy while being useless. Average precision is the number to
read: its no-skill value is the positive rate itself, so a model that has
learned nothing cannot look good.

Every baseline is fitted on the training genomes and scored on the held-out
genomes from lib_split, the same ones step 08 uses.
"""

import gzip
import json
from pathlib import Path

import numpy as np
from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import lib_split

ROOT = Path(__file__).resolve().parent.parent
FEATURES = ROOT / "data" / "interim" / "features.tsv.gz"
METRICS = ROOT / "results" / "metrics" / "07_baselines.json"


def load_table(path):
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh]
    return header, rows


def evaluate(name, scores, y_true, note=""):
    """A single scored baseline.

    roc_auc is undefined when one class is absent, which can legitimately
    happen if a held-out genome has no positive at all.
    """
    out = {
        "baseline": name,
        "average_precision": round(float(average_precision_score(y_true, scores)), 5),
        "roc_auc": (round(float(roc_auc_score(y_true, scores)), 5)
                    if len(set(y_true)) > 1 else None),
    }
    if note:
        out["note"] = note
    return out


def main():
    if not FEATURES.exists():
        raise SystemExit(f"{FEATURES} absent; run 06 first")

    header, rows = load_table(FEATURES)
    idx = {c: i for i, c in enumerate(header)}

    genome_gc = {}
    for r in rows:
        genome_gc.setdefault(r[idx["genome"]], float(r[idx["genome_gc"]]))
    train_genomes, test_genomes, _folds = lib_split.split_genomes(genome_gc)
    train_set, test_set = set(train_genomes), set(test_genomes)

    def subset(keep):
        sel = [r for r in rows if r[idx["genome"]] in keep]
        y = np.array([int(r[idx["label"]]) for r in sel])
        return sel, y

    train_rows, y_train = subset(train_set)
    test_rows, y_test = subset(test_set)

    def column(sel, name):
        return np.array([[float(r[idx[name]])] for r in sel])

    results = []

    # Majority class. Predicts the constant negative; its average precision is
    # the positive rate of the test set, which is the no-skill floor.
    dummy = DummyClassifier(strategy="most_frequent").fit(
        column(train_rows, "length_bp"), y_train)
    results.append(evaluate(
        "majority_class",
        dummy.predict_proba(column(test_rows, "length_bp"))[:, 1],
        y_test,
        "constant prediction; average_precision here is the no-skill floor",
    ))

    # Single-feature logistic regressions. Scaled, because raw length in bp
    # and GC as a fraction differ by four orders of magnitude.
    for feature in ("length_bp", "gc_content"):
        model = make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=1000, class_weight="balanced"),
        ).fit(column(train_rows, feature), y_train)
        results.append(evaluate(
            f"{feature}_only",
            model.predict_proba(column(test_rows, feature))[:, 1],
            y_test,
        ))

    payload = {
        "step": "07_baselines",
        "headline_metric": "average_precision",
        "headline_metric_note": (
            "Accuracy is not reported as a headline. At this positive rate a "
            "constant negative prediction scores well on accuracy and is "
            "useless. The no-skill average precision is the positive rate."
        ),
        "split": lib_split.describe(genome_gc),
        "n_train_rows": len(train_rows),
        "n_test_rows": len(test_rows),
        "train_positive_rate": round(float(y_train.mean()), 5),
        "test_positive_rate": round(float(y_test.mean()), 5),
        "no_skill_average_precision": round(float(y_test.mean()), 5),
        "baselines": results,
    }
    METRICS.write_text(json.dumps(payload, indent=2) + "\n")

    print(f"train {len(train_rows):,} rows / {len(train_genomes)} genomes "
          f"(positive rate {y_train.mean():.4f})")
    print(f"test  {len(test_rows):,} rows / {len(test_genomes)} genomes "
          f"(positive rate {y_test.mean():.4f})")
    print(f"no-skill average precision = {y_test.mean():.4f}")
    for r in results:
        print(f"  {r['baseline']:<18} AP {r['average_precision']:.4f}  "
              f"ROC-AUC {r['roc_auc'] if r['roc_auc'] is not None else 'n/a'}")
    print(f"wrote {METRICS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

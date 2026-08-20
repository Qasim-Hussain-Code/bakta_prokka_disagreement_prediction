#!/usr/bin/env python3
"""Fitting procedure shared by 08 (all features) and 10 (caller group removed).

Step 10's claim is that dropping the gene-caller-derived features changes the
score. That claim only means something if everything else is held identical:
same folds, same depth grid, same selection rule, same seed. Putting the
procedure here rather than writing it twice is what makes the two numbers
comparable.
"""

import gzip
import json
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score

import lib_split

ROOT = Path(__file__).resolve().parent.parent
FEATURES_TSV = ROOT / "data" / "interim" / "features.tsv.gz"
FEATURES_JSON = ROOT / "results" / "metrics" / "06_features.json"

DEPTH_GRID = [4, 6, 8, 12, 16, None]
N_ESTIMATORS = 300
RANDOM_STATE = 0


def load():
    """-> (rows, idx, manifest). Rows stay as strings; callers pick columns."""
    if not FEATURES_TSV.exists():
        raise SystemExit(f"{FEATURES_TSV} absent; run scripts/06_features.py first")
    with gzip.open(FEATURES_TSV, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh]
    manifest = json.loads(FEATURES_JSON.read_text())
    return rows, {c: i for i, c in enumerate(header)}, manifest


def feature_names(manifest, exclude_groups=()):
    groups = manifest["features_by_group"]
    return [
        name
        for group, names in sorted(groups.items())
        if group not in exclude_groups
        for name in names
    ]


def matrix(rows, idx, names):
    return np.array([[float(r[idx[n]]) for n in names] for r in rows],
                    dtype=np.float64)


def genome_gc_map(rows, idx):
    out = {}
    for r in rows:
        out.setdefault(r[idx["genome"]], float(r[idx["genome_gc"]]))
    return out


def new_forest(max_depth):
    return RandomForestClassifier(
        n_estimators=N_ESTIMATORS,
        max_depth=max_depth,
        # The label is heavily imbalanced. Without reweighting the forest can
        # reach a near-perfect accuracy by never predicting a disagreement.
        class_weight="balanced_subsample",
        n_jobs=-1,
        random_state=RANDOM_STATE,
    )


def sweep_depth(X, y, fold_of_row):
    """Select max_depth on the training folds only.

    The test genomes are not touched here. That is the whole point of doing
    the sweep in a separate function: the number reported in step 08 is from
    a single evaluation on genomes that played no part in choosing the depth.
    """
    folds = sorted(set(fold_of_row))
    results = []
    for depth in DEPTH_GRID:
        scores = []
        for held in folds:
            tr = fold_of_row != held
            te = fold_of_row == held
            if len(set(y[tr])) < 2 or not te.any():
                continue
            model = new_forest(depth).fit(X[tr], y[tr])
            proba = model.predict_proba(X[te])[:, 1]
            scores.append(average_precision_score(y[te], proba))
        results.append({
            "max_depth": depth,
            "mean_average_precision": round(float(np.mean(scores)), 5) if scores else None,
            "per_fold": [round(float(s), 5) for s in scores],
        })
    scored = [r for r in results if r["mean_average_precision"] is not None]
    if not scored:
        raise SystemExit("depth sweep produced no usable folds")
    best = max(scored, key=lambda r: r["mean_average_precision"])
    return best["max_depth"], results


def evaluate(model, X, y, genomes):
    """Test-set metrics, overall and per genome."""
    proba = model.predict_proba(X)[:, 1]
    overall = {
        "n_rows": int(len(y)),
        "n_positive": int(y.sum()),
        "positive_rate": round(float(y.mean()), 5),
        "no_skill_average_precision": round(float(y.mean()), 5),
        "average_precision": round(float(average_precision_score(y, proba)), 5),
        "roc_auc": (round(float(roc_auc_score(y, proba)), 5)
                    if len(set(y)) > 1 else None),
    }
    per_genome = []
    for g in sorted(set(genomes)):
        m = genomes == g
        yg = y[m]
        per_genome.append({
            "genome": g,
            "n_rows": int(m.sum()),
            "n_positive": int(yg.sum()),
            "positive_rate": round(float(yg.mean()), 5),
            "average_precision": (round(float(average_precision_score(yg, proba[m])), 5)
                                  if len(set(yg)) > 1 else None),
            "roc_auc": (round(float(roc_auc_score(yg, proba[m])), 5)
                        if len(set(yg)) > 1 else None),
        })
    return overall, per_genome, proba


def prepare(exclude_groups=()):
    """Everything 08 and 10 need, assembled the same way for both."""
    rows, idx, manifest = load()
    names = feature_names(manifest, exclude_groups)
    gc = genome_gc_map(rows, idx)
    train_genomes, test_genomes, folds = lib_split.split_genomes(gc)
    train_set = set(train_genomes)

    train_rows = [r for r in rows if r[idx["genome"]] in train_set]
    test_rows = [r for r in rows if r[idx["genome"]] not in train_set]

    return {
        "names": names,
        "manifest": manifest,
        "split": lib_split.describe(gc),
        "X_train": matrix(train_rows, idx, names),
        "y_train": np.array([int(r[idx["label"]]) for r in train_rows]),
        "fold_train": np.array([folds[r[idx["genome"]]] for r in train_rows]),
        "X_test": matrix(test_rows, idx, names),
        "y_test": np.array([int(r[idx["label"]]) for r in test_rows]),
        "genomes_test": np.array([r[idx["genome"]] for r in test_rows]),
        "train_genomes": train_genomes,
        "test_genomes": test_genomes,
    }

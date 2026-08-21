#!/usr/bin/env python3
"""Fitting and scoring procedure shared by 15-19.

Part one's claim in step 10 was that dropping the caller features changes the
score. That claim only means something if everything else is held identical.
The same applies here across three arms and three model families, so the
folds, the grids, the selection rule and the seed live in one place.

THE SELECTION RULE, declared before the first fit:

    mean Matthews correlation coefficient across five genome-grouped folds of
    the training genomes, then the one-standard-error rule: among all
    candidates whose mean MCC is within one standard error of the best, take
    the simplest.

A note on provenance, because the brief described this as part one's rule
reused verbatim and it is not. scripts/lib_model.py selects max_depth by mean
AVERAGE PRECISION with a plain argmax and no one-standard-error step. The rule
implemented here is the one specified for this half, and it is better suited
to it: the naming label is close to balanced (51.7% positive) where part one's
was 0.08%, and MCC is the headline metric requested. The discrepancy is
recorded in every sweep file rather than resolved silently in either
direction. Scripts 01-11 are untouched.

THE SPLIT: fold 0 of lib_split is the test set, exactly as in part one, so
both halves of the project are scored on the same five held-out genomes. The
remaining 20 training genomes are dealt into five cross-validation folds by
the same deterministic GC-rank round-robin. No genome appears in both training
and test, and that is asserted rather than assumed.
"""

import gzip
import json
from pathlib import Path

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import (accuracy_score, confusion_matrix, f1_score,
                             matthews_corrcoef, precision_score, recall_score,
                             roc_auc_score)
from sklearn.tree import DecisionTreeClassifier

import lib_split
from lib_model import MIN_POSITIVES_FOR_STABLE_AP

ROOT = Path(__file__).resolve().parent.parent
FEATURES_TSV = ROOT / "data" / "interim" / "content_features.tsv.gz"
FEATURES_JSON = ROOT / "results" / "metrics" / "14_content_features.json"

RANDOM_STATE = 42
N_CV_FOLDS = 5

# --- grids, declared before the first fit ---------------------------------
# Each entry is (label, kwargs). `complexity` below defines what "simplest"
# means for the one-standard-error tie-break, ascending.
GRIDS = {
    "decision_tree": [
        {"max_depth": d, "min_samples_leaf": leaf}
        for d in (2, 4, 6, 8, 12, 16, None)
        for leaf in (1, 20)
    ],
    "random_forest": [
        {"max_depth": d, "min_samples_leaf": leaf}
        for d in (4, 8, 12, 16, None)
        for leaf in (1, 5)
    ],
    "gradient_boosting": [
        {"max_depth": d, "learning_rate": lr, "max_iter": n}
        for d in (2, 3, 4)
        for lr in (0.05, 0.1)
        for n in (200,)
    ],
}

# Gradient boosting uses HistGradientBoostingClassifier rather than
# GradientBoostingClassifier. The change was made before any gradient-boosting
# result existed, for tractability alone, and is recorded rather than left
# implicit: on this training set a single GradientBoostingClassifier fit takes
# 54.1 s against 0.7 s for the histogram implementation, so one sweep would be
# 27 minutes against 0.4. The grid above is the same grid in the histogram
# implementation's parameter names -- max_iter is its n_estimators.
GRADIENT_BOOSTING_IMPLEMENTATION = {
    "class": "sklearn.ensemble.HistGradientBoostingClassifier",
    "instead_of": "sklearn.ensemble.GradientBoostingClassifier",
    "reason": ("computational tractability; measured 54.1 s versus 0.7 s per "
               "fit on this training set, a 75x difference"),
    "decided": "before any gradient-boosting result was computed",
}

N_ESTIMATORS_FOREST = 300


def complexity(params):
    """Ascending = simpler. Used only for the one-standard-error tie-break."""
    depth = params.get("max_depth")
    depth_rank = 999 if depth is None else depth
    return (
        depth_rank,
        params.get("n_estimators", params.get("max_iter", 0)),
        params.get("learning_rate", 0.0),
        -params.get("min_samples_leaf", 1),
    )


def build(model, params, oob=False):
    if model == "decision_tree":
        return DecisionTreeClassifier(
            class_weight="balanced", random_state=RANDOM_STATE, **params)
    if model == "random_forest":
        return RandomForestClassifier(
            n_estimators=N_ESTIMATORS_FOREST, class_weight="balanced_subsample",
            n_jobs=-1, random_state=RANDOM_STATE, oob_score=oob,
            bootstrap=True, **params)
    if model == "gradient_boosting":
        return HistGradientBoostingClassifier(
            class_weight="balanced", random_state=RANDOM_STATE, **params)
    raise ValueError(f"unknown model {model!r}")


# --- data -----------------------------------------------------------------

def load():
    if not FEATURES_TSV.exists():
        raise SystemExit(
            f"{FEATURES_TSV} absent; run scripts/14_content_features.py first")
    with gzip.open(FEATURES_TSV, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        rows = [line.rstrip("\n").split("\t") for line in fh]
    manifest = json.loads(FEATURES_JSON.read_text())
    # The primary analysis set is the only place the label is defined.
    rows = [r for r in rows if r[idx["in_primary_set"]] == "1"]
    if not rows:
        raise SystemExit("FATAL: primary analysis set is empty")
    return rows, idx, manifest


def feature_names(manifest, drop_caller=False, drop_db=False):
    """Read the flags from the manifest. Never hard-code a feature list."""
    out = []
    for name, meta in manifest["manifest"].items():
        if drop_caller and meta["caller_derived"]:
            continue
        if drop_db and meta["db_derived"]:
            continue
        out.append(name)
    return sorted(out)


def cv_folds(train_genomes, genome_gc):
    """Deal the training genomes into N_CV_FOLDS by GC rank, round-robin.

    Same algorithm as lib_split.assign_folds, applied to the training genomes
    only, so every fold spans the GC range instead of a random draw handing
    one fold all the high-GC genomes.
    """
    ordered = sorted(train_genomes, key=lambda g: (genome_gc[g], g))
    return {g: i % N_CV_FOLDS for i, g in enumerate(ordered)}


def prepare(drop_caller=False, drop_db=False):
    rows, idx, manifest = load()
    names = feature_names(manifest, drop_caller, drop_db)
    if not names:
        raise SystemExit("FATAL: no features left after applying the flags")

    gc = {}
    for r in rows:
        gc.setdefault(r[idx["genome"]], float(r[idx["genome_gc"]]))
    train_genomes, test_genomes, _ = lib_split.split_genomes(gc)

    overlap = sorted(set(train_genomes) & set(test_genomes))
    if overlap:
        raise SystemExit(
            f"FATAL: {len(overlap)} genome(s) in both training and test: {overlap}")

    folds = cv_folds(train_genomes, gc)
    train_set = set(train_genomes)
    tr = [r for r in rows if r[idx["genome"]] in train_set]
    te = [r for r in rows if r[idx["genome"]] not in train_set]

    def matrix(subset):
        return np.array([[float(r[idx[n]]) for n in names] for r in subset],
                        dtype=np.float64)

    def labels(subset):
        return np.array([int(r[idx["label"]]) for r in subset])

    return {
        "names": names,
        "manifest": manifest,
        "split": lib_split.describe(gc),
        "split_assertion": {
            "train_genomes": train_genomes,
            "test_genomes": test_genomes,
            "n_overlapping_genomes": 0,
            "assertion": "set(train_genomes) & set(test_genomes) == empty",
            "result": "passed",
        },
        "cv_folds": {"n_folds": N_CV_FOLDS,
                     "assignment": {g: folds[g] for g in sorted(folds)}},
        "X_train": matrix(tr), "y_train": labels(tr),
        "fold_train": np.array([folds[r[idx["genome"]]] for r in tr]),
        "genomes_train": np.array([r[idx["genome"]] for r in tr]),
        "X_test": matrix(te), "y_test": labels(te),
        "genomes_test": np.array([r[idx["genome"]] for r in te]),
        "train_genomes": train_genomes, "test_genomes": test_genomes,
    }


# --- selection ------------------------------------------------------------

def sweep(model, X, y, fold_of_row, grid=None):
    """Mean MCC over the grouped CV folds, then the one-standard-error rule.

    The test genomes are not touched. Selection happens entirely inside the
    training folds, which is what makes the single test evaluation meaningful.
    """
    grid = grid if grid is not None else GRIDS[model]
    folds = sorted(set(fold_of_row.tolist()))
    results = []
    for params in grid:
        scores = []
        for held in folds:
            tr, te = fold_of_row != held, fold_of_row == held
            if len(set(y[tr])) < 2 or not te.any():
                continue
            fitted = build(model, params).fit(X[tr], y[tr])
            scores.append(matthews_corrcoef(y[te], fitted.predict(X[te])))
        if not scores:
            continue
        mean = float(np.mean(scores))
        # standard error of the mean across folds
        se = float(np.std(scores, ddof=1) / np.sqrt(len(scores))) if len(scores) > 1 else 0.0
        results.append({
            "params": dict(params),
            "mean_mcc": round(mean, 5),
            "std_error": round(se, 5),
            "per_fold_mcc": [round(float(s), 5) for s in scores],
        })
    if not results:
        raise SystemExit(f"FATAL: {model} sweep produced no usable folds")

    best = max(results, key=lambda r: r["mean_mcc"])
    threshold = best["mean_mcc"] - best["std_error"]
    within = [r for r in results if r["mean_mcc"] >= threshold]
    chosen = min(within, key=lambda r: complexity(r["params"]))

    selection = {
        "rule": ("mean MCC across five genome-grouped CV folds, then the "
                 "one-standard-error rule: among candidates within one "
                 "standard error of the best mean, take the simplest"),
        "grid_declared_before_first_fit": True,
        "n_candidates": len(results),
        "best_mean_mcc": best["mean_mcc"],
        "best_params": best["params"],
        "one_se_threshold": round(threshold, 5),
        "n_within_one_se": len(within),
        "chosen_params": chosen["params"],
        "chosen_mean_mcc": chosen["mean_mcc"],
        "chosen_is_simpler_than_best": chosen["params"] != best["params"],
        "implementation_note": (
            GRADIENT_BOOSTING_IMPLEMENTATION if model == "gradient_boosting"
            else None),
        "differs_from_part_one_rule": (
            "scripts/lib_model.py selects by mean average precision with a "
            "plain argmax and no one-standard-error step. This half uses MCC "
            "with the one-standard-error rule, as specified for it. Recorded "
            "rather than resolved silently in either direction."),
    }
    return chosen["params"], results, selection


# --- scoring --------------------------------------------------------------

def evaluate(model, X, y, genomes, proba=None):
    """Test-set metrics. Counts, not just rates.

    The reliability warning from part one's lib_model.evaluate applies here
    unchanged: if the held-out set carries fewer than 30 positives the caveat
    is attached to the number rather than left in a notebook.
    """
    pred = model.predict(X)
    if proba is None:
        proba = (model.predict_proba(X)[:, 1]
                 if hasattr(model, "predict_proba") else pred.astype(float))
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    out = {
        "n_rows": int(len(y)),
        "n_positive": int(y.sum()),
        "n_negative": int((y == 0).sum()),
        "positive_rate": round(float(y.mean()), 5),
        "accuracy": round(float(accuracy_score(y, pred)), 5),
        "precision": round(float(precision_score(y, pred, zero_division=0)), 5),
        "recall": round(float(recall_score(y, pred, zero_division=0)), 5),
        "f1": round(float(f1_score(y, pred, zero_division=0)), 5),
        "mcc": round(float(matthews_corrcoef(y, pred)), 5),
        "roc_auc": (round(float(roc_auc_score(y, proba)), 5)
                    if len(set(y.tolist())) > 1 else None),
        "confusion_matrix_counts": {
            "true_negative": int(tn), "false_positive": int(fp),
            "false_negative": int(fn), "true_positive": int(tp),
        },
        "majority_class_accuracy": round(float(max(y.mean(), 1 - y.mean())), 5),
    }
    if int(y.sum()) < MIN_POSITIVES_FOR_STABLE_AP:
        out["reliability_warning"] = (
            f"only {int(y.sum())} positives in the held-out set. Metrics on "
            "this few positives carry very wide uncertainty -- a single "
            "differently-classified positive moves them substantially. Do not "
            "quote these as point estimates."
        )
    per_genome = []
    for g in sorted(set(genomes.tolist())):
        m = genomes == g
        yg, pg = y[m], pred[m]
        per_genome.append({
            "genome": g,
            "n_rows": int(m.sum()),
            "n_positive": int(yg.sum()),
            "positive_rate": round(float(yg.mean()), 5),
            "accuracy": round(float(accuracy_score(yg, pg)), 5),
            "mcc": (round(float(matthews_corrcoef(yg, pg)), 5)
                    if len(set(yg.tolist())) > 1 else None),
        })
    return out, per_genome, proba


def oob_report(forest, y_train):
    """OOB alongside grouped CV, with the interpretation stated in the file."""
    if not hasattr(forest, "oob_decision_function_"):
        return None
    oob_proba = forest.oob_decision_function_[:, 1]
    oob_pred = (oob_proba >= 0.5).astype(int)
    return {
        "oob_accuracy": round(float(forest.oob_score_), 5),
        "oob_error": round(1.0 - float(forest.oob_score_), 5),
        "oob_mcc": round(float(matthews_corrcoef(y_train, oob_pred)), 5),
        "interpretation": (
            "OOB bootstraps rows, not genomes. Rows from the same genome land "
            "both in-bag and out-of-bag, so an OOB estimate is scored against "
            "near-neighbours from the same annotation run. The gap between OOB "
            "and the genome-grouped CV mean is therefore a measurement of "
            "genome-level leakage, not a discrepancy to be explained away. OOB "
            "above grouped CV is the expected direction and its size is the "
            "quantity of interest."),
    }

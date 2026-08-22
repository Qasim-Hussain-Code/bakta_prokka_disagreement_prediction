#!/usr/bin/env python3
"""Every headline number for the content experiment, in one file.

Reads:  results/metrics/12..19_*.json
Writes: results/metrics/20_content_summary.json

This exists so that prose can be checked against a single file rather than
against nine. Every value here is copied from another metrics file at run
time; nothing is computed here and nothing is typed in. If a number in the
README is not in this file, it was not verified.
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
METRICS = ROOT / "results" / "metrics"
OUT = METRICS / "20_content_summary.json"

MODELS = ("decision_tree", "random_forest", "gradient_boosting")


def load(name, required=True):
    path = METRICS / name
    if not path.exists():
        if required:
            raise SystemExit(f"{path} absent; run the earlier steps first")
        return None
    return json.loads(path.read_text())


def main():
    cohort = load("12_content_cohort.json")
    rules = load("13_name_rules.json")
    feats = load("14_content_features.json")
    base = load("15_content_baselines.json")
    circ = load("19_content_circularity.json")
    mech = load("22_content_mechanism.json", required=False)

    finals, cvs, imps = {}, {}, {}
    for m in MODELS:
        finals[m] = load(f"17_content_final_{m}.json", required=False)
        cvs[m] = load(f"16_content_cv_{m}.json", required=False)
        imps[m] = load(f"18_content_importance_{m}.json", required=False)
    if not any(finals.values()):
        raise SystemExit("no 17_content_final_*.json files; run step 17 first")

    prim = cohort["primary_analysis_set"]
    ct = cohort["cross_tabulation"]

    model_table = {}
    for m, f in finals.items():
        if f is None:
            continue
        t = f["test"]
        model_table[m] = {
            "hyperparameters": f["hyperparameters"],
            "grouped_cv_mean_mcc": f["grouped_cv_mean_mcc"],
            "test_accuracy": t["accuracy"],
            "test_precision": t["precision"],
            "test_recall": t["recall"],
            "test_f1": t["f1"],
            "test_mcc": t["mcc"],
            "test_roc_auc": t["roc_auc"],
            "confusion_matrix_counts": t["confusion_matrix_counts"],
            "reliability_warning": t.get("reliability_warning"),
        }
        if "oob" in f:
            model_table[m]["oob_mcc"] = f["oob"]["oob_mcc"]
            model_table[m]["oob_error"] = f["oob"]["oob_error"]
            model_table[m]["oob_minus_grouped_cv_mcc"] = f["oob_vs_grouped_cv"]["gap"]

    baseline_table = {}
    for name, b in base["baselines"].items():
        if not b.get("computed", True):
            baseline_table[name] = {"computed": False, "reason": b["reason"]}
            continue
        t = b["test"]
        baseline_table[name] = {
            "computed": True,
            "test_accuracy": t["accuracy"], "test_f1": t["f1"],
            "test_mcc": t["mcc"],
        }

    best_model = max(
        (m for m in model_table),
        key=lambda m: model_table[m]["test_mcc"])
    best_baseline = max(
        (n for n, b in baseline_table.items() if b.get("computed")),
        key=lambda n: baseline_table[n]["test_mcc"])

    top_perm = {}
    for m, imp in imps.items():
        if imp is None:
            continue
        top_perm[m] = [
            {"feature": r["feature"], "caller_derived": r["caller_derived"],
             "db_derived": r["db_derived"],
             "permutation_importance_mean": r["permutation_importance_mean"]}
            for r in imp["ranked_by_permutation"][:8]
        ]

    payload = {
        "step": "20_content_summary",
        "experiment": (
            "Given a region both Bakta and Prokka called as a CDS, do the two "
            "tools give it the same product name, and is that predictable?"),
        "what_the_label_is": rules["what_the_label_is"],
        "seed": 42,

        "cohort": {
            "bakta_cds_total": cohort["bakta_cds_total"],
            "prokka_cds_total": cohort["prokka_cds_total"],
            "n_paired": cohort["n_paired"],
            "n_unpaired": cohort["n_unpaired"],
            "cross_tabulation": ct,
            "primary_analysis_set_n": prim["n"],
            "primary_analysis_set_share": prim["share_of_pairs"],
            "why_conditional": prim["why_conditional"],
        },
        "db_light_confound": cohort["asymmetric_cells_interpretation"],

        "name_disagreement": {
            "primary_rule": prim["name_disagreement"]["primary"],
            "strict_sensitivity_check": prim["name_disagreement"]["strict"],
            "loose_sensitivity_check": prim["name_disagreement"]["loose"],
            "modelled": "primary",
            "interpretation": rules["interpretation"],
        },
        "separate_columns": cohort["separate_columns_not_part_of_the_label"],
        "headline_pair_with_and_without_fallback_families": (
            {
                "whole_primary_set": mech["fallback_naming_families"]["whole_primary_set"],
                "excluding_fallback_families":
                    mech["fallback_naming_families"]["remainder"],
                "fallback_families":
                    mech["fallback_naming_families"]["families"],
                "statement": mech["fallback_naming_families"]["headline_pair"],
            } if mech else
            {"missing": "run scripts/22_content_mechanism.py before this step"}),
        "mechanism": ({
            "bakta_evidence_depth": mech["bakta_evidence_depth_remainder_only"],
            "prokka_evidence_depth": mech["prokka_evidence_depth_remainder_only"],
            "conclusion": mech["conclusion"],
        } if mech else None),
        "declared_generic_families": cohort[
            "declared_generic_families_within_primary_set"]["families"],

        "features": {
            "n_features": feats["n_features"],
            "counts": feats["counts"],
            "excluded_as_constant": feats.get("excluded_as_constant_in_primary_set"),
            "near_constant": feats["constant_feature_check"]["near_constant"],
        },

        "split": {
            "strategy": "grouped by genome; fold 0 of lib_split is the test set",
            "test_genomes": base["test_genomes"],
            "n_test_rows": base["n_test_rows"],
            "n_test_positive": base["n_test_positive"],
            "test_positive_rate": base["test_positive_rate"],
            "assertion": base["split"]["assertion"],
            "assertion_result": base["split"]["result"],
            "n_overlapping_genomes": base["split"]["n_overlapping_genomes"],
        },

        "baselines_test_set": baseline_table,
        "models_test_set": model_table,
        "best_model_by_test_mcc": best_model,
        "best_baseline_by_test_mcc": best_baseline,
        "best_model_minus_best_baseline_mcc": round(
            model_table[best_model]["test_mcc"]
            - baseline_table[best_baseline]["test_mcc"], 5),

        "circularity_audit": {
            "arms": circ["side_by_side"],
            "gaps_in_test_mcc": circ["gaps_in_test_mcc"],
            "reference_points": circ["reference_points_same_test_rows"],
        },

        "permutation_importance_top_8": top_perm,

        "hyperparameter_selection_rule": (
            "mean MCC across five genome-grouped folds of the training "
            "genomes, then the one-standard-error rule: among candidates "
            "within one standard error of the best mean, take the simplest. "
            "Grids declared before the first fit."),
        "selection_rule_provenance": (
            cvs[best_model]["selection"]["differs_from_part_one_rule"]
            if cvs.get(best_model) else None),

        "source_files": {
            "cohort": "12_content_cohort.json",
            "rules": "13_name_rules.json",
            "features": "14_content_features.json",
            "baselines": "15_content_baselines.json",
            "cv": [f"16_content_cv_{m}.json" for m in MODELS if cvs.get(m)],
            "final": [f"17_content_final_{m}.json" for m in MODELS if finals.get(m)],
            "importance": [f"18_content_importance_{m}.json" for m in MODELS if imps.get(m)],
            "circularity": "19_content_circularity.json",
            "mechanism": "22_content_mechanism.json" if mech else None,
        },
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n")

    print(f"primary analysis set {prim['n']:,}  "
          f"name_disagreement {prim['name_disagreement']['primary']['rate']:.1%}")
    print(f"test set {base['n_test_rows']:,} rows, "
          f"{base['n_test_positive']:,} positive")
    print("\nbaselines (test MCC):")
    for n, b in baseline_table.items():
        print(f"  {n:34s} {b['test_mcc']:+.4f}" if b.get("computed")
              else f"  {n:34s} not computable")
    print("models (test MCC):")
    for m, t in model_table.items():
        print(f"  {m:34s} {t['test_mcc']:+.4f}")
    print("circularity arms (test MCC):")
    for arm, a in circ["side_by_side"].items():
        print(f"  {arm:34s} {a['test_mcc']:+.4f}  ({a['n_features']} features)")
    print(f"\nwrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

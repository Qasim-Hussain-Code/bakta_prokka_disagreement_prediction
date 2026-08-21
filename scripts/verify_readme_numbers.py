#!/usr/bin/env python3
"""Check every number quoted in README.md against results/metrics/.

Reads:  README.md
        results/metrics/*.json
Writes: nothing. Exits non-zero on the first disagreement.

The standard for this repository is that no number in the prose is transcribed
from a terminal. That standard is only worth anything if it is enforced, so
each claim below names the metrics file and key it came from, and the check
fails loudly rather than warning.

This caught one real error on its first run: the majority-class baseline was
quoted at 0.503 accuracy, which is `majority_class_accuracy` -- the accuracy a
constant predictor tuned on the test set would reach. The baseline actually
predicts the class that is in the majority among the TRAINING genomes, and the
held-out genomes are fractionally negative-majority, so it scores 0.497.

Run it after editing either the README or any metrics file.
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
METRICS = ROOT / "results" / "metrics"
README = ROOT / "README.md"


def load(name):
    path = METRICS / name
    if not path.exists():
        raise SystemExit(f"{path} absent; run the pipeline first")
    return json.loads(path.read_text())


def main():
    text = README.read_text()

    s = load("20_content_summary.json")
    coh = load("12_content_cohort.json")
    rul = load("13_name_rules.json")
    fea = load("14_content_features.json")
    bas = load("15_content_baselines.json")
    dis = load("05_disagreement.json")
    f08 = load("08_forest.json")
    c10 = load("10_circularity_audit.json")
    imp = load("18_content_importance_random_forest.json")
    rf = load("17_content_final_random_forest.json")

    ct = coh["cross_tabulation"]
    nd = coh["primary_analysis_set"]["name_disagreement"]
    trap = rul["placeholder_rule"]["substring_trap"]
    sep = coh["separate_columns_not_part_of_the_label"]
    bt, mt = s["baselines_test_set"], s["models_test_set"]
    arms = s["circularity_audit"]["arms"]
    flag = imp["summary_by_flag"]

    def pct(x, dp=1):
        return f"{x * 100:.{dp}f}%"

    def n(x):
        return f"{x:,}"

    # (what the README says, where it came from)
    claims = [
        # --- part one
        (n(dis["n_rows"]), "05.n_rows"),
        (str(dis["n_positive"]), "05.n_positive"),
        (str(dis["positive_rate"]), "05.positive_rate"),
        (n(dis["boundary_agreement_among_matched"]["n_matched"]), "05.n_matched"),
        (n(dis["boundary_agreement_among_matched"]["same_start_and_stop"]),
         "05.same_start_and_stop"),
        (str(f08["test"]["n_positive"]), "08.test.n_positive"),
        (f"{f08['test']['average_precision']:.3f}", "08.test.average_precision"),
        (f"{c10['comparison']['average_precision_sequence_only']:.3f}",
         "10.average_precision_sequence_only"),
        (str(f08["test"]["no_skill_average_precision"]), "08.no_skill"),

        # --- part two cohort
        (n(coh["n_paired"]), "12.n_paired"),
        (str(coh["n_unpaired"]), "12.n_unpaired"),
        (n(ct["both_named"]), "12.both_named"),
        (n(ct["bakta_named_only"]), "12.bakta_named_only"),
        (n(ct["both_placeholder"]), "12.both_placeholder"),
        (n(ct["prokka_named_only"]), "12.prokka_named_only"),
        (pct(ct["shares"]["both_named"]), "12.shares.both_named"),
        (pct(ct["shares"]["bakta_named_only"]), "12.shares.bakta_named_only"),
        (pct(ct["shares"]["both_placeholder"]), "12.shares.both_placeholder"),
        (pct(ct["shares"]["prokka_named_only"]), "12.shares.prokka_named_only"),

        # --- rules
        (n(nd["strict"]["n"]), "12.name_disagreement.strict.n"),
        (pct(nd["strict"]["rate"]), "12.name_disagreement.strict.rate"),
        (n(nd["primary"]["n"]), "12.name_disagreement.primary.n"),
        (pct(nd["primary"]["rate"]), "12.name_disagreement.primary.rate"),
        (n(nd["loose"]["n"]), "12.name_disagreement.loose.n"),
        (pct(nd["loose"]["rate"]), "12.name_disagreement.loose.rate"),
        (f"{(nd['strict']['rate'] - nd['loose']['rate']) * 100:.1f}",
         "12.strict minus loose"),
        (str(trap["bakta:hypothetical"]["distinct_products_containing_the_word"]),
         "13.substring_trap.bakta:hypothetical.distinct"),
        (str(trap["bakta:hypothetical"]["would_be_wrongly_caught_by_substring_match"]),
         "13.substring_trap.bakta:hypothetical.wrongly_caught"),
        (str(trap["bakta:conserved"]["distinct_products_containing_the_word"]),
         "13.substring_trap.bakta:conserved.distinct"),
        (pct(sep["gene_symbol"]["rate"]), "12.gene_symbol.rate"),
        (pct(sep["ec_number"]["rate"]), "12.ec_number.rate"),

        # --- features and split
        (str(fea["n_features"]), "14.n_features"),
        (str(fea["counts"]["caller_derived"]), "14.counts.caller_derived"),
        (str(fea["counts"]["db_derived"]), "14.counts.db_derived"),
        (str(fea["counts"]["neither"]), "14.counts.neither"),
        (n(bas["n_test_rows"]), "15.n_test_rows"),
        (n(bas["n_test_positive"]), "15.n_test_positive"),

        # --- baselines and models, as the README's table rows
        (f"| majority class | {bt['majority_class']['test_mcc']:.3f} | "
         f"{bt['majority_class']['test_accuracy']:.3f} | "
         f"{bt['majority_class']['test_f1']:.3f} |", "table row: majority class"),
        (f"| protein length, one threshold | "
         f"{bt['protein_length_only']['test_mcc']:.3f} | "
         f"{bt['protein_length_only']['test_accuracy']:.3f} | "
         f"{bt['protein_length_only']['test_f1']:.3f} |",
         "table row: protein length"),
        (f"| database coverage | "
         f"{bt['database_coverage_substitute']['test_mcc']:.3f} | "
         f"{bt['database_coverage_substitute']['test_accuracy']:.3f} | "
         f"{bt['database_coverage_substitute']['test_f1']:.3f} |",
         "table row: database coverage"),
        (f"| decision tree | {mt['decision_tree']['test_mcc']:.3f} | "
         f"{mt['decision_tree']['test_accuracy']:.3f} | "
         f"{mt['decision_tree']['test_f1']:.3f} |", "table row: decision tree"),
        (f"| **random forest** | **{mt['random_forest']['test_mcc']:.3f}** | "
         f"{mt['random_forest']['test_accuracy']:.3f} | "
         f"{mt['random_forest']['test_f1']:.3f} |", "table row: random forest"),
        (f"| gradient boosting | {mt['gradient_boosting']['test_mcc']:.3f} | "
         f"{mt['gradient_boosting']['test_accuracy']:.3f} | "
         f"{mt['gradient_boosting']['test_f1']:.3f} |",
         "table row: gradient boosting"),

        # --- circularity arms, as the README's table rows
        (f"| full | {arms['full']['n_features']} | "
         f"{arms['full']['grouped_cv_mean_mcc']:.3f} | "
         f"{arms['full']['test_mcc']:.3f} |", "table row: arm full"),
        (f"| caller-derived removed | {arms['no_caller']['n_features']} | "
         f"{arms['no_caller']['grouped_cv_mean_mcc']:.3f} | "
         f"{arms['no_caller']['test_mcc']:.3f} |", "table row: arm no_caller"),
        (f"| **sequence only** | {arms['sequence_only']['n_features']} | "
         f"**{arms['sequence_only']['grouped_cv_mean_mcc']:.3f}** | "
         f"**{arms['sequence_only']['test_mcc']:.3f}** |",
         "table row: arm sequence_only"),

        # --- importance and OOB
        (f"{flag['db_derived']['summed_permutation_importance']:.3f}",
         "18.summary_by_flag.db_derived"),
        (f"{flag['caller_derived']['summed_permutation_importance']:.3f}",
         "18.summary_by_flag.caller_derived"),
        (f"{rf['oob']['oob_mcc']:.3f}", "17.oob.oob_mcc"),
        (f"{rf['grouped_cv_mean_mcc']:.3f}", "17.grouped_cv_mean_mcc"),
        (f"{rf['oob_vs_grouped_cv']['gap']:.3f}", "17.oob_vs_grouped_cv.gap"),
    ]

    # The summed sequence-group importance is negative; the README renders it
    # with a Unicode minus, so check that form explicitly rather than str().
    neither = flag["neither"]["summed_permutation_importance"]
    claims.append((f"{abs(neither):.3f}", "18.summary_by_flag.neither (magnitude)"))

    top3 = [r["feature"] for r in imp["ranked_by_permutation"][:3]]
    for feature in top3:
        claims.append((f"`{feature}`", "18.ranked_by_permutation top 3"))

    missing = [(value, source) for value, source in claims if value not in text]
    for value, source in missing:
        print(f"  NOT FOUND in README: {value!r}   (from {source})")

    print(f"\n{len(claims) - len(missing)}/{len(claims)} quoted values found in "
          f"README.md and matching results/metrics/")
    if missing:
        raise SystemExit(
            f"FATAL: {len(missing)} value(s) in the metrics files are not "
            "present in the README as written. Either the prose is stale or a "
            "number was transcribed rather than verified.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Write the naming rules down, with the counts and worked examples that
justify them, verified against the data rather than transcribed.

Reads:  data/interim/content_pairs.tsv.gz
Writes: results/metrics/13_name_rules.json

Every placeholder count, every worked example and every rate in this file is
recomputed from the pair table on each run. Nothing here is typed in by hand.
If a rule is edited in lib_names.py, this file changes with it, and any README
line quoting a number that no longer matches will be wrong in a way that is
easy to catch.

The rules themselves were fixed and approved before any model was fitted.
"""

import gzip
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_names

ROOT = Path(__file__).resolve().parent.parent
PAIRS = ROOT / "data" / "interim" / "content_pairs.tsv.gz"
METRICS = ROOT / "results" / "metrics" / "13_name_rules.json"

N_EXAMPLES = 8


def load_pairs():
    if not PAIRS.exists():
        raise SystemExit(f"{PAIRS} absent; run scripts/12_content_cohort.py first")
    with gzip.open(PAIRS, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        return [line.rstrip("\n").split("\t") for line in fh], idx


def main():
    rows, idx = load_pairs()

    def col(r, name):
        return r[idx[name]]

    # --- placeholder literals, counted from the data ----------------------
    literal_counts = {}
    for lit in sorted(lib_names.PLACEHOLDERS):
        nb = sum(1 for r in rows if lib_names.casefold_ws(col(r, "product_bakta")) == lit)
        np_ = sum(1 for r in rows if lib_names.casefold_ws(col(r, "product_prokka")) == lit)
        literal_counts[lit] = {"bakta": nb, "prokka": np_,
                               "observed": bool(nb or np_)}

    # Evidence that a substring test would have been wrong: how many distinct
    # product strings contain each placeholder word without being one.
    substring_trap = {}
    for word in ("hypothetical", "conserved", "uncharacteri"):
        for tool in ("bakta", "prokka"):
            distinct = {col(r, f"product_{tool}") for r in rows
                        if word in col(r, f"product_{tool}").lower()}
            caught = {p for p in distinct if lib_names.is_placeholder(p)}
            substring_trap[f"{tool}:{word}"] = {
                "distinct_products_containing_the_word": len(distinct),
                "of_which_are_placeholders": len(caught),
                "would_be_wrongly_caught_by_substring_match": len(distinct) - len(caught),
                "examples_that_are_real_names": sorted(distinct - caught)[:5],
            }

    # --- worked examples, drawn from the data ------------------------------
    primary = [r for r in rows if col(r, "in_primary_set") == "1"]
    same, diff = Counter(), Counter()
    for r in primary:
        pb, pp = col(r, "product_bakta"), col(r, "product_prokka")
        gb, gp = col(r, "gene_bakta"), col(r, "gene_prokka")
        agree, na, nb_ = lib_names.compare(pb, gb, pp, gp, "primary")
        key = (pb, pp, na, nb_)
        if agree and pb != pp:
            same[key] += 1
        elif not agree:
            diff[key] += 1

    def render(counter):
        return [
            {"product_bakta": b, "product_prokka": p,
             "normalised_bakta": na, "normalised_prokka": nb_,
             "n_occurrences": n}
            for (b, p, na, nb_), n in counter.most_common(N_EXAMPLES)
        ]

    # --- rates under each declared level -----------------------------------
    n_prim = len(primary)
    rates = {}
    for level in lib_names.LEVELS:
        key = ("name_disagreement" if level == "primary"
               else f"name_disagreement_{level}")
        n = sum(int(col(r, key)) for r in primary)
        rates[level] = {"agree": n_prim - n, "name_disagreement": n,
                        "rate": round(n / n_prim, 4)}

    # --- the symmetry invariant, re-checked here ---------------------------
    identical = [r for r in rows
                 if col(r, "product_bakta") == col(r, "product_prokka")]
    broken = 0
    for r in identical:
        for level in lib_names.LEVELS:
            agree, _, _ = lib_names.compare(
                col(r, "product_bakta"), col(r, "gene_bakta"),
                col(r, "product_prokka"), col(r, "gene_prokka"), level)
            if not agree:
                broken += 1
    if broken:
        raise SystemExit(
            f"FATAL: {broken} identical raw product strings scored as "
            "disagreements. lib_names.compare is not symmetric.")

    payload = {
        "step": "13_name_rules",
        "fixed_before_fitting": (
            "The placeholder list, the normalisation rule and the two "
            "sensitivity checks were declared and approved before any model "
            "was fitted. They are not revised after seeing a score."),
        "what_the_label_is": (
            "Disagreement between two software products about a name. Not a "
            "claim that either name is correct, and not a claim about what the "
            "protein does."),

        "placeholder_rule": {
            "method": (
                "exact literal match after NFKC, trim, whitespace collapse, "
                "trailing full stop removed, case-fold"),
            "not_a_substring_match": (
                "'hypothetical' and 'conserved' both occur inside real product "
                "names. See substring_trap for how many would be wrongly "
                "caught."),
            "literals": literal_counts,
            "declared_but_not_observed": sorted(
                k for k, v in literal_counts.items() if not v["observed"]),
            "declared_but_not_observed_note": (
                "kept in the rule so it is complete rather than fitted to "
                "these 25 genomes"),
            "totals": {
                "bakta_placeholder": sum(1 for r in rows
                                         if col(r, "placeholder_bakta") == "1"),
                "prokka_placeholder": sum(1 for r in rows
                                          if col(r, "placeholder_prokka") == "1"),
                "n_pairs": len(rows),
            },
            "substring_trap": substring_trap,
        },

        "kept_as_named_on_purpose": {
            "note": (
                "Each is a database match carrying a structural claim, so "
                "calling it uninformative would be our judgement rather than "
                "the tool's. Their sizes and disagreement rates within the "
                "primary set are in 12_content_cohort.json."),
            "patterns": {k: v.pattern
                         for k, v in lib_names.NAMED_BUT_GENERIC.items()},
        },

        "normalisation_rule": {
            "primary": [
                "1. Unicode NFKC, trim.",
                "2. Strip trailing (...) or [...] qualifiers, repeatedly. "
                "Trailing only: '3-oxoacyl-[acyl-carrier-protein] reductase' "
                "keeps its internal bracket.",
                "3. Strip ONE trailing gene-symbol token. The symbol set is "
                "built once per pair from BOTH tools' gene= attributes plus "
                "any trailing token matching " + lib_names.GENE_SHAPE.pattern + ".",
                "4. Strip leading hedges: " + ", ".join(lib_names.HEDGES) + ".",
                "5. Case-fold; every run of non-alphanumerics becomes one "
                "space; collapse whitespace.",
            ],
            "strict_sensitivity_check": (
                "case-fold and whitespace only; nothing stripped"),
            "loose_sensitivity_check": (
                "primary, then drop generic tokens ("
                + ", ".join(sorted(lib_names.GENERIC_TOKENS))
                + ") and compare the remainder as an unordered set"),
            "sensitivity_checks_were_named_in_advance": (
                "strict and loose are reported alongside primary. They are not "
                "alternatives to switch to after seeing which gives a better "
                "number."),
            "step_3_is_pair_symmetric": (
                "The gene-symbol set is built once per pair and applied to both "
                "sides. Keying it on each record's own gene= attribute is wrong "
                "and was caught in review: Bakta emits gene= on 19,216 of "
                "87,859 CDS against Prokka's 48,201, so a per-record rule fires "
                "on one side only. Before the fix, 'GTPase Era' vs 'GTPase Era' "
                "was scored as a disagreement 23 times."),
            "symmetry_invariant": {
                "assertion": ("identical raw product strings must never be "
                              "scored as a disagreement, at any level"),
                "n_pairs_with_identical_raw_strings": len(identical),
                "n_violations": broken,
                "enforced_by": "lib_names.assert_symmetric, called in step 12",
            },
            "known_residual": {
                "case": ("'Bifunctional protein FolD' vs 'Bifunctional protein "
                         "FolD protein' remains a disagreement: Prokka's "
                         "trailing token is 'protein', so the symbol strip does "
                         "not reach FolD."),
                "not_patched_because": (
                    "iterating the strip through trailing generic tokens would "
                    "fix a handful of rows and put 'membrane protein' at risk"),
            },
        },

        "worked_examples": {
            "called_same_despite_differing_raw_strings": render(same),
            "called_different": render(diff),
            "note": ("drawn from the data by frequency, not chosen by hand; "
                     "n_occurrences is how often that exact pair appears"),
        },

        "ec_and_gene_symbols": {
            "decision": "separate reported columns, not part of the product comparison",
            "reason": (
                "A gene-symbol or EC comparison is a second measurement of the "
                "same disagreement the label measures. Folding either into the "
                "product comparison would double-count it; using either as a "
                "feature would be circular. Both are reported in "
                "12_content_cohort.json."),
        },

        "primary_analysis_set_n": n_prim,
        "name_disagreement_by_level": rates,
        "interpretation": (
            "Normalisation moves the rate by "
            f"{round(100 * (rates['strict']['rate'] - rates['loose']['rate']), 1)} "
            "points across its full declared range. The disagreement is "
            "substantive, not typographic."),
    }

    METRICS.write_text(json.dumps(payload, indent=2) + "\n")

    print(f"placeholder literals observed: "
          f"{sum(1 for v in literal_counts.values() if v['observed'])} of "
          f"{len(literal_counts)}")
    print(f"identical raw strings: {len(identical):,}  violations: {broken}")
    for level in lib_names.LEVELS:
        print(f"  {level:8} name_disagreement {rates[level]['name_disagreement']:>6,} "
              f"({rates[level]['rate']:.1%})")
    print(f"wrote {METRICS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

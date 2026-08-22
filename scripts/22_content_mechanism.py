#!/usr/bin/env python3
"""Does thin database evidence explain the naming disagreement?

Reads:  data/interim/content_pairs.tsv.gz
        data/interim/content_features.tsv.gz
Writes: results/metrics/22_content_mechanism.json

Step 18 found the forest reads db-derived features almost exclusively, and
step 19 found sequence-only performance near the floor. The obvious reading is
a causal chain: weak reference coverage -> Bakta falls back on a structural
name -> Prokka names it differently -> disagreement.

This step tests that reading rather than asserting it, because it is cheap and
because the chain has an observable consequence: disagreement should fall as
database evidence deepens. It does on Prokka's side. It does the opposite on
Bakta's, and the reversal is the interesting part.

Nothing here is a model. These are conditional rates over the primary set.
"""

import gzip
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_names

ROOT = Path(__file__).resolve().parent.parent
PAIRS = ROOT / "data" / "interim" / "content_pairs.tsv.gz"
FEATURES = ROOT / "data" / "interim" / "content_features.tsv.gz"
METRICS = ROOT / "results" / "metrics" / "22_content_mechanism.json"


def read(path):
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            yield line.rstrip("\n").split("\t"), idx


def main():
    label, product = {}, {}
    for f, i in read(PAIRS):
        if f[i["in_primary_set"]] != "1":
            continue
        key = (f[i["genome"]], f[i["seqid"]], f[i["start"]], f[i["end"]])
        label[key] = int(f[i["name_disagreement"]])
        product[key] = f[i["product_bakta"]]
    if not label:
        raise SystemExit("FATAL: primary analysis set is empty")

    # --- the three declared Bakta fallback-naming families -----------------
    # Matched independently, then unioned. Their per-pattern counts overlap:
    # every DUF name in this panel also ends in 'domain-containing protein',
    # so summing the three double-counts and the union is the honest total.
    per_pattern = {
        name: {k for k, p in product.items() if pattern.search(p.strip())}
        for name, pattern in lib_names.NAMED_BUT_GENERIC.items()
    }
    union = set().union(*per_pattern.values())
    remainder = set(label) - union

    def rate(keys):
        keys = list(keys)
        return {"n": len(keys),
                "n_disagree": sum(label[k] for k in keys),
                "rate": round(sum(label[k] for k in keys) / len(keys), 4)} \
            if keys else {"n": 0, "n_disagree": 0, "rate": None}

    # --- evidence depth, on the remainder only -----------------------------
    # The families are excluded here because they are ~100% disagreement by
    # construction; leaving them in would let them drive every bucket.
    buckets = defaultdict(list)
    prokka_ec = defaultdict(list)
    prokka_sim = defaultdict(list)
    fam_ev, rem_ev = Counter(), Counter()
    fam_flags = defaultdict(int)
    rem_flags = defaultdict(int)

    for f, i in read(FEATURES):
        if f[i["in_primary_set"]] != "1":
            continue
        key = (f[i["genome"]], f[i["seqid"]], f[i["start"]], f[i["end"]])
        if key not in label:
            continue
        n_dbxref = int(f[i["bakta_n_dbxref"]])
        in_family = key in union
        (fam_ev if in_family else rem_ev)[n_dbxref] += 1
        target = fam_flags if in_family else rem_flags
        target["n"] += 1
        for flag in ("prokka_has_ec", "prokka_has_similarity_hit",
                     "prokka_has_protein_motif"):
            target[flag] += int(f[i[flag]] == "1")
        if in_family:
            continue
        bucket = ("2_the_minimum" if n_dbxref <= 2
                  else "3" if n_dbxref == 3 else "4_or_more")
        buckets[bucket].append(key)
        prokka_ec["with_ec" if f[i["prokka_has_ec"]] == "1"
                  else "without_ec"].append(key)
        prokka_sim["with_similarity_hit"
                   if f[i["prokka_has_similarity_hit"]] == "1"
                   else "without_similarity_hit"].append(key)

    def shares(d):
        n = d["n"]
        return {k: round(v / n, 4) for k, v in d.items() if k != "n"} | {"n": n}

    payload = {
        "step": "22_content_mechanism",
        "question": (
            "Is the naming disagreement explained by thin reference-database "
            "evidence? Tested rather than asserted."),
        "not_a_model": "conditional rates over the primary analysis set",

        "fallback_naming_families": {
            "note": (
                "Three Bakta naming conventions declared before fitting in "
                "12_content_cohort.json. Matched independently they total "
                f"{sum(len(v) for v in per_pattern.values()):,} rows, but the "
                "patterns overlap -- every DUF name in this panel also ends in "
                "'domain-containing protein' -- so the honest total is the "
                f"distinct union, {len(union):,}."),
            "per_pattern_counts_overlapping": {
                k: len(v) for k, v in sorted(per_pattern.items())},
            "sum_of_patterns_double_counts": sum(
                len(v) for v in per_pattern.values()),
            "distinct_union": len(union),
            "duf_entirely_inside_domain_containing": len(
                per_pattern["duf"] & per_pattern["domain_containing"]
            ) == len(per_pattern["duf"]),
            "families": rate(union),
            "remainder": rate(remainder),
            "whole_primary_set": rate(label),
            "headline_pair": (
                f"{rate(label)['rate']:.1%} of both-named regions are named "
                f"differently by the two tools; "
                f"{rate(remainder)['rate']:.1%} once three Bakta "
                "fallback-naming conventions that disagree by construction are "
                "excluded. Both numbers should be published together."),
        },

        "bakta_evidence_depth_remainder_only": {
            "excludes": ("the fallback-naming families, which are ~100% "
                         "disagreement by construction and would drive every "
                         "bucket"),
            "buckets": {k: rate(v) for k, v in sorted(buckets.items())},
            "finding": (
                "Disagreement RISES with Bakta cross-reference count, from "
                f"{rate(buckets['2_the_minimum'])['rate']:.1%} at the minimum "
                f"two to {rate(buckets['3'])['rate']:.1%} at three. This is "
                "the opposite of the thin-evidence chain. A third "
                "cross-reference is typically an EC number, a BlastRules hit, "
                "a virulence-factor or insertion-sequence match, and it comes "
                "with a MORE specific Bakta name -- which is then more likely "
                "to differ from Prokka's more generic one."),
        },

        "prokka_evidence_depth_remainder_only": {
            "by_ec": {k: rate(v) for k, v in sorted(prokka_ec.items())},
            "by_similarity_hit": {k: rate(v) for k, v in sorted(prokka_sim.items())},
            "finding": (
                "On Prokka's side the thin-evidence direction does hold: "
                f"{rate(prokka_ec['without_ec'])['rate']:.1%} disagreement "
                "without an EC number against "
                f"{rate(prokka_ec['with_ec'])['rate']:.1%} with one."),
        },

        "evidence_profile_families_vs_remainder": {
            "bakta_n_dbxref_distribution": {
                "families": dict(sorted(fam_ev.items())),
                "remainder": dict(sorted(rem_ev.items())),
                "note": (
                    "Family rows carry exactly the minimum two "
                    "cross-references almost without exception; the remainder "
                    "carries three or more about one time in ten."),
            },
            "prokka_flag_shares": {
                "families": shares(dict(fam_flags)),
                "remainder": shares(dict(rem_flags)),
            },
        },

        "conclusion": (
            "The chain holds for the fallback-naming families and on Prokka's "
            "side, and reverses on Bakta's. There are two mechanisms pulling "
            "in opposite directions, not one. bakta_n_dbxref is the forest's "
            "strongest feature because extra Bakta evidence predicts a more "
            "specific Bakta name that Prokka does not match -- not because "
            "thin evidence predicts a fallback name. The db-derived features "
            "and the naming families are related, but they are not two views "
            "of a single mechanism."),
        "what_this_does_not_show": (
            "Nothing here identifies which tool is right. The label remains "
            "disagreement between two software products about a name."),
    }
    METRICS.write_text(json.dumps(payload, indent=2) + "\n")

    fam = payload["fallback_naming_families"]
    print(f"families {fam['families']['n']:,} rows "
          f"({fam['families']['rate']:.1%})   "
          f"remainder {fam['remainder']['n']:,} rows "
          f"({fam['remainder']['rate']:.1%})   "
          f"overall {fam['whole_primary_set']['rate']:.1%}")
    print("\nremainder, by Bakta cross-reference count:")
    for k, v in sorted(payload["bakta_evidence_depth_remainder_only"]["buckets"].items()):
        print(f"  {k:14s} {v['n']:6,} rows  {v['rate']:.1%}")
    print("remainder, by Prokka EC:")
    for k, v in sorted(payload["prokka_evidence_depth_remainder_only"]["by_ec"].items()):
        print(f"  {k:14s} {v['n']:6,} rows  {v['rate']:.1%}")
    print(f"\nwrote {METRICS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

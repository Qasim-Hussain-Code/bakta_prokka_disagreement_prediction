#!/usr/bin/env python3
"""Find intervals called by one tool and not the other. This is the label.

Reads:  data/interim/calls.tsv.gz
Writes: data/interim/disagreement.tsv.gz   (gitignored)
        results/metrics/05_disagreement.json

The unit of analysis is one Bakta CDS call. The label is:

    1  no Prokka feature overlaps this interval on the same strand
    0  some Prokka feature does

Three choices are doing real work here, and each could reasonably have gone
the other way:

  Same strand only. Two genes on opposite strands can occupy the same
  coordinates and still be different genes. Counting an opposite-strand
  overlap as agreement would quietly mark real disagreements as matched.

  Any Prokka feature type, not only CDS. If Prokka called a tRNA where Bakta
  called a CDS, Prokka did not leave the region empty -- it read the region
  differently. That is a different phenomenon from silence, and folding it
  into the positive class would put two unlike things in one label. Those
  cases are labelled 0 and kept identifiable via overlap_category.

  One base pair of overlap is enough to count as a match. Gene callers
  routinely pick different start codons for the same gene; requiring
  reciprocal overlap would relabel start-codon disagreements as missing
  genes. Boundary disagreement is recorded separately, as same_start and
  same_stop, rather than being mixed into the label.

The label is tool output. A positive case means Bakta made a call where
Prokka did not. It does not mean a gene is present.
"""

import gzip
import json
from bisect import bisect_left
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CALLS = ROOT / "data" / "interim" / "calls.tsv.gz"
ACCESSIONS = ROOT / "data" / "accessions.tsv"
OUT_TSV = ROOT / "data" / "interim" / "disagreement.tsv.gz"
METRICS = ROOT / "results" / "metrics" / "05_disagreement.json"

COLUMNS = [
    "genome", "species", "seqid", "start", "end", "strand", "length_bp",
    "locus_tag", "gene", "product",
    "label", "overlap_category", "n_overlaps",
    "best_overlap_bp", "best_overlap_frac", "best_overlap_ftype",
    "same_start", "same_stop",
]


def species_map():
    out = {}
    for line in ACCESSIONS.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        p = line.split("\t")
        if p[0] == "species":
            continue
        out[p[3]] = p[0]
    return out


def load_calls():
    """calls.tsv.gz -> {(genome, tool)} -> {(seqid, strand): [intervals]}"""
    by_tool = defaultdict(lambda: defaultdict(list))
    with gzip.open(CALLS, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            by_tool[(f[idx["genome"]], f[idx["tool"]])][
                (f[idx["seqid"]], f[idx["strand"]])
            ].append({
                "start": int(f[idx["start"]]),
                "end": int(f[idx["end"]]),
                "ftype": f[idx["ftype"]],
                "locus_tag": f[idx["locus_tag"]],
                "gene": f[idx["gene"]],
                "product": f[idx["product"]],
            })
    return by_tool


def index_intervals(intervals):
    """Sort once and keep a parallel list of starts for bisect."""
    ordered = sorted(intervals, key=lambda r: (r["start"], r["end"]))
    return ordered, [r["start"] for r in ordered]


def overlaps(query, ordered, starts, max_len):
    """Every interval overlapping query, 1-based inclusive on both ends.

    Scanning back by the longest feature on the contig bounds the window:
    nothing starting before query.start - max_len can still reach query.
    """
    if not ordered:
        return []
    lo = bisect_left(starts, query["start"] - max_len)
    hits = []
    for r in ordered[lo:]:
        if r["start"] > query["end"]:
            break
        if r["end"] >= query["start"]:
            hits.append(r)
    return hits


def main():
    if not CALLS.exists():
        raise SystemExit(f"{CALLS} absent; run scripts/04_extract_calls.py first")

    species = species_map()
    by_tool = load_calls()
    genomes = sorted({g for (g, _t) in by_tool})

    rows = []
    per_genome = []

    for genome in genomes:
        bakta = by_tool[(genome, "bakta")]
        prokka = by_tool[(genome, "prokka")]

        prokka_idx = {}
        for key, intervals in prokka.items():
            ordered, starts = index_intervals(intervals)
            longest = max((r["end"] - r["start"] + 1) for r in ordered)
            prokka_idx[key] = (ordered, starts, longest)

        counts = Counter()
        for (seqid, strand), intervals in bakta.items():
            ordered, starts, longest = prokka_idx.get((seqid, strand), ([], [], 1))
            for cds in intervals:
                if cds["ftype"] != "CDS":
                    continue
                hits = overlaps(cds, ordered, starts, longest)
                length = cds["end"] - cds["start"] + 1

                best_bp, best_ftype, best = 0, "", None
                for h in hits:
                    bp = min(cds["end"], h["end"]) - max(cds["start"], h["start"]) + 1
                    if bp > best_bp:
                        best_bp, best_ftype, best = bp, h["ftype"], h

                if not hits:
                    label, category = 1, "none"
                elif any(h["ftype"] == "CDS" for h in hits):
                    label, category = 0, "cds"
                else:
                    label, category = 0, "non_cds"
                counts[category] += 1

                rows.append({
                    "genome": genome,
                    "species": species.get(genome, ""),
                    "seqid": seqid,
                    "start": cds["start"],
                    "end": cds["end"],
                    "strand": strand,
                    "length_bp": length,
                    "locus_tag": cds["locus_tag"],
                    "gene": cds["gene"],
                    "product": cds["product"],
                    "label": label,
                    "overlap_category": category,
                    "n_overlaps": len(hits),
                    "best_overlap_bp": best_bp,
                    "best_overlap_frac": round(best_bp / length, 4) if length else 0,
                    "best_overlap_ftype": best_ftype,
                    "same_start": int(bool(best) and best["start"] == cds["start"]),
                    "same_stop": int(bool(best) and best["end"] == cds["end"]),
                })

        # The reverse direction is not the label, but reporting only one
        # direction would let an asymmetry in the tools pass unnoticed.
        bakta_idx = {}
        for key, intervals in bakta.items():
            ordered, starts = index_intervals(intervals)
            longest = max((r["end"] - r["start"] + 1) for r in ordered)
            bakta_idx[key] = (ordered, starts, longest)

        prokka_unique = 0
        prokka_cds = 0
        for (seqid, strand), intervals in prokka.items():
            ordered, starts, longest = bakta_idx.get((seqid, strand), ([], [], 1))
            for cds in intervals:
                if cds["ftype"] != "CDS":
                    continue
                prokka_cds += 1
                if not overlaps(cds, ordered, starts, longest):
                    prokka_unique += 1

        n_bakta_cds = sum(counts.values())
        per_genome.append({
            "genome": genome,
            "species": species.get(genome, ""),
            "bakta_cds": n_bakta_cds,
            "prokka_cds": prokka_cds,
            "unique_to_bakta": counts["none"],
            "unique_to_bakta_rate": round(counts["none"] / n_bakta_cds, 4) if n_bakta_cds else 0,
            "matched_cds": counts["cds"],
            "matched_non_cds": counts["non_cds"],
            "unique_to_prokka": prokka_unique,
        })

    with gzip.open(OUT_TSV, "wt", newline="") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")

    pos = sum(r["label"] for r in rows)
    matched = [r for r in rows if r["label"] == 0]
    payload = {
        "step": "05_build_disagreement_set",
        "unit": "one Bakta CDS call",
        "label_definition": (
            "1 = no Prokka feature of any type overlaps this interval by at "
            "least 1 bp on the same strand; 0 = one does"
        ),
        "label_is_tool_output": (
            "A positive case means Bakta made a call where Prokka did not. It "
            "does not mean a gene is present."
        ),
        "overlap_rule": "same strand, >=1 bp, any Prokka feature type",
        "n_rows": len(rows),
        "n_positive": pos,
        "positive_rate": round(pos / len(rows), 4) if rows else 0,
        "category_counts": dict(Counter(r["overlap_category"] for r in rows)),
        "boundary_agreement_among_matched": {
            "n_matched": len(matched),
            "same_start_and_stop": sum(
                1 for r in matched if r["same_start"] and r["same_stop"]),
            "same_stop_only": sum(
                1 for r in matched if r["same_stop"] and not r["same_start"]),
            "same_start_only": sum(
                1 for r in matched if r["same_start"] and not r["same_stop"]),
            "neither": sum(
                1 for r in matched if not r["same_start"] and not r["same_stop"]),
            "note": (
                "Start-codon choice is where these tools disagree most often "
                "while still calling the same gene. These are label 0."
            ),
        },
        "table": str(OUT_TSV.relative_to(ROOT)),
        "per_genome": per_genome,
    }
    METRICS.write_text(json.dumps(payload, indent=2) + "\n")

    print(f"{len(rows):,} Bakta CDS calls across {len(genomes)} genomes")
    print(f"  positive (Prokka silent): {pos:,} "
          f"({payload['positive_rate'] * 100:.2f}%)")
    print(f"  matched by Prokka CDS:    {payload['category_counts'].get('cds', 0):,}")
    print(f"  matched by non-CDS:       {payload['category_counts'].get('non_cds', 0):,}")
    print(f"wrote {OUT_TSV.relative_to(ROOT)}")
    print(f"wrote {METRICS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

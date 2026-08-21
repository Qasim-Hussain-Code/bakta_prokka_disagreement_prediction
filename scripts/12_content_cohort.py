#!/usr/bin/env python3
"""Pair each Bakta CDS with the Prokka CDS at the same locus, and decide which
pairs the naming question can be asked of.

Reads:  data/annotations/bakta/<acc>/<acc>.gff3
        data/annotations/prokka/<acc>/<acc>.gff
Writes: data/interim/content_pairs.tsv.gz   (gitignored)
        results/metrics/12_content_cohort.json

No new annotation and no new compute beyond parsing. Part one established that
these two tools call almost exactly the same intervals (87,788 of 87,888
matched calls are coordinate-identical). This step takes that as given and
asks the next question: given a region both tools called, do they name it the
same thing?

Why this re-parses the GFFs instead of reusing data/interim/calls.tsv.gz:
04_extract_calls.py reads the gene symbol as `a.get("gene", a.get("Name", ""))`.
Bakta sets Name=product on every CDS, so that fallback filled the gene column
with product strings and `has_gene_symbol` came out constant 1 for all 87,960
part-one rows. A constant feature can never be split on, so no part-one number
moves -- but the column is unusable here, where the true gene symbol is needed
both for normalisation and as a separately reported comparison. Bakta's real
gene= is present on 19,216 CDS, not all of them. Scripts 01-11 are left exactly
as they are; this step simply does not depend on that column.

The pairing rule, fixed in advance:

  One row per Bakta CDS. Its partner is the Prokka CDS on the same sequence
  and the same strand with the largest overlap in bp. A Bakta CDS whose only
  overlap is with a non-CDS Prokka feature, or with nothing at all, has no
  second name to compare and is out of scope -- recorded, not silently
  dropped.

The primary analysis set is the pairs where BOTH tools assigned a product that
is not a placeholder. That restriction is imposed by what can be measured: you
cannot compare two names when one tool did not produce one. The four cells of
the placeholder cross-tabulation are all reported.
"""

import gzip
import json
import re
import sys
from bisect import bisect_left
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_names

ROOT = Path(__file__).resolve().parent.parent
ACCESSIONS = ROOT / "data" / "accessions.tsv"
ANN = ROOT / "data" / "annotations"
OUT_TSV = ROOT / "data" / "interim" / "content_pairs.tsv.gz"
METRICS = ROOT / "results" / "metrics" / "12_content_cohort.json"

TOOLS = {"bakta": "gff3", "prokka": "gff"}

COLUMNS = [
    "genome", "species", "seqid", "strand", "start", "end", "length_bp",
    "overlap_bp", "same_start", "same_stop",
    "locus_tag_bakta", "locus_tag_prokka",
    "product_bakta", "product_prokka",
    "gene_bakta", "gene_prokka",
    "ec_bakta", "ec_prokka",
    "dbxref_bakta", "dbxref_prokka", "inference_prokka",
    "placeholder_bakta", "placeholder_prokka", "in_primary_set",
    "generic_family",
    "name_disagreement", "name_disagreement_strict", "name_disagreement_loose",
    "symbol_comparable", "symbol_disagreement",
    "ec_comparable", "ec_disagreement",
]


def read_panel():
    rows = []
    for line in ACCESSIONS.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        p = line.split("\t")
        if p[0] == "species":
            continue
        rows.append({"species": p[0], "accession": p[3]})
    if not rows:
        raise SystemExit("no accessions; run scripts/01_fetch_genomes.py first")
    return rows


def parse_attributes(field):
    out = {}
    for chunk in field.split(";"):
        chunk = chunk.strip()
        if not chunk or "=" not in chunk:
            continue
        k, _, v = chunk.partition("=")
        out[unquote(k)] = unquote(v)
    return out


def read_cds(path, tool):
    """CDS features only, with the attributes the naming question needs."""
    out = defaultdict(list)
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != "CDS":
                continue
            a = parse_attributes(f[8])
            n += 1
            out[(f[0], f[6])].append({
                "start": int(f[3]), "end": int(f[4]),
                "locus_tag": a.get("locus_tag", ""),
                "product": a.get("product", ""),
                # raw gene=, with no Name fallback. See the module docstring.
                "gene": a.get("gene", ""),
                "ec": a.get("eC_number", ""),
                "dbxref": a.get("Dbxref", a.get("db_xref", "")),
                "inference": a.get("inference", ""),
            })
    if n == 0:
        raise SystemExit(f"FATAL: no CDS features parsed from {path}")
    return out, n


def best_partner(cds, ordered, starts, longest):
    """Largest-overlap Prokka CDS, or None."""
    if not ordered:
        return None, 0
    lo = bisect_left(starts, cds["start"] - longest)
    best, best_bp = None, 0
    for r in ordered[lo:]:
        if r["start"] > cds["end"]:
            break
        if r["end"] < cds["start"]:
            continue
        bp = min(cds["end"], r["end"]) - max(cds["start"], r["start"]) + 1
        if bp > best_bp:
            best, best_bp = r, bp
    return best, best_bp


def bakta_ec(dbxref):
    return sorted(x[3:] for x in dbxref.split(",") if x.startswith("EC:"))


def generic_family(product):
    for name, pattern in lib_names.NAMED_BUT_GENERIC.items():
        if pattern.search(product.strip()):
            return name
    return ""


def clean(value):
    """GFF fields are tab-delimited already, so a tab cannot appear inside one.
    Guard anyway: a stray tab would silently shift every later column."""
    return str(value).replace("\t", " ").replace("\n", " ")


def main():
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    METRICS.parent.mkdir(parents=True, exist_ok=True)

    panel = read_panel()
    pairs, per_genome = [], []
    unpaired_total = 0
    tool_cds_totals = Counter()

    for entry in panel:
        acc, species = entry["accession"], entry["species"]
        cds = {}
        for tool, ext in TOOLS.items():
            path = ANN / tool / acc / f"{acc}.{ext}"
            if not path.exists():
                raise SystemExit(
                    f"missing {tool} annotation for {acc}: {path}\n"
                    f"part two runs on the annotation already on disk; it does "
                    f"not re-annotate."
                )
            cds[tool], n = read_cds(path, tool)
            tool_cds_totals[tool] += n

        index = {}
        for key, intervals in cds["prokka"].items():
            ordered = sorted(intervals, key=lambda r: (r["start"], r["end"]))
            index[key] = (ordered, [r["start"] for r in ordered],
                          max(r["end"] - r["start"] + 1 for r in ordered))

        n_bakta = n_paired = n_unpaired = 0
        for (seqid, strand), intervals in cds["bakta"].items():
            ordered, starts, longest = index.get((seqid, strand), ([], [], 1))
            for c in intervals:
                n_bakta += 1
                partner, overlap = best_partner(c, ordered, starts, longest)
                if partner is None:
                    n_unpaired += 1
                    continue
                n_paired += 1

                ph_b = lib_names.is_placeholder(c["product"])
                ph_p = lib_names.is_placeholder(partner["product"])
                primary = not ph_b and not ph_p

                row = {
                    "genome": acc, "species": species, "seqid": seqid,
                    "strand": strand, "start": c["start"], "end": c["end"],
                    "length_bp": c["end"] - c["start"] + 1,
                    "overlap_bp": overlap,
                    "same_start": int(partner["start"] == c["start"]),
                    "same_stop": int(partner["end"] == c["end"]),
                    "locus_tag_bakta": c["locus_tag"],
                    "locus_tag_prokka": partner["locus_tag"],
                    "product_bakta": c["product"],
                    "product_prokka": partner["product"],
                    "gene_bakta": c["gene"], "gene_prokka": partner["gene"],
                    "ec_bakta": ";".join(bakta_ec(c["dbxref"])),
                    "ec_prokka": partner["ec"],
                    "dbxref_bakta": c["dbxref"], "dbxref_prokka": partner["dbxref"],
                    "inference_prokka": partner["inference"],
                    "placeholder_bakta": int(ph_b),
                    "placeholder_prokka": int(ph_p),
                    "in_primary_set": int(primary),
                    "generic_family": generic_family(c["product"]) if primary else "",
                }

                # Labels at all three declared levels. Only the primary level
                # is modelled; strict and loose were named as sensitivity
                # checks before any of these numbers existed.
                for level in lib_names.LEVELS:
                    agree, _, _ = lib_names.compare(
                        c["product"], c["gene"], partner["product"],
                        partner["gene"], level)
                    key = ("name_disagreement" if level == "primary"
                           else f"name_disagreement_{level}")
                    row[key] = int(not agree) if primary else ""

                # Gene symbol and EC number are separate reported columns, not
                # part of the product comparison and never used as features. A
                # symbol comparison is a second measurement of the same thing
                # the label measures; feeding it in would be circular.
                sym_ok = bool(c["gene"].strip()) and bool(partner["gene"].strip())
                row["symbol_comparable"] = int(sym_ok and primary)
                row["symbol_disagreement"] = (
                    int(c["gene"].strip().lower() != partner["gene"].strip().lower())
                    if (sym_ok and primary) else "")
                ec_b, ec_p = bakta_ec(c["dbxref"]), partner["ec"]
                ec_ok = bool(ec_b) and bool(ec_p)
                row["ec_comparable"] = int(ec_ok and primary)
                row["ec_disagreement"] = (
                    int(set(ec_b) != set(x for x in ec_p.split(",") if x))
                    if (ec_ok and primary) else "")

                pairs.append(row)

        unpaired_total += n_unpaired
        prim = [p for p in pairs if p["genome"] == acc and p["in_primary_set"]]
        per_genome.append({
            "genome": acc, "species": species,
            "bakta_cds": n_bakta, "paired": n_paired, "unpaired": n_unpaired,
            "primary_set": len(prim),
            "name_disagreement": sum(p["name_disagreement"] for p in prim),
            "name_disagreement_rate": (
                round(sum(p["name_disagreement"] for p in prim) / len(prim), 4)
                if prim else None),
        })
        print(f"  {acc} {n_bakta:>6,} CDS  {n_paired:>6,} paired  "
              f"{len(prim):>6,} primary", flush=True)

    if not pairs:
        raise SystemExit("FATAL: no pairs built")

    # The rule must never split a string from itself. See lib_names.
    lib_names.assert_symmetric(pairs)

    with gzip.open(OUT_TSV, "wt", newline="") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in pairs:
            fh.write("\t".join(clean(r[c]) for c in COLUMNS) + "\n")

    # --- the cross-tabulation ---------------------------------------------
    cross = Counter((r["placeholder_bakta"], r["placeholder_prokka"]) for r in pairs)
    n_pairs = len(pairs)
    primary = [r for r in pairs if r["in_primary_set"]]
    n_prim = len(primary)
    if n_prim == 0:
        raise SystemExit("FATAL: primary analysis set is empty")

    fam = Counter(r["generic_family"] for r in primary)
    fam_rates = {}
    for name in list(lib_names.NAMED_BUT_GENERIC) + [""]:
        sub = [r for r in primary if r["generic_family"] == name]
        if sub:
            fam_rates[name or "none_of_these"] = {
                "n": len(sub),
                "share_of_primary_set": round(len(sub) / n_prim, 4),
                "name_disagreement_rate": round(
                    sum(r["name_disagreement"] for r in sub) / len(sub), 4),
            }

    sym = [r for r in primary if r["symbol_comparable"]]
    ec = [r for r in primary if r["ec_comparable"]]

    payload = {
        "step": "12_content_cohort",
        "unit": "one Bakta CDS call paired with the Prokka CDS at the same locus",
        "label_definition": (
            "name_disagreement = 1 when the two tools' product strings differ "
            "after the primary normalisation rule in lib_names.py. This is "
            "disagreement between two software products about a name. It is "
            "not a claim that either name is correct, and it is not a claim "
            "about what the protein does."),
        "pairing_rule": (
            "same sequence, same strand, largest overlap in bp, Prokka feature "
            "must be a CDS"),
        "no_new_annotation": (
            "runs on the annotation output already on disk; the Bakta database "
            "is not required and has been removed from this machine"),
        "bakta_cds_total": tool_cds_totals["bakta"],
        "prokka_cds_total": tool_cds_totals["prokka"],
        "n_paired": n_pairs,
        "n_unpaired": unpaired_total,
        "unpaired_note": (
            "Bakta CDS with no overlapping same-strand Prokka CDS. There is no "
            "second name to compare, so the naming question cannot be asked of "
            "them. They are part one's label, not part two's."),
        "coordinate_agreement_among_pairs": {
            "identical_start_and_stop": sum(
                1 for r in pairs if r["same_start"] and r["same_stop"]),
            "note": ("carried over from part one for context: these tools agree "
                     "on where genes are almost always, which is why the "
                     "naming question is the one with signal in it"),
        },
        "cross_tabulation": {
            "both_named": cross[(0, 0)],
            "both_placeholder": cross[(1, 1)],
            "bakta_named_only": cross[(0, 1)],
            "prokka_named_only": cross[(1, 0)],
            "total": n_pairs,
            "shares": {
                "both_named": round(cross[(0, 0)] / n_pairs, 4),
                "both_placeholder": round(cross[(1, 1)] / n_pairs, 4),
                "bakta_named_only": round(cross[(0, 1)] / n_pairs, 4),
                "prokka_named_only": round(cross[(1, 0)] / n_pairs, 4),
            },
        },
        "asymmetric_cells_interpretation": {
            "finding": (
                "The large asymmetric cell is 'Bakta named, Prokka not' "
                f"({cross[(0, 1)]:,}), not the reverse ({cross[(1, 0)]:,}). "
                "Bakta on db-light leaves "
                f"{round(100 * (cross[(1, 1)] + cross[(1, 0)]) / n_pairs, 1)}% "
                "of paired CDS as a placeholder against Prokka's "
                f"{round(100 * (cross[(1, 1)] + cross[(0, 1)]) / n_pairs, 1)}%."),
            "db_light_direction": (
                "db-light was used because db-full does not fit on this "
                "machine. It has fewer reference proteins, so it can only "
                "make Bakta name FEWER regions than db-full would. It "
                "therefore works AGAINST the observed asymmetry rather than "
                "producing it: a db-full run would move rows out of "
                "'both_placeholder' and 'prokka_named_only' into "
                "'bakta_named_only' and make the imbalance larger."),
            "therefore": (
                "bakta_named_only is a LOWER BOUND under db-light. The "
                "asymmetry is a difference in reference-database breadth and "
                "naming policy between the two tools, of which the db-light "
                "constraint is one component pushing in a known direction. It "
                "is not an artefact that a db-full run would remove."),
            "what_is_still_confounded": (
                "The absolute size of 'both_placeholder' and "
                "'prokka_named_only' is affected by db-light and should not be "
                "read as a property of Bakta in general."),
        },
        "primary_analysis_set": {
            "n": n_prim,
            "share_of_pairs": round(n_prim / n_pairs, 4),
            "definition": "both tools assigned a non-placeholder product",
            "why_conditional": (
                "Comparing two names requires both tools to have produced one. "
                "This restriction is imposed by what is measurable under "
                "db-light, not chosen for scientific convenience."),
            "name_disagreement": {
                level: {
                    "n": sum(r[key] for r in primary),
                    "rate": round(sum(r[key] for r in primary) / n_prim, 4),
                }
                for level, key in (
                    ("primary", "name_disagreement"),
                    ("strict", "name_disagreement_strict"),
                    ("loose", "name_disagreement_loose"))
            },
            "modelled_level": "primary",
        },
        "declared_generic_families_within_primary_set": {
            "note": (
                "Declared before fitting. These Bakta naming styles are kept as "
                "NAMED because each is a database match carrying a structural "
                "claim, but they are near-deterministic disagreements and their "
                "size is recorded here so the positive class is not read as "
                "uniform."),
            "families": fam_rates,
        },
        "separate_columns_not_part_of_the_label": {
            "note": (
                "Gene symbols and EC numbers are reported, never used as "
                "features. Either would be a second measurement of the label."),
            "gene_symbol": {
                "n_comparable": len(sym),
                "n_disagree": sum(r["symbol_disagreement"] for r in sym),
                "rate": (round(sum(r["symbol_disagreement"] for r in sym) / len(sym), 4)
                         if sym else None),
            },
            "ec_number": {
                "n_comparable": len(ec),
                "n_disagree": sum(r["ec_disagreement"] for r in ec),
                "rate": (round(sum(r["ec_disagreement"] for r in ec) / len(ec), 4)
                         if ec else None),
                "source": "Bakta: Dbxref EC:; Prokka: eC_number attribute",
            },
        },
        "per_genome": sorted(per_genome, key=lambda r: r["genome"]),
        "table": str(OUT_TSV.relative_to(ROOT)),
    }
    METRICS.write_text(json.dumps(payload, indent=2) + "\n")

    ct = payload["cross_tabulation"]
    print(f"\n{n_pairs:,} paired regions, {unpaired_total:,} unpaired")
    print(f"  both named        {ct['both_named']:>7,}")
    print(f"  both placeholder  {ct['both_placeholder']:>7,}")
    print(f"  Bakta named only  {ct['bakta_named_only']:>7,}")
    print(f"  Prokka named only {ct['prokka_named_only']:>7,}")
    nd = payload["primary_analysis_set"]["name_disagreement"]
    print(f"\nprimary analysis set {n_prim:,}")
    for level in ("strict", "primary", "loose"):
        print(f"  {level:8} name_disagreement {nd[level]['n']:>6,} "
              f"({nd[level]['rate']:.1%})")
    print(f"\nwrote {OUT_TSV.relative_to(ROOT)}")
    print(f"wrote {METRICS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

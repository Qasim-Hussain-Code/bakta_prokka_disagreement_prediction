#!/usr/bin/env python3
"""Parse both GFF3 sets into a common interval table.

Reads:  data/annotations/bakta/<acc>/<acc>.gff3
        data/annotations/prokka/<acc>/<acc>.gff
Writes: data/interim/calls.tsv.gz   (gitignored)
        results/metrics/04_calls.json

Every feature from both tools is kept, not only CDS. A region where Bakta
calls a CDS and Prokka calls a tRNA is a different thing from a region Prokka
leaves empty, and collapsing the two here would make that distinction
unrecoverable downstream. Step 05 decides what counts as a disagreement; this
step only reads what each tool actually emitted.

Coordinates stay exactly as GFF3 gives them: 1-based, both ends inclusive.
No conversion happens here, so nothing can be half-converted later.
"""

import gzip
import json
import re
from collections import Counter
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parent.parent
ACCESSIONS = ROOT / "data" / "accessions.tsv"
ANN = ROOT / "data" / "annotations"
OUT_TSV = ROOT / "data" / "interim" / "calls.tsv.gz"
METRICS = ROOT / "results" / "metrics" / "04_calls.json"

TOOLS = {"bakta": "gff3", "prokka": "gff"}

COLUMNS = [
    "genome", "tool", "seqid", "source", "ftype",
    "start", "end", "strand", "length_bp",
    "feature_id", "locus_tag", "gene", "product",
]


def read_accessions():
    rows = []
    for line in ACCESSIONS.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if parts[0] == "species":
            continue
        rows.append({"species": parts[0], "accession": parts[3]})
    if not rows:
        raise SystemExit("no accessions; run 01 first")
    return rows


def parse_attributes(field):
    """GFF3 attributes are URL-encoded; product names routinely contain
    commas, equals signs and semicolons that are escaped there."""
    out = {}
    for chunk in field.split(";"):
        chunk = chunk.strip()
        if not chunk or "=" not in chunk:
            continue
        key, _, value = chunk.partition("=")
        out[unquote(key)] = unquote(value)
    return out


def parse_gff(path, genome, tool):
    """Yield one record per feature line.

    Both tools append the input sequence after a ##FASTA pragma. Parsing must
    stop there: FASTA lines have no tabs, so they would not survive the column
    check below, but stopping explicitly is cheaper and states the intent.
    """
    rows = []
    seqids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            seqid, source, ftype, start, end, _score, strand, _phase, attrs = f
            a = parse_attributes(attrs)
            start, end = int(start), int(end)
            seqids.append(seqid)
            rows.append({
                "genome": genome,
                "tool": tool,
                "seqid": seqid,
                "source": source,
                "ftype": ftype,
                "start": start,
                "end": end,
                "strand": strand,
                "length_bp": end - start + 1,
                "feature_id": a.get("ID", ""),
                "locus_tag": a.get("locus_tag", ""),
                "gene": a.get("gene", a.get("Name", "")),
                "product": a.get("product", ""),
            })
    return rows, set(seqids)


def fasta_seqids(path):
    ids = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def main():
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    METRICS.parent.mkdir(parents=True, exist_ok=True)

    panel = read_accessions()
    all_rows = []
    per_genome = []
    seqid_problems = []

    for entry in panel:
        acc = entry["accession"]
        ref_ids = fasta_seqids(ROOT / "data" / "genomes" / f"{acc}.fna")
        record = {"accession": acc, "species": entry["species"],
                  "contigs": len(ref_ids), "tools": {}}
        tool_ids = {}

        for tool, ext in TOOLS.items():
            path = ANN / tool / acc / f"{acc}.{ext}"
            if not path.exists():
                raise SystemExit(
                    f"missing {tool} annotation for {acc}: {path}\n"
                    f"run scripts/0{'2' if tool == 'bakta' else '3'}_annotate_{tool}.sh"
                )
            rows, ids = parse_gff(path, acc, tool)
            all_rows.extend(rows)
            tool_ids[tool] = ids
            counts = Counter(r["ftype"] for r in rows)
            record["tools"][tool] = {
                "features": len(rows),
                "cds": counts.get("CDS", 0),
                "feature_counts": dict(sorted(counts.items())),
            }

        # The contig-naming trap. Bakta renames sequences to contig_N unless
        # --keep-contig-headers is passed. If that ever regresses, the two
        # tables share no seqid, every region looks like a disagreement, and
        # nothing downstream would notice. Fail here instead.
        for tool, ids in tool_ids.items():
            unknown = ids - ref_ids
            if unknown:
                seqid_problems.append({
                    "accession": acc, "tool": tool,
                    "unknown_seqids": sorted(unknown)[:5],
                    "n_unknown": len(unknown),
                })
        shared = tool_ids["bakta"] & tool_ids["prokka"]
        record["shared_seqids"] = len(shared)
        if not shared:
            seqid_problems.append({
                "accession": acc, "tool": "both",
                "error": "bakta and prokka share no sequence identifier",
                "bakta": sorted(tool_ids["bakta"])[:3],
                "prokka": sorted(tool_ids["prokka"])[:3],
            })
        per_genome.append(record)

    if seqid_problems:
        print(json.dumps(seqid_problems, indent=2))
        raise SystemExit(
            "FATAL: sequence identifiers do not line up between the tools and "
            "the input FASTA. Overlap comparison would be meaningless. Check "
            "that scripts/02_annotate_bakta.sh passes --keep-contig-headers."
        )

    with gzip.open(OUT_TSV, "wt", newline="") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in all_rows:
            fh.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")

    totals = {}
    for tool in TOOLS:
        rows = [r for r in all_rows if r["tool"] == tool]
        counts = Counter(r["ftype"] for r in rows)
        totals[tool] = {
            "features": len(rows),
            "cds": counts.get("CDS", 0),
            "feature_counts": dict(sorted(counts.items())),
        }

    METRICS.write_text(json.dumps({
        "step": "04_extract_calls",
        "coordinate_convention": "GFF3 as given: 1-based, both ends inclusive",
        "kept": "all feature types from both tools, not only CDS",
        "kept_note": (
            "A region one tool calls CDS and the other calls tRNA is not the "
            "same as a region the other leaves empty. Step 05 draws that line; "
            "this step does not."
        ),
        "n_genomes": len(per_genome),
        "table": str(OUT_TSV.relative_to(ROOT)),
        "totals": totals,
        "per_genome": per_genome,
    }, indent=2) + "\n")

    print(f"{len(all_rows):,} features from {len(per_genome)} genomes")
    for tool, t in totals.items():
        print(f"  {tool:7} {t['features']:>7,} features  {t['cds']:>7,} CDS")
    print(f"wrote {OUT_TSV.relative_to(ROOT)}")
    print(f"wrote {METRICS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

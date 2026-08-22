#!/usr/bin/env python3
"""Download genome assemblies and record accessions + checksums.

Reads:   data/species.tsv          curated panel, one species per line
Writes:  data/accessions.tsv       resolved accession per species (committed)
         data/genomes/<acc>.fna    sequence (gitignored, reproducible)
         results/metrics/01_genomes.json

Two stages, deliberately separable:

  resolve   species name -> a specific assembly accession, via NCBI
  fetch     accession -> bytes on disk, with a checksum

Resolution is done once and pinned into data/accessions.tsv. Every later run
reads the pins and never asks NCBI what "the reference" is again, because that
answer changes over time: NCBI can re-designate a species' reference genome,
and a panel that silently drifts is not a panel. Delete accessions.tsv to
re-resolve on purpose.

Only the genome FASTA is downloaded. NCBI's own annotation is not fetched --
this experiment compares Bakta against Prokka, and having a third annotation
sitting in the working tree invites accidentally treating it as truth.

    python3 scripts/01_fetch_genomes.py [--resolve-only] [--force]
"""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SPECIES_TSV = ROOT / "data" / "species.tsv"
ACCESSIONS_TSV = ROOT / "data" / "accessions.tsv"
GENOME_DIR = ROOT / "data" / "genomes"
METRICS = ROOT / "results" / "metrics" / "01_genomes.json"

ACC_COLUMNS = [
    "species",
    "phylum",
    "expected_gc",
    "accession",
    "assembly_name",
    "assembly_level",
    "taxid",
    "resolved_by",
]


def run(cmd):
    """Run a command, return stdout. stderr is left visible for diagnosis."""
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(cmd)}")
    return proc.stdout


def read_species():
    rows = []
    for line in SPECIES_TSV.read_text().splitlines():
        line = line.rstrip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) != 3:
            raise SystemExit(f"malformed line in species.tsv: {line!r}")
        rows.append({
            "species": parts[0],
            "phylum": parts[1],
            "expected_gc": parts[2],
        })
    if not rows:
        raise SystemExit("species.tsv contains no species")
    return rows


def summarise(extra):
    """One JSON report per line from `datasets summary genome taxon ...`."""
    out = run(["datasets", "summary", "genome", "taxon", *extra, "--as-json-lines"])
    reports = []
    for line in out.splitlines():
        line = line.strip()
        if line.startswith("{"):
            reports.append(json.loads(line))
    return reports


def contig_count(report):
    stats = report.get("assembly_stats", {}) or {}
    for key in ("number_of_contigs", "number_of_scaffolds"):
        if stats.get(key):
            return int(stats[key])
    return 10**9


def resolve(species):
    """Pick one complete assembly for a species.

    Prefers NCBI's designated reference, but the panel is complete genomes
    only, and a designated reference is not always one. Assembly level is not
    cosmetic here: both annotators stop an ORF at a contig boundary, so a
    scaffold-level genome carries a fragmentation-driven disagreement rate
    that no complete genome in the panel shares. Letting one in would put a
    per-genome artefact into a per-genome grouped split.

    Which path produced the accession is recorded, because "NCBI's reference"
    and "a complete genome we picked" are different provenance claims.
    """
    reports = summarise([species, "--reference"])
    resolved_by = "ncbi_reference"

    if reports:
        level = (reports[0].get("assembly_info", {}) or {}).get("assembly_level", "")
        if level != "Complete Genome":
            reports = []
            resolved_by = "reference_not_complete"

    if not reports:
        reports = summarise([
            species,
            "--assembly-level", "complete",
            "--assembly-source", "RefSeq",
            "--limit", "40",
        ])
        if resolved_by == "ncbi_reference":
            resolved_by = "no_reference_designated"
        # Deterministic: current assemblies first, then fewest contigs, then
        # accession as a stable tie-break. No "whatever NCBI listed first".
        reports.sort(key=lambda r: (
            (r.get("assembly_info", {}).get("assembly_status") != "current"),
            contig_count(r),
            r.get("accession", ""),
        ))

    if not reports:
        raise RuntimeError(f"no complete assembly found for {species!r}")

    r = reports[0]
    info = r.get("assembly_info", {}) or {}
    org = r.get("organism", {}) or {}
    level = info.get("assembly_level", "")
    if level != "Complete Genome":
        raise RuntimeError(
            f"{species!r}: best candidate {r['accession']} is {level!r}, "
            "not a complete genome"
        )
    return {
        "accession": r["accession"],
        "assembly_name": info.get("assembly_name", ""),
        "assembly_level": level,
        "taxid": str(org.get("tax_id", "")),
        "resolved_by": resolved_by,
    }


def write_accessions(rows):
    lines = [
        "# Pinned assembly accessions. Generated by scripts/01_fetch_genomes.py.",
        "# Delete this file to re-resolve against NCBI; otherwise it is the",
        "# authority and NCBI is not consulted for which genome to use.",
        "#",
        "# The panel is complete genomes only. resolved_by records how each",
        "# accession was arrived at:",
        "#   ncbi_reference          NCBI's designated reference, complete",
        "#   reference_not_complete  designated reference was below complete",
        "#                           assembly level and was rejected",
        "#   no_reference_designated species has no designated reference",
        "# The latter two fall back to the current complete RefSeq assembly",
        "# with the fewest contigs, accession as tie-break.",
        "#",
        "\t".join(ACC_COLUMNS),
    ]
    for row in rows:
        lines.append("\t".join(str(row[c]) for c in ACC_COLUMNS))
    ACCESSIONS_TSV.write_text("\n".join(lines) + "\n")


def read_accessions():
    if not ACCESSIONS_TSV.exists() or not ACCESSIONS_TSV.read_text().strip():
        return None
    rows = []
    for line in ACCESSIONS_TSV.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if parts[0] == "species":
            continue
        rows.append(dict(zip(ACC_COLUMNS, parts)))
    return rows or None


def download(accession, dest):
    """Fetch one assembly's genomic FASTA."""
    with tempfile.TemporaryDirectory(dir=str(GENOME_DIR)) as tmp:
        zip_path = Path(tmp) / "d.zip"
        run([
            "datasets", "download", "genome", "accession", accession,
            "--include", "genome",
            "--filename", str(zip_path),
            "--no-progressbar",
        ])
        with zipfile.ZipFile(zip_path) as zf:
            members = [
                n for n in zf.namelist()
                if n.endswith((".fna", ".fasta")) and "/data/" in n
            ]
            if not members:
                raise RuntimeError(f"no FASTA inside download for {accession}")
            # Some assemblies ship chromosome and plasmids as separate files.
            # Concatenate in sorted order so the result is byte-reproducible.
            with open(dest, "wb") as out:
                for name in sorted(members):
                    with zf.open(name) as src:
                        shutil.copyfileobj(src, out)


def sequence_stats(path):
    """Length, contig count, GC and checksum.

    Counted with str.count rather than a per-base loop: the panel is ~100 Mbp
    and a Python-level loop over every base is minutes of nothing.
    """
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)

    a = c = g = t = total = contigs = 0
    with open(path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                contigs += 1
                continue
            seq = line.strip().lower()
            a += seq.count("a")
            c += seq.count("c")
            g += seq.count("g")
            t += seq.count("t")
            total += len(seq)

    acgt = a + c + g + t
    gc = (c + g) / acgt * 100 if acgt else 0.0
    return {
        "sha256": h.hexdigest(),
        "bytes": path.stat().st_size,
        "contigs": contigs,
        "length_bp": total,
        "ambiguous_bp": total - acgt,
        "observed_gc": round(gc, 2),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--resolve-only", action="store_true",
                    help="resolve accessions and stop before downloading")
    ap.add_argument("--force", action="store_true",
                    help="re-download genomes that are already present")
    args = ap.parse_args()

    GENOME_DIR.mkdir(parents=True, exist_ok=True)
    METRICS.parent.mkdir(parents=True, exist_ok=True)

    rows = read_accessions()
    if rows is None:
        panel = read_species()
        print(f"resolving {len(panel)} species against NCBI", flush=True)
        rows = []
        for entry in panel:
            hit = resolve(entry["species"])
            row = {**entry, **hit}
            rows.append(row)
            print(f"  {entry['species']:<32} {hit['accession']:<18} "
                  f"{hit['assembly_level']:<16} {hit['resolved_by']}", flush=True)
        write_accessions(rows)
        print(f"pinned -> {ACCESSIONS_TSV.relative_to(ROOT)}")
    else:
        print(f"using {len(rows)} pinned accessions from "
              f"{ACCESSIONS_TSV.relative_to(ROOT)}")

    accs = [r["accession"] for r in rows]
    if len(set(accs)) != len(accs):
        dupes = sorted({a for a in accs if accs.count(a) > 1})
        raise SystemExit(f"duplicate accessions in panel: {dupes}")

    if args.resolve_only:
        return

    genomes = []
    for row in rows:
        dest = GENOME_DIR / f"{row['accession']}.fna"
        if dest.exists() and dest.stat().st_size > 0 and not args.force:
            print(f"  {row['accession']} present, skipping", flush=True)
        else:
            print(f"  {row['accession']} downloading", flush=True)
            download(row["accession"], dest)
        stats = sequence_stats(dest)
        genomes.append({**row, "path": str(dest.relative_to(ROOT)), **stats})
        print(f"    {stats['length_bp']:>10,} bp  {stats['contigs']:>3} contigs  "
              f"GC {stats['observed_gc']:>5.2f}%", flush=True)

    gcs = sorted(g["observed_gc"] for g in genomes)
    payload = {
        "step": "01_fetch_genomes",
        "n_genomes": len(genomes),
        "total_bp": sum(g["length_bp"] for g in genomes),
        "gc_min": gcs[0],
        "gc_max": gcs[-1],
        "gc_median": gcs[len(gcs) // 2],
        "resolved_by_counts": {
            m: sum(1 for g in genomes if g["resolved_by"] == m)
            for m in sorted({g["resolved_by"] for g in genomes})
        },
        "genomes": genomes,
    }
    METRICS.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"\n{len(genomes)} genomes, {payload['total_bp']:,} bp, "
          f"GC {payload['gc_min']}-{payload['gc_max']}%")
    print(f"wrote {METRICS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

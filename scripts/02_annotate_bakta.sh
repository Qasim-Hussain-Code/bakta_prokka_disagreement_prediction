#!/usr/bin/env bash
# Run Bakta over every genome. Record tool and DB version.
#
# The database release is part of the result, not metadata: the same genome
# against a newer Bakta database produces different calls. Tool version, DB
# version and DOI all land in results/metrics/02_annotation_versions.json,
# which is the file any published number about Bakta's calls must be read
# against.
#
# Neither annotator is given --genus/--species; see 03_annotate_prokka.sh for
# why. Both run at translation table 11 with defaults otherwise.
#
# --skip-plot only suppresses the cosmetic circular genome figure. It changes
# nothing about which features are called.
#
# --keep-contig-headers is load-bearing, not cosmetic. Bakta renames sequences
# to contig_1, contig_2, ... by default, while Prokka keeps the accession from
# the FASTA header. Two tables of intervals keyed on different sequence names
# cannot be overlapped, and the failure is silent: every region looks like a
# disagreement because no interval ever shares a seqid. The flag makes both
# tools emit the same NC_/NZ_ identifiers the input uses.
#
# Writes: data/annotations/bakta/<accession>/   (gitignored)
#         results/metrics/02_bakta_runs.json
#         results/metrics/02_annotation_versions.json

set -euo pipefail
source "$(dirname "$0")/lib_env.sh"
cd "$BPDP_ROOT"

THREADS="${BPDP_THREADS:-8}"
OUT_ROOT="data/annotations/bakta"
ACCESSIONS="data/accessions.tsv"
RUNS_JSON="results/metrics/02_bakta_runs.json"
VERSIONS_JSON="results/metrics/02_annotation_versions.json"

[[ -f "$ACCESSIONS" ]] || { echo "FATAL: $ACCESSIONS absent; run 01 first" >&2; exit 1; }
if [[ ! -f "$BAKTA_DB/version.json" ]]; then
  echo "FATAL: no Bakta database at $BAKTA_DB" >&2
  echo "       fetch it with: bash provenance/fetch_bakta_db.sh" >&2
  exit 1
fi

mkdir -p "$OUT_ROOT" results/metrics
: > "$OUT_ROOT/.runs.tsv"

mapfile -t accs < <(grep -v '^#' "$ACCESSIONS" | awk -F'\t' 'NR>1 {print $4}' | sort)
echo "bakta ${THREADS} threads over ${#accs[@]} genomes, db $BAKTA_DB"

for acc in "${accs[@]}"; do
  genome="data/genomes/${acc}.fna"
  out="$OUT_ROOT/$acc"
  gff="$out/${acc}.gff3"

  [[ -f "$genome" ]] || { echo "FATAL: $genome absent" >&2; exit 1; }

  if [[ -s "$gff" ]]; then
    echo "  $acc present, skipping"
  else
    echo "  $acc annotating"
    start=$(date +%s)
    run_bakta --db "$BAKTA_DB" --output "$out" --prefix "$acc" \
              --threads "$THREADS" --skip-plot --keep-contig-headers --force \
              "$genome" > "$OUT_ROOT/${acc}.log" 2>&1 \
      || { echo "FATAL: bakta failed on $acc; see $OUT_ROOT/${acc}.log" >&2
           tail -20 "$OUT_ROOT/${acc}.log" >&2; exit 1; }
    echo "$acc	$(( $(date +%s) - start ))" >> "$OUT_ROOT/.runs.tsv"
  fi

  [[ -s "$gff" ]] || { echo "FATAL: $gff missing after run" >&2; exit 1; }
  printf '    %s CDS\n' "$(grep -cP '\tCDS\t' "$gff" || true)"
done

echo "=== recording bakta run provenance ==="
python3 - "$OUT_ROOT" "$RUNS_JSON" "$VERSIONS_JSON" "$THREADS" "$BAKTA_DB" <<'PY'
import hashlib, json, subprocess, sys
from pathlib import Path

out_root, runs_json, versions_json = (Path(a) for a in sys.argv[1:4])
threads, bakta_db = int(sys.argv[4]), Path(sys.argv[5])

seconds = {}
runs_tsv = out_root / ".runs.tsv"
if runs_tsv.exists():
    for line in runs_tsv.read_text().splitlines():
        if line.strip():
            acc, sec = line.split("\t")
            seconds[acc] = int(sec)

records = []
for d in sorted(p for p in out_root.iterdir() if p.is_dir()):
    gff = d / f"{d.name}.gff3"
    if not gff.exists():
        continue
    # Bakta appends the input sequence after ##FASTA, as Prokka does.
    body = gff.read_text().split("\n##FASTA", 1)[0]
    kinds = {}
    for ln in body.splitlines():
        if ln and not ln.startswith("#") and ln.count("\t") >= 8:
            kind = ln.split("\t")[2]
            kinds[kind] = kinds.get(kind, 0) + 1
    records.append({
        "accession": d.name,
        "gff": str(gff),
        "sha256": hashlib.sha256(gff.read_bytes()).hexdigest(),
        "feature_counts": dict(sorted(kinds.items())),
        "cds": kinds.get("CDS", 0),
        "seconds": seconds.get(d.name),
    })

runs_json.write_text(json.dumps({
    "step": "02_annotate_bakta",
    "threads": threads,
    "n_genomes": len(records),
    "total_cds": sum(r["cds"] for r in records),
    "runs": records,
}, indent=2) + "\n")
print(f"wrote {runs_json}: {len(records)} genomes, "
      f"{sum(r['cds'] for r in records):,} CDS")


def tool_version(fn):
    out = subprocess.run(
        ["bash", "-c", f"source scripts/lib_env.sh && {fn} --version 2>&1 | head -1"],
        capture_output=True, text=True)
    return out.stdout.strip()


db_files = sorted(p for p in bakta_db.rglob("*") if p.is_file())

db_version = json.loads((bakta_db / "version.json").read_text())

# AMRFinderPlus is versioned separately from the Bakta database that ships it.
# AMRFinderPlus 4.x refuses to run against a database older than a version it
# hardcodes, and the release bundled in bakta db-light v6.0 is older than the
# binary available here demands, so it had to be updated in place. The result
# is that the AMR annotations come from a newer reference set than the Bakta
# release declares. That affects annotation *content* -- which CDS get an AMR
# gene name -- and not which regions are called as CDS, so it does not move
# the label. Recorded rather than reconciled, because a reader comparing this
# run against a stock db-light v6.0 needs to know.
amr_latest = bakta_db / "amrfinderplus-db" / "latest"
amr_used = None
if (amr_latest / "version.txt").exists():
    amr_used = (amr_latest / "version.txt").read_text().strip()
amr_bundled = next(
    (d.get("release") for d in db_version.get("dependencies", [])
     if d.get("name") == "AMRFinderPlus"),
    None,
)
versions_json.write_text(json.dumps({
    "step": "02_annotation_versions",
    "note": (
        "The database release determines which calls Bakta makes and is part "
        "of the result, not metadata. Any published Bakta number must be read "
        "against this file."
    ),
    "bakta": {
        "tool": tool_version("run_bakta"),
        "db_path": str(bakta_db),
        "db_type": "light",
        "db_type_reason": (
            "db-full does not fit in available disk on this machine. "
            "db-light has fewer reference proteins, so some CDS that db-full "
            "would name are left hypothetical. This affects annotation "
            "content; it does not change which regions are called as CDS."
        ),
        "db_version": db_version,
        "amrfinderplus_db_version_used": amr_used,
        "amrfinderplus_db_version_bundled": amr_bundled,
        "amrfinderplus_db_note": (
            "AMRFinderPlus 4.2.7 refuses any database older than a version it "
            "hardcodes, and the release bundled in bakta db-light v6.0 was "
            "older, so it was updated in place. AMR annotation content "
            "therefore comes from a newer reference set than the Bakta "
            "release declares. It does not change which regions are called "
            "as CDS, so the label is unaffected."
        ),
        "db_bytes": sum(p.stat().st_size for p in db_files),
        "db_file_count": len(db_files),
    },
    "prokka": {"tool": tool_version("run_prokka")},
    "shared_parameters": {
        "organism_hints": None,
        "organism_hints_note": (
            "Neither tool given --genus/--species. Both change output when "
            "hinted; hinting one and not the other would make the "
            "disagreement partly a function of the hint."
        ),
        "translation_table": 11,
        "kingdom": "Bacteria",
        "contig_headers": "original",
        "contig_headers_note": (
            "Bakta run with --keep-contig-headers so both tools key intervals "
            "on the same sequence identifiers. Without it Bakta emits "
            "contig_N and no interval could ever be matched to Prokka's."
        ),
    },
}, indent=2) + "\n")
print(f"wrote {versions_json}")
PY

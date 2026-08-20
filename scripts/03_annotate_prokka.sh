#!/usr/bin/env bash
# Run Prokka over the same genomes, same order.
#
# Neither annotator is given --genus/--species. Both tools accept organism
# hints and both change their output when given them. Supplying hints to one
# and not the other, or supplying them unevenly across a 25-genome panel,
# would make the disagreement partly a function of what we told the tools
# rather than what they inferred. Defaults on both sides, kingdom Bacteria,
# translation table 11 -- the way a user runs these out of the box.
#
# Writes: data/annotations/prokka/<accession>/   (gitignored)
#         results/metrics/03_prokka_runs.json

set -euo pipefail
source "$(dirname "$0")/lib_env.sh"
cd "$BPDP_ROOT"

THREADS="${BPDP_THREADS:-8}"
OUT_ROOT="data/annotations/prokka"
ACCESSIONS="data/accessions.tsv"
RUNS_JSON="results/metrics/03_prokka_runs.json"

[[ -f "$ACCESSIONS" ]] || { echo "FATAL: $ACCESSIONS absent; run 01 first" >&2; exit 1; }

mkdir -p "$OUT_ROOT" results/metrics
: > "$OUT_ROOT/.runs.tsv"

# Deterministic order: accession, not directory listing order.
mapfile -t accs < <(grep -v '^#' "$ACCESSIONS" | awk -F'\t' 'NR>1 {print $4}' | sort)
echo "prokka ${THREADS} threads over ${#accs[@]} genomes"

for acc in "${accs[@]}"; do
  genome="data/genomes/${acc}.fna"
  out="$OUT_ROOT/$acc"
  gff="$out/${acc}.gff"

  [[ -f "$genome" ]] || { echo "FATAL: $genome absent" >&2; exit 1; }

  if [[ -s "$gff" ]]; then
    echo "  $acc present, skipping"
  else
    echo "  $acc annotating"
    start=$(date +%s)
    # Prokka refuses to write into a non-empty directory without --force.
    run_prokka --outdir "$out" --prefix "$acc" --cpus "$THREADS" --force \
               "$genome" > "$OUT_ROOT/${acc}.log" 2>&1 \
      || { echo "FATAL: prokka failed on $acc; see $OUT_ROOT/${acc}.log" >&2
           tail -20 "$OUT_ROOT/${acc}.log" >&2; exit 1; }
    echo "$acc	$(( $(date +%s) - start ))" >> "$OUT_ROOT/.runs.tsv"
  fi

  [[ -s "$gff" ]] || { echo "FATAL: $gff missing after run" >&2; exit 1; }
  printf '    %s CDS\n' "$(grep -cP '\tCDS\t' "$gff" || true)"
done

echo "=== recording prokka run provenance ==="
python3 - "$OUT_ROOT" "$RUNS_JSON" "$THREADS" <<'PY'
import hashlib, json, subprocess, sys
from pathlib import Path

out_root, runs_json, threads = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])

seconds = {}
runs_tsv = out_root / ".runs.tsv"
if runs_tsv.exists():
    for line in runs_tsv.read_text().splitlines():
        if line.strip():
            acc, sec = line.split("\t")
            seconds[acc] = int(sec)

records = []
for d in sorted(p for p in out_root.iterdir() if p.is_dir()):
    gff = d / f"{d.name}.gff"
    if not gff.exists():
        continue
    text = gff.read_text()
    # Prokka appends the input sequence after ##FASTA. Feature counting must
    # stop there or the sequence lines get scanned as annotation.
    body = text.split("\n##FASTA", 1)[0]
    feats = [ln.split("\t") for ln in body.splitlines()
             if ln and not ln.startswith("#") and ln.count("\t") >= 8]
    kinds = {}
    for f in feats:
        kinds[f[2]] = kinds.get(f[2], 0) + 1
    records.append({
        "accession": d.name,
        "gff": str(gff),
        "sha256": hashlib.sha256(gff.read_bytes()).hexdigest(),
        "feature_counts": dict(sorted(kinds.items())),
        "cds": kinds.get("CDS", 0),
        "seconds": seconds.get(d.name),
    })

version = subprocess.run(["bash", "-c",
    'source scripts/lib_env.sh && run_prokka --version 2>&1 | head -1'],
    capture_output=True, text=True).stdout.strip()

runs_json.write_text(json.dumps({
    "step": "03_annotate_prokka",
    "tool": version,
    "threads": threads,
    "organism_hints": None,
    "organism_hints_note": (
        "No --genus/--species given, matching Bakta. Organism hints change "
        "both tools' output; supplying them would make disagreement partly a "
        "function of the hint rather than of the sequence."
    ),
    "kingdom": "Bacteria",
    "translation_table": 11,
    "n_genomes": len(records),
    "total_cds": sum(r["cds"] for r in records),
    "runs": records,
}, indent=2) + "\n")
print(f"wrote {runs_json}: {len(records)} genomes, "
      f"{sum(r['cds'] for r in records):,} CDS")
PY

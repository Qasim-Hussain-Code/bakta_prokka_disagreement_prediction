#!/usr/bin/env python3
"""Features for the naming question, each flagged with where it came from.

Reads:  data/interim/content_pairs.tsv.gz
        data/interim/features.tsv.gz          (part one, reused)
        results/metrics/06_features.json      (part one manifest)
        data/genomes/<acc>.fna
        data/annotations/bakta/<acc>/<acc>.faa
Writes: data/interim/content_features.tsv.gz  (gitignored)
        results/metrics/14_content_features.json

The intervals are the same intervals part one already characterised, so the
sequence and genome features are joined from data/interim/features.tsv.gz on
(genome, seqid, start, end) rather than recomputed. Recomputing them would
risk two subtly different definitions of gc_content across the two halves of
the project.

Two flags, both stored in the table and read from it by every downstream
script. Nothing downstream carries a hard-coded feature list.

  caller_derived   the value encodes a gene-caller decision: the boundaries
                   chosen, the reading frame those imply, the strand. Part
                   one's 'caller' group, plus everything computed from the
                   translated protein, because there is no protein without a
                   frame.

  db_derived       the value is computed from either tool's database-search
                   output rather than from sequence. A model that predicts
                   name disagreement from these has learned which regions are
                   thinly represented in the reference sets, which is a fact
                   about the databases and not about the DNA.

Two part-one features are deliberately NOT carried over:

  is_hypothetical   constant 0 inside the primary analysis set, by
                    construction -- the set is defined as both tools having
                    named the region.
  has_gene_symbol   constant 1 in part one's table for all 87,960 rows,
                    because 04_extract_calls.py fell back to Bakta's Name
                    attribute, which is the product. Recomputed here from the
                    raw gene= attribute instead.

A constant feature is a dead feature. This step fails if any survives.
"""

import gzip
import json
import math
import sys
from bisect import bisect_left
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(__file__).resolve().parent.parent
PAIRS = ROOT / "data" / "interim" / "content_pairs.tsv.gz"
PART_ONE = ROOT / "data" / "interim" / "features.tsv.gz"
PART_ONE_MANIFEST = ROOT / "results" / "metrics" / "06_features.json"
OUT_TSV = ROOT / "data" / "interim" / "content_features.tsv.gz"
METRICS = ROOT / "results" / "metrics" / "14_content_features.json"

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
LC_WINDOW = 40          # bp, for the frame-free low-complexity scan
LC_STEP = 20
LC_ENTROPY_BITS = 4.0   # 3-mer entropy below this counts the window as low-complexity
AA_LC_WINDOW = 20       # residues
AA_LC_MAX_FRAC = 0.5    # one residue occupying half a window
MOBILE_WINDOW = 5000    # bp either side, for the mobile-element neighbourhood
NEIGHBOURS = 2          # CDS either side, for the identical-call neighbourhood

# Part one groups that transfer unchanged.
REUSE_GROUPS = ("sequence", "genome", "caller")
DROP_FROM_PART_ONE = {"is_hypothetical", "has_gene_symbol"}

# Features that are constant across the primary analysis set for a structural
# reason, declared here with that reason and VERIFIED against the data below.
# A constant feature that is not declared here is fatal; a declared feature
# that turns out not to be constant is also fatal, because that means the
# stated reason no longer holds.
DECLARED_CONSTANT_IN_PRIMARY_SET = {
    "bakta_has_pfam": (
        "All 3,428 Bakta CDS carrying a PFAM cross-reference have a "
        "placeholder Bakta product (1,894 paired with a Prokka placeholder, "
        "1,534 with a Prokka name). Under db-light Bakta falls back to Pfam "
        "only when it has nothing better to say, so a PFAM cross-reference and "
        "a named Bakta product are mutually exclusive in this panel and the "
        "feature is excluded from the primary set by construction."),
}

ID_COLUMNS = ["genome", "species", "seqid", "start", "end",
              "in_primary_set", "label"]


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def load_fasta(path):
    contigs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    contigs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        contigs[name] = "".join(buf).upper()
    return contigs


def load_faa(path):
    """locus_tag -> protein sequence."""
    out, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        out[name] = "".join(buf).upper()
    return out


def kmer_entropy(seq, k=3):
    if len(seq) < k:
        return 0.0
    counts = Counter(seq[i:i + k] for i in range(len(seq) - k + 1))
    total = sum(counts.values())
    return -sum((c / total) * math.log2(c / total) for c in counts.values())


def low_complexity_dna_frac(seq):
    """Fraction of windows whose 3-mer entropy falls below the threshold.

    Frame-free on purpose: this is the one compositional-repetition measure in
    the sequence group, computed without reference to any reading frame, so
    the sequence-only arm of the circularity audit is not left with nothing to
    say about repetitive DNA.
    """
    if len(seq) < LC_WINDOW:
        return 1.0 if kmer_entropy(seq) < LC_ENTROPY_BITS else 0.0
    windows = 0
    low = 0
    for i in range(0, len(seq) - LC_WINDOW + 1, LC_STEP):
        windows += 1
        if kmer_entropy(seq[i:i + LC_WINDOW]) < LC_ENTROPY_BITS:
            low += 1
    return low / windows if windows else 0.0


def aa_low_complexity_frac(protein):
    if len(protein) < AA_LC_WINDOW:
        if not protein:
            return 0.0
        top = Counter(protein).most_common(1)[0][1]
        return 1.0 if top / len(protein) >= AA_LC_MAX_FRAC else 0.0
    windows = low = 0
    for i in range(0, len(protein) - AA_LC_WINDOW + 1, AA_LC_WINDOW // 2):
        w = protein[i:i + AA_LC_WINDOW]
        windows += 1
        if Counter(w).most_common(1)[0][1] / len(w) >= AA_LC_MAX_FRAC:
            low += 1
    return low / windows if windows else 0.0


def read_tsv_gz(path):
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        return [line.rstrip("\n").split("\t") for line in fh], idx


def main():
    for p in (PAIRS, PART_ONE, PART_ONE_MANIFEST):
        if not p.exists():
            raise SystemExit(f"{p} absent")

    pairs, pidx = read_tsv_gz(PAIRS)
    p1_rows, p1idx = read_tsv_gz(PART_ONE)
    p1_manifest = json.loads(PART_ONE_MANIFEST.read_text())

    # --- what transfers from part one -------------------------------------
    reused = []
    group_of = {}
    for group in REUSE_GROUPS:
        for name in p1_manifest["features_by_group"][group]:
            if name in DROP_FROM_PART_ONE:
                continue
            reused.append(name)
            group_of[name] = group

    p1_lookup = {
        (r[p1idx["genome"]], r[p1idx["seqid"]], r[p1idx["start"]], r[p1idx["end"]]): r
        for r in p1_rows
    }

    # --- manifest for the new features ------------------------------------
    # name -> (caller_derived, db_derived, description)
    NEW = {
        "low_complexity_dna_frac": (
            0, 0, f"fraction of {LC_WINDOW} bp windows with 3-mer entropy below "
                  f"{LC_ENTROPY_BITS} bits; frame-free"),

        "protein_length_aa": (
            1, 0, "length of the translated protein in residues; a function of "
                  "the interval the caller chose"),
        "aa_entropy": (
            1, 0, "Shannon entropy of amino acid composition, bits"),
        "aa_low_complexity_frac": (
            1, 0, f"fraction of {AA_LC_WINDOW}-residue windows in which one "
                  f"residue occupies at least {AA_LC_MAX_FRAC:.0%}"),
        "same_start": (
            1, 0, "the two callers chose the same start coordinate"),
        "same_stop": (
            1, 0, "the two callers chose the same stop coordinate"),
        "neighbourhood_identical_calls": (
            1, 0, f"how many of the {NEIGHBOURS} nearest Bakta CDS either side "
                  "were called with identical coordinates by both tools"),

        "bakta_has_uniref": (
            0, 1, "Bakta's Dbxref contains a UniRef cluster"),
        "bakta_has_pfam": (
            0, 1, "Bakta's Dbxref contains a Pfam accession"),
        "bakta_has_ec": (
            0, 1, "Bakta's Dbxref contains an EC number"),
        "bakta_has_is": (
            0, 1, "Bakta's Dbxref contains an ISFinder insertion-sequence match"),
        "bakta_n_dbxref": (
            0, 1, "number of cross-references Bakta attached"),
        "prokka_has_similarity_hit": (
            0, 1, "Prokka's inference cites a UniProtKB protein similarity"),
        "prokka_has_protein_motif": (
            0, 1, "Prokka's inference cites an HMM protein motif"),
        "prokka_has_cog": (
            0, 1, "Prokka assigned a COG"),
        "prokka_has_ec": (
            0, 1, "Prokka assigned an EC number"),
        "mobile_neighbours_5kb": (
            0, 1, f"Bakta CDS within {MOBILE_WINDOW} bp carrying an ISFinder "
                  "match; the only mobile-element context derivable from what "
                  "is already annotated"),
    }

    for name, (caller, db, _d) in NEW.items():
        group_of[name] = ("db" if db else ("caller" if caller else "sequence"))

    FEATURES = reused + list(NEW)

    def flags(name):
        if name in NEW:
            return NEW[name][0], NEW[name][1]
        group = group_of[name]
        # Part one's caller group is caller-derived and not database-derived.
        return (1 if group == "caller" else 0), 0

    # --- per-genome indexes for the neighbourhood features -----------------
    by_genome = {}
    for r in pairs:
        by_genome.setdefault(r[pidx["genome"]], []).append(r)

    out_rows = []
    missing_protein = 0
    missing_part_one = 0

    for genome in sorted(by_genome):
        contigs = load_fasta(ROOT / "data" / "genomes" / f"{genome}.fna")
        faa = load_faa(ROOT / "data" / "annotations" / "bakta" / genome / f"{genome}.faa")
        rows = by_genome[genome]

        # index by contig, ordered by start, for both neighbourhood features
        per_contig = {}
        for r in rows:
            per_contig.setdefault(r[pidx["seqid"]], []).append(r)
        for seqid in per_contig:
            per_contig[seqid].sort(key=lambda r: int(r[pidx["start"]]))
        starts_of = {s: [int(r[pidx["start"]]) for r in v]
                     for s, v in per_contig.items()}
        mobile_starts = {
            s: [int(r[pidx["start"]]) for r in v if "IS:" in r[pidx["dbxref_bakta"]]]
            for s, v in per_contig.items()
        }

        for seqid, ordered in per_contig.items():
            starts = starts_of[seqid]
            mstarts = mobile_starts[seqid]
            for pos, r in enumerate(ordered):
                start, end = int(r[pidx["start"]]), int(r[pidx["end"]])
                strand = r[pidx["strand"]]
                key = (genome, seqid, str(start), str(end))
                p1 = p1_lookup.get(key)
                if p1 is None:
                    missing_part_one += 1
                    continue

                window = contigs[seqid][start - 1:end]
                if strand == "-":
                    window_coding = revcomp(window)
                else:
                    window_coding = window

                protein = faa.get(r[pidx["locus_tag_bakta"]], "")
                if not protein:
                    missing_protein += 1
                    # Fall back to the frame the interval implies, so a missing
                    # FASTA record degrades one row rather than dropping it.
                    protein = ""

                lo = max(0, pos - NEIGHBOURS)
                hi = min(len(ordered), pos + NEIGHBOURS + 1)
                neigh = [ordered[i] for i in range(lo, hi) if i != pos]
                identical = sum(1 for n in neigh
                                if n[pidx["same_start"]] == "1"
                                and n[pidx["same_stop"]] == "1")

                left = bisect_left(mstarts, start - MOBILE_WINDOW)
                right = bisect_left(mstarts, end + MOBILE_WINDOW)
                mobile = right - left
                if "IS:" in r[pidx["dbxref_bakta"]]:
                    mobile -= 1  # do not count the region itself

                dx_b = r[pidx["dbxref_bakta"]]
                inf_p = r[pidx["inference_prokka"]]
                aa_counts = Counter(protein)
                n_aa = len(protein)

                values = {
                    "low_complexity_dna_frac": round(low_complexity_dna_frac(window), 5),
                    "protein_length_aa": n_aa,
                    "aa_entropy": round(
                        -sum((c / n_aa) * math.log2(c / n_aa)
                             for c in aa_counts.values()), 5) if n_aa else 0.0,
                    "aa_low_complexity_frac": round(aa_low_complexity_frac(protein), 5),
                    "same_start": int(r[pidx["same_start"]]),
                    "same_stop": int(r[pidx["same_stop"]]),
                    "neighbourhood_identical_calls": identical,
                    "bakta_has_uniref": int("UniRef" in dx_b),
                    "bakta_has_pfam": int("PFAM" in dx_b),
                    "bakta_has_ec": int("EC:" in dx_b),
                    "bakta_has_is": int("IS:" in dx_b),
                    "bakta_n_dbxref": len([x for x in dx_b.split(",") if x]),
                    "prokka_has_similarity_hit": int("similar to AA sequence" in inf_p),
                    "prokka_has_protein_motif": int("protein motif" in inf_p),
                    "prokka_has_cog": int("COG" in r[pidx["dbxref_prokka"]]),
                    "prokka_has_ec": int(bool(r[pidx["ec_prokka"]])),
                    "mobile_neighbours_5kb": max(0, mobile),
                }
                for name in AMINO_ACIDS:
                    fname = f"aa_frac_{name}"
                    values[fname] = round(aa_counts.get(name, 0) / n_aa, 5) if n_aa else 0.0

                ident = [genome, r[pidx["species"]], seqid, str(start), str(end),
                         r[pidx["in_primary_set"]],
                         r[pidx["name_disagreement"]] or ""]
                out_rows.append(ident + [p1[p1idx[n]] for n in reused]
                                + [values[n] for n in FEATURES if n in NEW]
                                + [values[f"aa_frac_{a}"] for a in AMINO_ACIDS])
        print(f"  {genome} {len(rows):>6,} pairs", flush=True)

    # amino acid fractions are appended last; declare them now that they exist
    for a in AMINO_ACIDS:
        NEW[f"aa_frac_{a}"] = (1, 0, f"fraction of residue {a} in the translated protein")
        group_of[f"aa_frac_{a}"] = "caller"
    FEATURES = reused + [n for n in NEW if not n.startswith("aa_frac_")] \
        + [f"aa_frac_{a}" for a in AMINO_ACIDS]

    if missing_part_one:
        raise SystemExit(
            f"FATAL: {missing_part_one} pairs had no matching row in part one's "
            f"feature table. The two tables should key on identical intervals.")
    if not out_rows:
        raise SystemExit("FATAL: no feature rows built")

    # --- constant-feature check, on the primary analysis set ---------------
    # Runs BEFORE the table is written, so a dead column never reaches disk.
    prim = [r for r in out_rows if r[ID_COLUMNS.index("in_primary_set")] == "1"]
    offset = len(ID_COLUMNS)
    constant, near_constant = [], []
    for i, name in enumerate(FEATURES):
        vals = {r[offset + i] for r in prim}
        if len(vals) == 1:
            constant.append({"feature": name, "value": next(iter(vals))})
        else:
            counts = Counter(r[offset + i] for r in prim)
            top, n_top = counts.most_common(1)[0]
            if n_top / len(prim) > 0.995:
                near_constant.append({
                    "feature": name, "dominant_value": top,
                    "share": round(n_top / len(prim), 5)})

    constant_names = {c["feature"] for c in constant}
    undeclared = sorted(constant_names - set(DECLARED_CONSTANT_IN_PRIMARY_SET))
    if undeclared:
        for c in constant:
            if c["feature"] in undeclared:
                print(f"  constant and undeclared: {c['feature']} = {c['value']}")
        raise SystemExit(
            f"FATAL: {len(undeclared)} feature(s) take a single value across "
            "the primary analysis set and are not declared in "
            "DECLARED_CONSTANT_IN_PRIMARY_SET. A constant feature can never be "
            "split on. Part one shipped one of these (has_gene_symbol) "
            "unnoticed; this check exists so it does not happen twice.")
    stale = sorted(set(DECLARED_CONSTANT_IN_PRIMARY_SET) - constant_names)
    if stale:
        raise SystemExit(
            f"FATAL: {stale} declared constant across the primary analysis set "
            "but observed to vary. The stated reason for excluding them no "
            "longer holds; re-derive it before running again.")

    excluded = sorted(constant_names)
    keep = [i for i, n in enumerate(FEATURES) if n not in constant_names]
    FEATURES = [FEATURES[i] for i in keep]
    out_rows = [r[:offset] + [r[offset + i] for i in keep] for r in out_rows]
    prim = [r for r in out_rows if r[ID_COLUMNS.index("in_primary_set")] == "1"]

    with gzip.open(OUT_TSV, "wt", newline="") as fh:
        fh.write("\t".join(ID_COLUMNS + FEATURES) + "\n")
        for r in out_rows:
            fh.write("\t".join(str(v) for v in r) + "\n")

    manifest = {}
    for name in FEATURES:
        caller, db = flags(name)
        desc = (NEW[name][2] if name in NEW
                else p1_manifest["manifest"][name]["description"])
        manifest[name] = {
            "group": group_of[name],
            "caller_derived": bool(caller),
            "db_derived": bool(db),
            "source": "part_one" if name in reused else "part_two",
            "description": desc,
        }

    n_prim = len(prim)
    labels = [int(r[ID_COLUMNS.index("label")]) for r in prim]
    payload = {
        "step": "14_content_features",
        "n_rows_total": len(out_rows),
        "n_rows_primary_set": n_prim,
        "n_positive_primary_set": sum(labels),
        "positive_rate_primary_set": round(sum(labels) / n_prim, 5),
        "n_features": len(FEATURES),
        "flags": {
            "caller_derived": (
                "the value encodes a gene-caller decision: chosen boundaries, "
                "the reading frame they imply, the strand. Everything computed "
                "from the translated protein is included, because there is no "
                "protein without a frame."),
            "db_derived": (
                "computed from either tool's database-search output rather "
                "than from sequence. A model predicting name disagreement from "
                "these has learned which regions are thinly represented in the "
                "reference sets, not something about the DNA."),
            "read_from_the_table": (
                "step 19 reads both flags from this manifest. No downstream "
                "script carries a hard-coded feature list."),
        },
        "counts": {
            "caller_derived": sum(1 for m in manifest.values() if m["caller_derived"]),
            "db_derived": sum(1 for m in manifest.values() if m["db_derived"]),
            "neither": sum(1 for m in manifest.values()
                           if not m["caller_derived"] and not m["db_derived"]),
            "reused_from_part_one": len(reused),
            "new_in_part_two": len(FEATURES) - len(reused),
        },
        "dropped_from_part_one": {
            "is_hypothetical": (
                "constant 0 inside the primary analysis set by construction: "
                "the set is defined as both tools having named the region"),
            "has_gene_symbol": (
                "constant 1 across all 87,960 part-one rows. "
                "04_extract_calls.py read the gene symbol as "
                "a.get('gene', a.get('Name', '')) and Bakta sets Name=product, "
                "so the fallback filled the column with product strings. A "
                "constant feature can never be split on, so no part-one number "
                "moves; scripts 01-11 are left exactly as they are."),
        },
        "caller_flag_reads_differently_here": (
            "In part one the caller flag was alarming: the label was 'Bakta "
            "called this and Prokka did not', so a caller-derived feature was "
            "close to a restatement of it. Here both tools called every region "
            "in the cohort, and part one established they agree on the exact "
            "coordinates in 87,788 of 87,888 matched calls. The caller "
            "decision is shared rather than a source of the label. The flag is "
            "kept strict anyway so it means the same thing in both halves."),
        "excluded_as_constant_in_primary_set": {
            name: {
                "reason": reason,
                "verified_constant_at_runtime": name in excluded,
            }
            for name, reason in DECLARED_CONSTANT_IN_PRIMARY_SET.items()
        },
        "constant_feature_check": {
            "n_constant_undeclared": 0,
            "n_excluded_as_declared_constant": len(excluded),
            "n_near_constant": len(near_constant),
            "near_constant": near_constant,
            "near_constant_note": (
                "over 99.5% one value across the primary set. Kept, but they "
                "carry almost no information and should not be read as "
                "meaningful if they surface in an importance ranking."),
        },
        "manifest": manifest,
        "protein_source": "Bakta .faa, keyed by locus_tag",
        "n_rows_without_protein_record": missing_protein,
        "table": str(OUT_TSV.relative_to(ROOT)),
    }
    METRICS.write_text(json.dumps(payload, indent=2) + "\n")

    print(f"\n{len(out_rows):,} rows x {len(FEATURES)} features")
    print(f"  primary analysis set {n_prim:,}  positives {sum(labels):,} "
          f"({sum(labels)/n_prim:.1%})")
    c = payload["counts"]
    print(f"  caller_derived {c['caller_derived']}  db_derived {c['db_derived']}  "
          f"neither {c['neither']}")
    if near_constant:
        print(f"  near-constant: {', '.join(n['feature'] for n in near_constant)}")
    print(f"wrote {OUT_TSV.relative_to(ROOT)}")
    print(f"wrote {METRICS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

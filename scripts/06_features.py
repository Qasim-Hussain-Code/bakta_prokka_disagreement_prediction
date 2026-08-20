#!/usr/bin/env python3
"""Sequence features per interval. Flag which come from gene-caller logic.

Reads:  data/interim/disagreement.tsv.gz
        data/genomes/<acc>.fna
        results/metrics/01_genomes.json
Writes: data/interim/features.tsv.gz   (gitignored)
        results/metrics/06_features.json

The question is whether disagreement is predictable *from the sequence*. That
makes the provenance of each feature part of the experiment, not bookkeeping,
so every feature carries a group:

  sequence  computed from the DNA in the interval and its flanks. Knowing
            where Bakta drew the boundary is required to find the window, but
            nothing about Bakta's reasoning enters the value.

  genome    a property of the whole genome, identical for every interval in
            it. Legitimate under a grouped-by-genome split, where test
            genomes are unseen, and it lets "disagreement is higher in
            high-GC genomes" be learned as the generalisable claim it is.

  caller    derived from Bakta's own decision or its database lookup: the
            boundaries it chose, the frame that implies, whether it found a
            name for the product. These are the circular ones.

The caller group is where the interesting failure lives. `is_hypothetical` in
particular is close to a restatement of the label: a CDS Bakta could not name
is exactly the kind of call Prokka is likely not to make either. A model given
that feature can score well while having learned nothing about sequence.
Step 10 refits with the whole caller group dropped and reports both numbers.
Nothing here is removed on suspicion -- it is labelled, and the audit decides.
"""

import gzip
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DISAGREEMENT = ROOT / "data" / "interim" / "disagreement.tsv.gz"
GENOME_METRICS = ROOT / "results" / "metrics" / "01_genomes.json"
OUT_TSV = ROOT / "data" / "interim" / "features.tsv.gz"
METRICS = ROOT / "results" / "metrics" / "06_features.json"

FLANK = 100          # bp either side, for the upstream/downstream composition
EDGE_THRESHOLD = 500  # bp from a contig end to count as "near the edge"

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# name -> (group, description)
FEATURE_MANIFEST = {
    "gc_content":        ("sequence", "GC fraction of the called interval"),
    "gc_skew":           ("sequence", "(G-C)/(G+C) over the interval"),
    "at_skew":           ("sequence", "(A-T)/(A+T) over the interval"),
    "gc_minus_genome":   ("sequence", "interval GC minus whole-genome GC; a "
                                      "classic signal of horizontally acquired "
                                      "or mobile sequence"),
    "entropy_3mer":      ("sequence", "Shannon entropy of 3-mer composition, "
                                      "bits; low in repetitive sequence"),
    "max_homopolymer":   ("sequence", "longest single-base run"),
    "longest_orf_frac":  ("sequence", "longest stop-free stretch across all 6 "
                                      "frames, as a fraction of interval "
                                      "length; computed independently of the "
                                      "frame Bakta chose"),
    "min_stops_per_kb":  ("sequence", "stop codons per kb in the emptiest of "
                                      "the 6 frames"),
    "flank_gc_up":       ("sequence", f"GC of the {FLANK} bp before the interval"),
    "flank_gc_down":     ("sequence", f"GC of the {FLANK} bp after the interval"),
    "dist_to_contig_end": ("sequence", "bp to the nearer end of the contig"),
    "near_contig_edge":  ("sequence", f"within {EDGE_THRESHOLD} bp of a contig end"),
    "ambiguous_frac":    ("sequence", "fraction of non-ACGT bases in the interval"),

    "genome_gc":         ("genome", "GC of the whole genome"),
    "genome_length_bp":  ("genome", "total genome length"),
    "genome_n_contigs":  ("genome", "number of replicons/contigs"),

    "length_bp":         ("caller", "length of the interval Bakta chose"),
    "strand_plus":       ("caller", "Bakta called this on the forward strand"),
    "gc3":               ("caller", "GC at third codon positions; requires the "
                                    "reading frame Bakta chose"),
    "start_atg":         ("caller", "Bakta's chosen start codon is ATG"),
    "start_gtg":         ("caller", "Bakta's chosen start codon is GTG"),
    "start_ttg":         ("caller", "Bakta's chosen start codon is TTG"),
    "stop_taa":          ("caller", "stop codon of the called ORF is TAA"),
    "stop_tag":          ("caller", "stop codon of the called ORF is TAG"),
    "stop_tga":          ("caller", "stop codon of the called ORF is TGA"),
    "is_hypothetical":   ("caller", "Bakta's product is hypothetical/unknown; "
                                    "near-restatement of the label"),
    "has_gene_symbol":   ("caller", "Bakta assigned a gene symbol"),
}

FEATURES = list(FEATURE_MANIFEST)
ID_COLUMNS = ["genome", "species", "seqid", "start", "end", "label"]


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


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


def gc_fraction(seq):
    acgt = seq.count("A") + seq.count("C") + seq.count("G") + seq.count("T")
    return (seq.count("G") + seq.count("C")) / acgt if acgt else 0.0


def skew(seq, a, b):
    na, nb = seq.count(a), seq.count(b)
    return (na - nb) / (na + nb) if (na + nb) else 0.0


def entropy_3mer(seq):
    if len(seq) < 3:
        return 0.0
    counts = Counter(seq[i:i + 3] for i in range(len(seq) - 2))
    total = sum(counts.values())
    return -sum((c / total) * math.log2(c / total) for c in counts.values())


def max_homopolymer(seq):
    best = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        if run > best:
            best = run
    return best if seq else 0


def six_frame_orf_stats(seq):
    """Longest stop-free stretch and minimum stop density across 6 frames.

    Deliberately independent of the frame Bakta chose: a gene caller's frame
    is its conclusion, and reusing it would smuggle that conclusion in as a
    feature.
    """
    stops = {"TAA", "TAG", "TGA"}
    longest, min_density = 0, float("inf")
    for strand_seq in (seq, revcomp(seq)):
        for frame in range(3):
            since_stop, n_stops, codons = 0, 0, 0
            for i in range(frame, len(strand_seq) - 2, 3):
                codons += 1
                if strand_seq[i:i + 3] in stops:
                    n_stops += 1
                    longest = max(longest, since_stop * 3)
                    since_stop = 0
                else:
                    since_stop += 1
            longest = max(longest, since_stop * 3)
            if codons:
                min_density = min(min_density, n_stops / (codons * 3 / 1000))
    return longest, (0.0 if min_density == float("inf") else min_density)


def main():
    if not DISAGREEMENT.exists():
        raise SystemExit(f"{DISAGREEMENT} absent; run 05 first")

    genome_meta = {
        g["accession"]: g
        for g in json.loads(GENOME_METRICS.read_text())["genomes"]
    }

    with gzip.open(DISAGREEMENT, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        rows = [line.rstrip("\n").split("\t") for line in fh]

    # Group by genome so each FASTA is read once, not once per interval.
    by_genome = {}
    for r in rows:
        by_genome.setdefault(r[idx["genome"]], []).append(r)

    out_rows = []
    for genome in sorted(by_genome):
        contigs = load_fasta(ROOT / "data" / "genomes" / f"{genome}.fna")
        meta = genome_meta[genome]
        ggc = meta["observed_gc"] / 100.0

        for r in by_genome[genome]:
            seqid = r[idx["seqid"]]
            start, end = int(r[idx["start"]]), int(r[idx["end"]])
            strand = r[idx["strand"]]
            product = r[idx["product"]]
            contig = contigs[seqid]

            window = contig[start - 1:end]
            coding = window if strand == "+" else revcomp(window)
            length = len(window)

            up = contig[max(0, start - 1 - FLANK):start - 1]
            down = contig[end:end + FLANK]

            longest_orf, min_stops = six_frame_orf_stats(window)
            gc = gc_fraction(window)
            acgt = sum(window.count(b) for b in "ACGT")

            values = {
                "gc_content": round(gc, 5),
                "gc_skew": round(skew(window, "G", "C"), 5),
                "at_skew": round(skew(window, "A", "T"), 5),
                "gc_minus_genome": round(gc - ggc, 5),
                "entropy_3mer": round(entropy_3mer(window), 5),
                "max_homopolymer": max_homopolymer(window),
                "longest_orf_frac": round(longest_orf / length, 5) if length else 0,
                "min_stops_per_kb": round(min_stops, 5),
                "flank_gc_up": round(gc_fraction(up), 5) if up else 0,
                "flank_gc_down": round(gc_fraction(down), 5) if down else 0,
                "dist_to_contig_end": min(start - 1, len(contig) - end),
                "near_contig_edge": int(min(start - 1, len(contig) - end) < EDGE_THRESHOLD),
                "ambiguous_frac": round((length - acgt) / length, 5) if length else 0,

                "genome_gc": round(ggc, 5),
                "genome_length_bp": meta["length_bp"],
                "genome_n_contigs": meta["contigs"],

                "length_bp": length,
                "strand_plus": int(strand == "+"),
                "gc3": round(gc_fraction(coding[2::3]), 5),
                "start_atg": int(coding[:3] == "ATG"),
                "start_gtg": int(coding[:3] == "GTG"),
                "start_ttg": int(coding[:3] == "TTG"),
                "stop_taa": int(coding[-3:] == "TAA"),
                "stop_tag": int(coding[-3:] == "TAG"),
                "stop_tga": int(coding[-3:] == "TGA"),
                "is_hypothetical": int(
                    ("hypothetical" in product.lower())
                    or ("unknown" in product.lower())
                    or (product.strip() == "")
                ),
                "has_gene_symbol": int(bool(r[idx["gene"]].strip())),
            }
            out_rows.append(
                [r[idx[c]] for c in ID_COLUMNS] + [values[f] for f in FEATURES]
            )
        print(f"  {genome} {len(by_genome[genome]):>6,} intervals", flush=True)

    with gzip.open(OUT_TSV, "wt", newline="") as fh:
        fh.write("\t".join(ID_COLUMNS + FEATURES) + "\n")
        for r in out_rows:
            fh.write("\t".join(str(v) for v in r) + "\n")

    groups = {}
    for name, (group, _desc) in FEATURE_MANIFEST.items():
        groups.setdefault(group, []).append(name)

    label_i = len(ID_COLUMNS) - 1
    labels = [int(r[label_i]) for r in out_rows]
    METRICS.write_text(json.dumps({
        "step": "06_features",
        "n_rows": len(out_rows),
        "n_positive": sum(labels),
        "n_features": len(FEATURES),
        "feature_groups": {g: len(f) for g, f in groups.items()},
        "features_by_group": groups,
        "manifest": {
            name: {"group": grp, "description": desc}
            for name, (grp, desc) in FEATURE_MANIFEST.items()
        },
        "circularity_note": (
            "The 'caller' group is derived from Bakta's own decision or "
            "database lookup, not from sequence. is_hypothetical is close to a "
            "restatement of the label. Step 10 refits without the whole group "
            "and reports both numbers; the sequence-only number is the one "
            "that answers the question."
        ),
        "table": str(OUT_TSV.relative_to(ROOT)),
    }, indent=2) + "\n")

    print(f"\n{len(out_rows):,} rows x {len(FEATURES)} features")
    for g, f in sorted(groups.items()):
        print(f"  {g:9} {len(f):>2} features")
    print(f"wrote {OUT_TSV.relative_to(ROOT)}")
    print(f"wrote {METRICS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

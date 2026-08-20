# Lab notebook

Working notes. Not published. Numbers here are provisional until the matching
file lands in `results/metrics/`.

## 2026-08-20 — panel chosen, annotation started

### Genome panel

25 species, one genome each, complete genomes only, GC 28.45–72.12 %
(observed, from `results/metrics/01_genomes.json`). Total 95.2 Mbp.

Diversity is the point. GC is a candidate feature, so a panel drawn from one
clade would confound "GC" with "this lineage" and the forest would learn the
lineage. Spread across 8 phyla instead.

Two panel decisions worth remembering:

**Mollicutes excluded.** *Mycoplasma* and *Ureaplasma* use genetic code 4, not
11. Including them would mean per-genome translation tables in both
annotators, turning a uniform comparison into a conditional one. Recorded in
`data/species.tsv` rather than quietly dropped.

**Complete genomes only, enforced in code.** NCBI's designated reference for
*Enterococcus faecalis* (GCF_000393015.1) is scaffold-level. Both annotators
truncate an ORF at a contig boundary, so a fragmented assembly carries a
fragmentation-driven disagreement rate that no complete genome in the panel
shares — a per-genome artefact landing in a per-genome grouped split. The
resolver now rejects a non-complete reference and falls back to the best
complete RefSeq assembly, recording which path it took in
`resolved_by`. E. faecalis resolved to GCF_001598635.1 instead.

Accessions are pinned in `data/accessions.tsv` and NCBI is not consulted again
on later runs. NCBI re-designates references over time; a panel that silently
drifts is not a panel.

### Two environment bugs, both silent, both label-affecting

Neither of these announced itself. Both are recorded because both would have
changed the numbers.

**1. Prokka died before annotating anything.** Its `minced` dependency is
Java. An activated conda base env exports `JAVA_HOME=/home/qasim/miniconda3`
and `JAVA_LD_LIBRARY_PATH`; Prokka's bundled JVM honours them, looks for its
own runtime libraries under the wrong prefix, and dies with
`Could not resolve "ZIP_Open"`. Same class of leak as `PERL5LIB`, which
`lib_env.sh` already guarded. Both are now cleared in `run_prokka`.

Loud failure, easy to catch.

**2. Prokka found zero rRNAs in all 25 genomes — and still exited 0.** This
is the dangerous one. Prokka's cached env ships `libgsl.so.27`, but the
`nhmmer` binary in it was linked against `libgsl.so.25`. `barrnap` shells out
to `nhmmer`, so rRNA detection failed at runtime while Prokka reported
`Found 0 rRNAs` and finished successfully.

Why it matters: with no rRNA from Prokka, every rRNA locus Bakta calls looks
like a tool disagreement when it is really a missing shared library. That is
a broken install writing itself into the label.

Fixed by resolving the missing soname from a sibling cached env, appended
*after* Prokka's own lib dir so Prokka's libraries still win every other
lookup. First full Prokka run (87,901 CDS, 0 rRNA) was discarded and re-run.

The general lesson for this project: an annotator exiting 0 is not evidence
it did its job. Feature-type counts per tool are now recorded in
`results/metrics/` for exactly this reason.

### Fairness of the comparison

Neither tool is given `--genus`/`--species`. Both change their output when
hinted, and hinting one but not the other would make the disagreement partly
a function of the hint rather than of the sequence. Defaults on both sides,
kingdom Bacteria, translation table 11.

Bakta gets `--keep-contig-headers`, which is load-bearing rather than
cosmetic: by default Bakta renames sequences to `contig_1, contig_2, ...`
while Prokka keeps the accession. Two interval tables keyed on different
sequence names cannot be overlapped, and the failure mode is silent — every
region would look like a disagreement because no interval would ever share a
seqid. `04_extract_calls.py` now fails loudly if the identifiers stop lining
up.

### Bakta database

db-light, not db-full: ~13 GB free on this machine and db-full does not fit.
Recorded as a result-affecting constraint in
`results/metrics/02_annotation_versions.json`, not buried in an install note.
db-light has fewer reference proteins, so some CDS that db-full would name
are left hypothetical. It affects annotation *content*; it does not change
which regions are called as CDS. The label is about which regions are called,
so the effect on the label should be small — but "should be" is not measured,
and that is a limitation to state plainly rather than assume away.

Zenodo download is slow (~1.34 GB at 50–300 KB/s, one connection dropped at
48 %). Switched to a resumable multi-connection transfer.

### Still open

- Bakta annotation not yet run; waiting on the database.
- Nothing downstream of `04` has been executed on real data yet.

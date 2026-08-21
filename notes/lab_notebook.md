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

## 2026-08-20 (later) — the pipeline ran, and the premise mostly failed

Everything executed end to end. The result is a negative one, and it is more
interesting than the result the experiment was designed to find.

### The two tools barely disagree about where genes are

87,960 Bakta CDS calls. **72** have no overlapping Prokka feature. Positive
rate **0.0008**. Among the 87,888 matched calls, **87,788 share an identical
start *and* stop** — not merely overlapping, byte-identical coordinates.

The reason is visible in the GFF `source` column:

    bakta    87,917 Pyrodigal  +  43 Bakta (sORF module)
    prokka   87,860 Prodigal:v2.6

Both tools delegate CDS calling to the same algorithm. Pyrodigal is a Cython
reimplementation of Prodigal. Asking whether Bakta and Prokka disagree about
gene boundaries is largely asking whether Prodigal agrees with itself, and it
does: 38 divergences in 87,917 calls, 0.04%.

This was knowable before a single genome was downloaded, by reading what each
tool wraps. It was not knowable from their output, which is why it is written
down here rather than quietly fixed.

### What the 72 disagreements actually are

Of the 72 positives, **34 come from Bakta's sORF module** — a separate short-ORF
detection stage that Prokka has no equivalent of. That is 34 of the 43 sORF
calls Bakta made in total: nearly every sORF Bakta finds is a disagreement,
because Prokka structurally cannot find them. Prodigal's default minimum gene
length is 90 bp; the shortest positive here is 45 bp.

Positive lengths: median **118 bp** (~39 aa), 43 of 72 under 300 bp.
Negative lengths: median **816 bp**.

The products name the mechanism outright: Type I toxin-antitoxin system Ibs
family toxins, *trp* and *his* operon leader peptides. Classic short
regulatory peptides, exactly the class a dedicated sORF finder exists to catch.

So the label is not "these tools read the sequence differently". It is
"Bakta runs one extra module". That is a fact about software architecture, not
about DNA.

### The model numbers, and why they should not be quoted

Held-out set: 18,735 rows, **13 positives**.

    no-skill AP           0.00069
    length_bp only        0.542
    forest, all features  0.589
    forest, sequence only 0.470

Two things make these unquotable as point estimates. First, 13 positives:
re-ranking one of them moves average precision by tenths. `lib_model.evaluate`
now attaches a reliability warning to any metrics file computed on fewer than
30 positives, so the caveat travels with the number instead of living only in
this notebook.

Second, and worse: a single feature, `length_bp`, gets AP 0.542 against the
full forest's 0.589. The forest adds almost nothing over "is this call
short?". And `length_bp` sits in the **caller** group — it is the length of the
interval Bakta chose. The model is not reading sequence. It is reading the
footprint of the sORF module, which is precisely the circularity step 10 was
built to detect. The audit fires as designed: sequence-only drops to 0.470 and
keeps most of its apparent skill, because short ORFs also look different in
composition — but on 13 positives that residual is not distinguishable from
noise.

### Honest conclusion

The question as posed — can Bakta/Prokka CDS disagreement be predicted from
sequence — does not have enough disagreement to answer. Not because the
modelling failed, but because the premise does not hold for these two tools at
the interval level.

The comparison that *does* have signal is annotation **content**: both tools
call the same 87.9k regions and then differ in what they name them. That is a
different experiment, on the same annotation output already sitting on disk,
and it needs no new compute.

If the interval question is worth keeping, it needs a gene-caller pair that
does not share an implementation — Prodigal against GeneMarkS-2 or Glimmer,
say. Comparing two wrappers around one caller cannot produce disagreement to
model.

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

## 2026-08-22 — the content experiment. Predictable, but not from sequence.

Part one's negative result stands untouched. Scripts 01–11 and their metrics
files were not edited. This is a second experiment on the same annotation
output, no new genomes and no re-annotation.

Housekeeping first: the tree was already clean at `414f949`, and the 3.8 GB
Bakta database was deleted once it was certain nothing downstream needs it —
the whole content experiment reads GFF and `.faa` files. 3.9 GB free became
7.8 GB. `provenance/fetch_bakta_db.sh` re-downloads the pinned release.

### Two defects in part one's artefacts, found and reported, not fixed

`has_gene_symbol` is **constant 1** across all 87,960 part-one rows.
`04_extract_calls.py` reads the gene symbol as `a.get("gene", a.get("Name",
""))`, and Bakta sets `Name=product` on every CDS, so the fallback filled the
column with product strings. Bakta's real `gene=` is present on 19,216 CDS,
not all of them.

No part-one number moves: a constant feature can never be split on and its
permutation importance is exactly zero. But it was a dead feature shipped
without anyone noticing, which is the same class of failure as Prokka exiting
0 with no rRNAs. Step 14 now hard-fails on any constant feature. That check
fired on its first run, on `bakta_has_pfam`, for a real structural reason:
every one of the 3,428 Bakta CDS carrying a PFAM cross-reference has a
placeholder product, so a PFAM hit and a Bakta name are mutually exclusive
here. Declared as a verified exclusion rather than silently dropped, and the
declaration itself fails if it ever stops being true.

### The db-light confound runs the opposite way to what was assumed

The brief expected db-light to leave Bakta under-named relative to Prokka. It
does relative to db-full — but Bakta on db-light still leaves only 12.5% of
paired CDS as a placeholder against Prokka's 40.2%. The big asymmetric cell is
**Bakta named, Prokka did not** (26,699), against 2,297 the other way.

db-light can only make Bakta name *fewer* regions, so it works against that
asymmetry rather than producing it. A db-full run would push rows from
`both_placeholder` and `prokka_named_only` into `bakta_named_only` and make
the imbalance larger. So `bakta_named_only` is a **lower bound**, and calling
the asymmetric cells a "db-light artefact" would have been wrong in a way that
flattered the constraint. What db-light genuinely confounds is the absolute
size of `both_placeholder` and `prokka_named_only`.

### A normalisation rule that split strings from themselves

First draft of the gene-symbol strip keyed on each record's own `gene=`
attribute. Bakta populates it on 19,216 CDS and Prokka on 48,201, so the strip
fired on one side and not the other, and `'GTPase Era'` was scored as a
disagreement against `'GTPase Era'` 23 times. `'Co-chaperonin GroES'` and
`'Transcription termination factor Rho'` did the same.

Fixed by building the symbol set once per pair from both tools. The invariant
— identical raw strings can never be scored as a disagreement — is now an
assertion in `lib_names.assert_symmetric`, called by step 12, and it refuses
to write the label file rather than warning. 29,808 pairs have identical raw
strings; zero violations.

Worth recording that this bug inflated the disagreement rate in the direction
that made the experiment look more interesting. It would have survived any
amount of staring at summary statistics.

### Two declared deviations from the brief

**The selection rule.** The brief called for "part one's rule verbatim: mean
MCC across five grouped folds, 1-SE tie-break". That is not part one's rule —
`lib_model.sweep_depth` selects by mean average precision with a plain argmax
and no 1-SE step. Implemented what was specified, which also suits a
51.7%-positive problem better than AP does, and recorded the discrepancy in
every sweep file instead of quietly picking one.

**Gradient boosting.** `HistGradientBoostingClassifier` rather than
`GradientBoostingClassifier`: 54.1 s versus 0.7 s per fit here, so one sweep
is 27 minutes against 0.4. Decided before any boosting result existed and
recorded with the measurement.

### The specified baseline was degenerate

"Is either product hypothetical" cannot be computed on the primary analysis
set, because that set is *defined* as both tools having named the region. It
is constant 0 on every row it would be scored on. Reported as uncomputable
rather than swapped for something that produces a number; a substitute
carrying the same intent — did Prokka find a UniProtKB similarity hit —
was added and labelled as a substitute.

### The result

50,210 both-named regions, 51.7% name disagreement under the primary rule
(56.0% strict, 50.6% loose). A real, balanced problem, unlike part one's 13
positives.

    random forest, test MCC   0.303
    gradient boosting         0.288
    decision tree             0.261
    database coverage         0.141   <- baseline
    protein length            0.028
    majority class            0.000

And then the audit:

    full            61 features   CV 0.293   test 0.303
    no caller       26 features   CV 0.232   test 0.259
    sequence only   17 features   CV 0.070   test 0.142

Sequence-only collapses to the database-coverage baseline. Permutation
importance says it without ambiguity: 9 db-derived features are worth 0.165 in
MCC, 35 caller-derived features 0.019, and the 17 sequence and genome features
**−0.001**. Not "small" — negative, i.e. indistinguishable from noise.

So the answer is that Bakta/Prokka naming disagreement *is* predictable, at
MCC ~0.30, and what predicts it is how well each tool's reference database
covers the region. `bakta_n_dbxref`, `prokka_has_ec`,
`prokka_has_protein_motif`. That is a fact about reference databases, not
about DNA.

The honest reading: the second experiment did not fail the way the first did.
The label had plenty of signal and the models found it. But the question on
the tin was "can sequence predict this", and the sequence-only grouped-CV MCC
of 0.070 says no. The audit earned its place for a second time, in a different
way — part one it caught a module footprint, here it caught database coverage.

The `db_derived` flag was the load-bearing design decision of this half. Had
the features carried only part one's caller flag, the `no_caller` arm would
have scored 0.259 and looked like a sequence result.

### OOB against grouped CV

Forest OOB MCC 0.340 against a grouped-CV mean of 0.293. Gap +0.047, and in
the expected direction: OOB bootstraps rows, so rows from the same genome sit
both in-bag and out-of-bag. The gap is a measurement of genome-level leakage,
which is the reason the split is grouped by genome in the first place.

### Loose ends

- 8.4% of the primary set is Bakta's `... domain-containing protein` naming
  style, which disagrees ~100% of the time. Declared before fitting, counts in
  `12_content_cohort.json`. Not excluded, because it is a real Bakta name.
- Six features are near-constant (>99.5% one value) and recorded as such:
  `near_contig_edge`, `ambiguous_frac`, `same_start`, `same_stop`,
  `neighbourhood_identical_calls`, `bakta_has_uniref`. `same_start`/`same_stop`
  being near-constant is part one's finding restated.
- `scripts/verify_readme_numbers.py` checks every quoted number against the
  metrics files. It caught two errors in the first README draft: the
  majority-class accuracy (0.503 is `majority_class_accuracy`, the best
  possible constant; the baseline actually scores 0.497 because its class is
  fixed on the training genomes) and the third-ranked feature. Both were
  exactly the kind of transcription error the standard exists to prevent.

## 2026-08-22 — part two ran, and the premise held this time

The content experiment executed end to end on the annotation already on disk.
No new genomes, no re-annotation. The Bakta database was deleted first (4.0 GB,
root filesystem was at 98 %); `provenance/fetch_bakta_db.sh` re-downloads the
pinned release and nothing downstream needs it.

### The headline is the rate, not the model

50,210 regions where both tools produced a real product name. They disagree on
**25,977 of them, 51.7 %**. That is the finding. It needs no classifier.

The model numbers are real but modest and should not be allowed to displace
it: forest test MCC 0.303, accuracy 0.651. Reported, not led with.

### The objection, pre-empted

Three Bakta naming conventions cannot match Prokka by construction:
`… domain-containing protein`, `Uncharacterized protein <ORF-name>`, and DUF
numbers. Declared before fitting, so this is not a post-hoc carve-out.

Counting them needed care. Matched independently they total 4,650 rows, but
**every one of the 156 DUF names also ends in `domain-containing protein`**, so
the sum double-counts. The distinct union is **4,494**. Publishing 4,650 would
have been wrong by 156.

Those 4,494 disagree 99.9 % of the time. The remaining 45,716 disagree 47.0 %.
Both numbers go in the README. The finding barely moves, which is the point.

### The db-light confound runs the opposite way to what was assumed

Worth recording because the brief for this half assumed the other direction.

db-light has fewer reference proteins, so it can only make Bakta name *fewer*
regions than db-full would. But Bakta on db-light still leaves only 12.5 % of
paired CDS as a placeholder against Prokka's 40.2 %. The large asymmetric cell
is *Bakta named, Prokka did not* (26,699) versus the reverse (2,297).

So db-light works **against** the observed asymmetry rather than producing it.
26,699 is a lower bound, and 51.7 % is plausibly an underestimate: db-full
would add exactly the regions Bakta currently cannot name, which are
thin-evidence cases, and that is where disagreement concentrates.

### The mechanism is two mechanisms, and one of them inverts

The importance ranking makes a single chain look obvious — weak coverage →
Bakta fallback name → disagreement. It is cheap to test, so I tested it
(`22_content_mechanism.json`) rather than writing it up as a mechanism.

On Prokka's side it holds: 54.0 % disagreement without an EC number, 42.1 %
with one.

On Bakta's side it **reverses**. Disagreement rises with cross-reference
count: 44.2 % at the minimum two, 79.8 % at three, 62.8 % at four or more. A
third cross-reference is usually an EC number, a BlastRules hit or a
virulence-factor match, and it comes with a *more specific* Bakta name — which
is then more likely to differ from Prokka's generic one.

`bakta_n_dbxref` is the forest's top feature (+0.087) for that reason, not the
one that looked obvious. Two mechanisms pulling opposite ways, not one.

### The circularity audit fired again, harder

    full             61 features   CV MCC 0.293   test MCC 0.303
    no_caller        26 features   CV MCC 0.232   test MCC 0.259
    sequence_only    17 features   CV MCC 0.070   test MCC 0.142

Summed permutation importance: db-derived **0.165** across 9 features,
caller-derived 0.019 across 35, sequence and genome **−0.001** across 17.

The sequence features contribute nothing. The model reads database-search
output. That is the answer, and it is the answer the `db_derived` flag was
added to be able to give.

Note the sequence-only arm's CV/test gap (0.070 vs 0.142). The grouped-CV
number is the one to trust for generalisation; a five-genome test set moves.

### Things that went wrong, and the checks that caught them

**The normalisation rule was not symmetric.** First draft stripped a trailing
gene symbol using *each record's own* `gene=` attribute. Bakta emits `gene=` on
19,216 of 87,859 CDS and Prokka on 48,201, so the strip fired on one side only
and drove identical raw strings apart — `GTPase Era` vs `GTPase Era` scored as
a disagreement 23 times. The symbol set is now built once per pair from both
tools. `lib_names.assert_symmetric` fails the run if it ever regresses.

**`has_gene_symbol` was dead in part one.** `04_extract_calls.py` reads
`a.get("gene", a.get("Name", ""))` and Bakta sets `Name=product`, so the column
was filled with product strings and the feature was constant 1 across all
87,960 rows. A constant feature cannot be split on, so no part-one number
moves — scripts 01–11 are untouched — but it was dead weight nobody noticed.
Step 14 now fails on any constant feature. That check immediately caught
`bakta_has_pfam`: every Bakta CDS carrying a PFAM cross-reference has a
placeholder product, so it is excluded from the primary set by construction.
Handled as a *declared and verified* exclusion, not a silent drop — the run
also fails if a declared-constant feature turns out to vary.

**A specified baseline was degenerate.** "Is either product hypothetical"
cannot be computed on a cohort defined as both tools having named the region.
Reported as uncomputable rather than quietly swapped; the database-coverage
predictor is labelled as its substitute.

**A figure would have crashed.** `FLAG_COLOUR` had no entry for the `genome`
group. Now raises rather than falling back to a colour that reads as a
different group.

### Two documented deviations from the brief

Recorded in the metrics files rather than resolved silently.

1. The brief called the selection rule "part one's, verbatim". It is not —
   `lib_model.sweep_depth` uses mean average precision with a plain argmax and
   no one-standard-error step. Implemented as *specified* (MCC + 1-SE), which
   also suits a 51.7 %-positive problem better than AP does.
2. Gradient boosting uses `HistGradientBoostingClassifier`. Measured 54.1 s
   versus 0.7 s per fit here — one sweep would have been 27 minutes against
   0.4. Decided before any gradient-boosting result existed.

### Still open

- The interval question still needs a gene-caller pair that does not share an
  implementation. Nothing here changes that.
- `bakta_pfam_style_name` as an explicit flag, and a with/without-families
  model comparison, were considered and not run. The 47.0 % remainder rate
  covers what they would have shown.

# bakta_prokka_disagreement_prediction

Machine Learning for Biology | Chapter 3.

Two experiments on the same 25 annotated genomes. The first asks whether the
two annotators disagree about **where** genes are, and finds that they barely
do. The second asks whether they disagree about **what** those genes are, and
finds that they do, roughly half the time — but that the disagreement is
predictable from database coverage rather than from sequence.

Both results are negative in the sense that neither confirms the hypothesis it
was built to test. Both are reported in full.

## What the labels are, and are not

Neither label is biology. Both are tool output.

**Part one.** A positive case means Bakta made a CDS call where Prokka did
not. It does not mean a gene is present.

**Part two.** A positive case, `name_disagreement`, means the two tools wrote
different product names for the same region. It is not a claim that either
name is correct, and it is not a claim about what the protein does. The
columns are `product_bakta` and `product_prokka`; there is no column anywhere
in this repository for a correct answer, because the experiment does not
observe one.

---

## Part one — where genes are

Can Bakta/Prokka CDS disagreement be predicted from sequence alone?

**The premise does not hold for these two tools.** 87,960 Bakta CDS calls.
**72** have no overlapping Prokka feature, a positive rate of **0.0008**.
Among the 87,888 matched calls, **87,788** share an identical start *and*
stop — not merely overlapping, byte-identical coordinates.

The reason is visible in the GFF `source` column. Bakta delegates CDS calling
to Pyrodigal and Prokka to Prodigal, and Pyrodigal is a Cython
reimplementation of Prodigal. The comparison was largely asking whether
Prodigal agrees with itself. It does.

Of the 72 positives, **34** come from Bakta's sORF module, a short-ORF stage
Prokka has no equivalent of. The label is not "these tools read the sequence
differently"; it is "Bakta runs one extra module".

The held-out set carried **13 positives**, so the model numbers
(average precision 0.589 with all features, 0.470 sequence-only, against a
no-skill floor of 0.00069) are not quotable as point estimates. A reliability
warning is attached to any metrics file computed on fewer than 30 positives.
The circularity audit fired as designed: the strongest single feature was
`length_bp`, which is the length of the interval Bakta chose — the footprint
of the sORF module, not a property of DNA.

If the interval question is worth keeping, it needs a caller pair that does
not share an implementation.

## Part two — what genes are

Both tools call the same ~87.9k regions. Do they agree on what those regions
are? Product naming comes from database search and curation rules the two
tools do not share, so unlike gene-boundary calling it is not one algorithm
answering twice.

No new genomes, no re-annotation. This half runs on the annotation output
already on disk.

### The cohort

87,960 Bakta CDS, of which **87,859** pair with a same-strand Prokka CDS.
101 have no Prokka CDS to compare against and are out of scope here.

| | n | share |
|---|---:|---:|
| both tools named the region — **primary analysis set** | **50,210** | 57.1% |
| Bakta named only | 26,699 | 30.4% |
| both left a placeholder | 8,653 | 9.8% |
| Prokka named only | 2,297 | 2.6% |

**The primary analysis is conditional on both tools naming the region.** That
restriction is imposed by what is measurable, not chosen: you cannot compare
two names when one tool did not produce one.

**On db-light.** Bakta ran against db-light rather than db-full because the
machine had ~13 GB free. db-light has fewer reference proteins, so it can only
make Bakta name *fewer* regions than db-full would. It therefore works
*against* the asymmetry above rather than producing it — the large cell is
"Bakta named, Prokka did not", and a db-full run would make it larger.
`bakta_named_only` is a **lower bound** under db-light. What db-light does
confound is the absolute size of `both_placeholder` and `prokka_named_only`,
which should not be read as properties of Bakta in general.

### The rules, fixed before any model was fitted

**Placeholders** are matched as exact literals after case-folding and
whitespace normalisation, never as substrings — `hypothetical` occurs inside
12 distinct Bakta products of which 11 are real names, and `conserved` inside
37 of which nearly all are. The enumerated list and its counts are in
`13_name_rules.json`.

**Name normalisation** strips trailing bracketed qualifiers, one trailing
gene-symbol token, and leading hedges, then case-folds and reduces punctuation
to spaces. The gene-symbol strip is built once per pair from *both* tools'
`gene=` attributes: keyed per record it fires on one side only and drives
identical raw strings apart, which scored `'GTPase Era'` against
`'GTPase Era'` as a disagreement 23 times before it was caught. The code
asserts that identical raw strings can never be scored as a disagreement.

Two sensitivity checks were named in advance and are reported alongside, never
substituted for the primary rule:

| rule | name_disagreement |
|---|---:|
| strict — case-fold and whitespace only | 28,121 (56.0%) |
| **primary — modelled** | **25,977 (51.7%)** |
| loose — generic tokens dropped, token-set equality | 25,418 (50.6%) |

Normalisation moves the rate by 5.4 points across its full declared range. The
disagreement is substantive, not typographic.

**Gene symbols and EC numbers are separate reported columns**, never part of
the product comparison and never used as features — either would be a second
measurement of the label. Where both tools assign a gene symbol they differ on
41.1% of cases; where both assign an EC number, 16.7%.

### Result: predictable, but not from sequence

61 features, each carrying two flags: `caller_derived` (encodes a gene-caller
decision) and `db_derived` (computed from either tool's database-search
output). Split grouped by genome, five held out, overlap asserted to be zero.
Hyperparameters by mean MCC across five genome-grouped folds with the
one-standard-error tie-break to the simplest model.

Scored on the same 12,283 held-out rows, 6,107 positive:

| | test MCC | accuracy | F1 |
|---|---:|---:|---:|
| majority class | 0.000 | 0.497 | 0.664 |
| protein length, one threshold | 0.028 | 0.511 | 0.604 |
| database coverage | 0.141 | 0.533 | 0.145 |
| decision tree | 0.261 | 0.630 | 0.622 |
| **random forest** | **0.303** | 0.651 | 0.649 |
| gradient boosting | 0.288 | 0.643 | 0.657 |

Gradient boosting does not improve on the forest. The majority-class baseline
takes the highest F1 in the table while its MCC is exactly zero, because it
predicts the positive class on every row; MCC is the headline metric for that
reason. Its accuracy of 0.497 is below the 0.503 a constant predictor tuned on
the test set would reach: the class it predicts is fixed by the training
genomes, and the held-out genomes are fractionally negative-majority.

**The circularity audit is the result.**

| arm | features | grouped CV MCC | test MCC |
|---|---:|---:|---:|
| full | 61 | 0.293 | 0.303 |
| caller-derived removed | 26 | 0.232 | 0.259 |
| **sequence only** | 17 | **0.070** | **0.142** |

Removing the caller features costs little. Removing the database features
collapses the model to the database-coverage baseline it was supposed to beat.
Permutation importance says the same thing directly: summed over the held-out
genomes, the 9 `db_derived` features are worth **0.165** in MCC, the 35
`caller_derived` features **0.019**, and the 17 sequence and genome features
**−0.001**.

So naming disagreement between Bakta and Prokka *is* predictable — from how
well each tool's reference database covers the region, not from the DNA. The
top features are `bakta_n_dbxref`, `prokka_has_ec` and
`prokka_has_protein_motif`, with `prokka_has_similarity_hit` fourth. That is a fact about reference databases, not
about sequence, and the sequence-only arm's grouped-CV MCC of 0.070 is the
honest answer to the question actually asked.

**OOB against grouped CV.** The forest's OOB MCC is 0.340 against a grouped-CV
mean of 0.293, a gap of **+0.047**. OOB bootstraps rows, not genomes, so rows
from the same genome land both in-bag and out-of-bag. The gap measures
genome-level leakage; it is not a discrepancy to be explained away.

---

## Provenance

Annotation depends on database release as much as on tool version. Both are
recorded in `provenance/`. The same genome annotated against a newer Bakta
database will produce different calls. The database used here was db-light
v6.0 (2025-02-24, DOI 10.5281/zenodo.14916843); it has been deleted from this
machine to free disk and `provenance/fetch_bakta_db.sh` re-downloads it.

## Structure

    data/          genomes and annotation output (gitignored)
    scripts/       numbered pipeline, run in order
    results/       metrics as numbered JSON, one file per step
    figures/       PNG and PDF
    provenance/    tool versions, database versions, checksums
    notes/         working notes, not published

Steps 01–11 are part one, 12–21 part two. `results/metrics/20_content_summary.json`
carries every headline number for part two in one file.

## Reproducing

    conda env create -f environment.yml
    conda activate bpdp
    bash scripts/00_run_all.sh          # everything, including annotation
    bash scripts/run_content_only.sh    # part two only, on annotation already on disk

Seed 42. Every number quoted above maps to a file in `results/metrics/`.

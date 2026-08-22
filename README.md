# bakta_prokka_disagreement_prediction

Machine Learning for Biology | Chapter 3.

Two experiments on one panel of 25 complete bacterial genomes (95.2 Mbp,
8 phyla, GC 28.45–72.12 %), annotated once with Bakta and once with Prokka.

---

## The headline

**Bakta and Prokka give the same coding region a different product name about
half the time.**

Of 87,859 CDS regions both tools called, 50,210 were given a non-placeholder
product name by both. The two tools disagree about that name on **25,977 of
them — 51.7 %**.

That number needs no model to matter. Product names are what pangenome
analyses, functional enrichment tests and AMR reports are built on, and two
standard tools reading identical DNA return different names for half of it.

The obvious objection is that Bakta has naming conventions which cannot match
Prokka by construction — `… domain-containing protein`, `Uncharacterized
protein <ORF-name>`, and names carrying a DUF number. Those were declared
before any model was fitted. They are **4,494 regions (9.0 %) and they
disagree 99.9 % of the time**. Excluding them:

| set | n | name disagreement |
|---|---:|---:|
| all both-named regions | 50,210 | **51.7 %** |
| excluding three Bakta fallback-naming conventions | 45,716 | **47.0 %** |
| those three conventions alone | 4,494 | 99.9 % |

Both numbers are published together. The finding survives the objection.

Under the two declared sensitivity checks the rate moves from 56.0 % (strict:
case-folding and whitespace only) to 50.6 % (loose: generic tokens dropped,
compared as unordered token sets). The disagreement is substantive, not
typographic.

## The question, and the answer

Given a region both tools called, can sequence predict whether they will name
it differently?

**Largely no.** A random forest reaches MCC 0.303 on five held-out genomes
against a majority-class floor of 0.000, but the audit shows what it is
reading. With every database-derived and caller-derived feature removed, the
sequence-only arm falls to grouped-CV MCC 0.070. Summed permutation
importance across 61 features: **database-derived 0.165 across 9 features,
caller-derived 0.019 across 35, sequence and genome features −0.001 across
17.**

The model is reading database-search output, not DNA. That is reported as the
result rather than presented as a sequence model.

## What the label is, and is not

The label is **disagreement between two software products about a name**. It
is not a claim that either name is correct, and it is not a claim about what
any protein does. Columns are named `name_disagreement`, `product_bakta`,
`product_prokka` for that reason. Establishing which name is right would
require evidence this repository does not contain.

## The db-light constraint, and which way it pushes

Bakta ran against **db-light**, not db-full: this machine had ~13 GB free and
db-full does not fit. That is recorded as a result-affecting constraint in
`results/metrics/02_annotation_versions.json`.

The primary analysis is **conditional on both tools having named the region**.
That restriction is imposed by what is measurable — two names cannot be
compared when one tool produced none — and not chosen for convenience.

Two consequences, both stated rather than assumed away:

**The naming asymmetry is not a db-light artefact.** The large asymmetric cell
is *Bakta named, Prokka did not* (26,699), not the reverse (2,297). Bakta on
db-light leaves 12.5 % of paired CDS as a placeholder against Prokka's 40.2 %.
db-light can only make Bakta name *fewer* regions than db-full would, so it
works against this asymmetry rather than producing it. 26,699 is a lower
bound.

**51.7 % is plausibly an underestimate.** Under db-full the primary set would
grow, and the regions added would be exactly those Bakta currently cannot
name — thin-evidence cases, which is where disagreement concentrates. The
composition of the primary set is db-light-dependent even though its
definition dodges the asymmetry.

## The mechanism is two mechanisms, pulling opposite ways

The natural reading of the importance ranking is a single causal chain: weak
reference coverage → Bakta falls back on a structural name → Prokka names it
differently → disagreement. Tested rather than asserted
(`22_content_mechanism.json`), that chain holds in one place and **reverses in
another**.

On Prokka's side it holds: within the 45,716 remainder regions, disagreement
is 54.0 % without an EC number against 42.1 % with one.

On Bakta's side it inverts. Disagreement *rises* with Bakta's cross-reference
count — 44.2 % at the minimum two, **79.8 % at three**, 62.8 % at four or
more. A third cross-reference is typically an EC number, a BlastRules hit or a
virulence-factor match, and it comes with a *more specific* Bakta name, which
is then more likely to differ from Prokka's more generic one.

`bakta_n_dbxref` is the forest's strongest single feature (permutation
importance +0.087) because extra Bakta evidence predicts a more specific name
that Prokka does not match — not because thin evidence predicts a fallback
name.

## Part one: the interval experiment, and why it failed

The original question was whether *interval* disagreement — one tool calling a
CDS where the other leaves the region empty — is predictable from sequence.
It is not answerable with this tool pair.

87,960 Bakta CDS calls. **72** have no overlapping Prokka feature: a positive
rate of 0.0008. Among 87,888 matched calls, **87,788 share byte-identical
start and stop coordinates**.

The reason is architectural. Both tools delegate CDS calling to the same
algorithm — Pyrodigal is a Cython reimplementation of Prodigal — so the
comparison largely asked whether Prodigal agrees with itself. It does. Of the
72 positives, 34 come from Bakta's sORF module, which Prokka has no equivalent
of.

That negative result stands and is published in full as steps 01–11. It is
also what motivated part two: the tools agree on *where* genes are and
disagree on *what they are*.

## Decisions fixed before any result was seen

Written into the code and this file before the first model was fitted, and not
revised afterwards.

- **Placeholder strings** are matched as an enumerated list of exact literals
  after case-folding, never by substring. `hypothetical` occurs inside 12
  distinct Bakta products of which 11 are real names. The full list and its
  counts are in `13_name_rules.json`.
- **Name normalisation**: NFKC → strip trailing bracketed qualifiers → strip
  one trailing gene-symbol token → strip leading hedges → case-fold →
  punctuation to whitespace. The gene-symbol strip is **pair-symmetric**; a
  per-record version scored `GTPase Era` against `GTPase Era` as a
  disagreement 23 times. `lib_names.assert_symmetric` fails the run if
  identical raw strings are ever scored as different.
- **Sensitivity checks** (strict, loose) were named in advance, as checks and
  not as alternatives to switch to.
- **EC numbers and gene symbols** are separate reported columns, never
  features — either would be a second measurement of the label. Gene symbols
  differ on 5,656 of 13,760 comparable regions (41.1 %); EC numbers on 125 of
  750 (16.7 %).
- **The split** is grouped by genome, five held out, asserted disjoint at run
  time.
- **Hyperparameters**: mean MCC across five genome-grouped folds, then the
  one-standard-error rule to the simplest model. Grids declared before the
  first fit.

## Results at a glance

Test set: 12,283 regions from five held-out genomes, 6,107 positive.

| | test MCC | accuracy |
|---|---:|---:|
| majority class | 0.000 | 0.497 |
| protein length, one threshold | 0.028 | 0.511 |
| database coverage (substitute baseline) | 0.141 | 0.533 |
| decision tree | 0.261 | 0.630 |
| **random forest** | **0.303** | **0.651** |
| gradient boosting | 0.288 | 0.643 |

Circularity audit, same held-out genomes:

| arm | features | grouped CV MCC | test MCC |
|---|---:|---:|---:|
| full | 61 | 0.293 | 0.303 |
| caller-derived removed | 26 | 0.232 | 0.259 |
| caller- and db-derived removed | 17 | **0.070** | 0.142 |

Forest OOB MCC is 0.340 against a grouped-CV mean of 0.293. That gap of
**+0.047** is a measurement of genome-level leakage, not a discrepancy to
explain away: OOB bootstraps rows, so a held-out row is scored by trees that
saw its neighbours from the same genome and the same annotation run.

The majority-class baseline scores 0.497 accuracy rather than 0.503 because
the majority class of the training genomes is *positive* (52.4 %) while the
held-out genomes are 49.7 % positive — always predicting the training majority
lands just below a coin flip. That is the floor, and it is why MCC rather than
accuracy is the headline metric throughout.

One specified baseline, *is either product hypothetical*, is **degenerate by
construction** — the primary set is defined as both tools having named the
region, so it is constant. It is reported as uncomputable rather than quietly
replaced; the database-coverage row above is labelled as its substitute.

## Provenance

Annotation depends on database release as much as on tool version. Both are
recorded in `provenance/`. Bakta db-light v6.0 (2025-02-24, DOI
`10.5281/zenodo.14916843`). The same genome annotated against a newer database
will produce different calls.

The database itself is not in this repository — it is 4.0 GB and
`provenance/fetch_bakta_db.sh` re-downloads the pinned release.

## Structure

    data/          genomes and annotation output (gitignored)
    scripts/       numbered pipeline
    results/       metrics as numbered JSON, one file per step
    figures/       PNG and PDF
    provenance/    tool versions, database versions, checksums
    notes/         working notes

Steps 01–11 are the interval experiment; 12–22 are the content experiment.

## Reproducing

    conda env create -f environment.yml
    conda activate bpdp
    bash scripts/00_run_all.sh          # everything, from an empty checkout

    bash scripts/run_content_only.sh    # part two only, on annotation already on disk

Part two runs in dependency order rather than numeric order: step 22 computes
numbers step 20 quotes.

Every number in this file maps to a file in `results/metrics/`. Nothing here
was transcribed from a terminal.

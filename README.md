# bakta_prokka_disagreement_prediction

Machine Learning for Biology, Chapter 3.

## The question

Bakta and Prokka annotate the same genome and disagree. When one calls a
coding sequence and the other leaves the region empty, can that
disagreement be predicted from the sequence alone?

## What the label is, and is not

The label is tool output, not biology. A positive case means Bakta made a
call where Prokka did not. It does not mean a gene is present. Any claim
about real genes requires evidence this repository does not contain.

## Provenance

Annotation depends on database release as much as on tool version. Both
are recorded in `provenance/`. The same genome annotated against a newer
Bakta database will produce different calls.

## Structure

    data/          genomes and annotation output (gitignored)
    scripts/       numbered pipeline, run in order
    results/       metrics as numbered JSON, one file per step
    figures/       PNG and PDF
    provenance/    tool versions, database versions, checksums
    notes/         working notes, not published

## Reproducing

    conda env create -f environment.yml
    conda activate bpdp
    bash scripts/00_run_all.sh

Every number quoted in any post maps to a file in `results/metrics/`.

# Coming Soon

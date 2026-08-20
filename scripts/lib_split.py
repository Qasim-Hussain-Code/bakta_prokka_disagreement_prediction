#!/usr/bin/env python3
"""The train/test split, defined once and imported by 07, 08 and 10.

Splitting is by genome, never by interval. Two CDS calls from the same genome
share its GC, its codon usage, its assembly quirks and its annotation run; a
random interval split would put near-copies of the same thing on both sides
and report a score that says nothing about a genome the model has not seen.

Which genomes land in the test set is decided by GC rank, not at random. With
25 genomes and 5 held out, a random draw can easily hand you a test set that
is entirely high-GC, and the resulting number would be a statement about
high-GC genomes rather than about the panel. Ranking by GC and dealing the
genomes round-robin into 5 folds guarantees each fold spans the range.

Fold 0 is the test set and stays closed until the depth sweep is finished.
Folds 1-4 are the cross-validation folds used for model selection.

Deterministic: no seed, no shuffle. The same panel always yields the same
split, so numbers from separate runs are comparable.
"""

N_FOLDS = 5
TEST_FOLD = 0


def assign_folds(genome_gc):
    """{genome: gc} -> {genome: fold}, dealt round-robin by GC rank."""
    ordered = sorted(genome_gc, key=lambda g: (genome_gc[g], g))
    return {g: i % N_FOLDS for i, g in enumerate(ordered)}


def split_genomes(genome_gc):
    """-> (train_genomes, test_genomes, folds), all sorted/deterministic."""
    folds = assign_folds(genome_gc)
    test = sorted(g for g, f in folds.items() if f == TEST_FOLD)
    train = sorted(g for g, f in folds.items() if f != TEST_FOLD)
    return train, test, folds


def describe(genome_gc):
    """Human-readable split summary for the metrics files."""
    train, test, folds = split_genomes(genome_gc)
    return {
        "strategy": "grouped by genome, folds dealt round-robin by GC rank",
        "n_folds": N_FOLDS,
        "test_fold": TEST_FOLD,
        "n_train_genomes": len(train),
        "n_test_genomes": len(test),
        "test_genomes": test,
        "test_gc_range": [
            round(min(genome_gc[g] for g in test), 4),
            round(max(genome_gc[g] for g in test), 4),
        ] if test else [],
        "train_gc_range": [
            round(min(genome_gc[g] for g in train), 4),
            round(max(genome_gc[g] for g in train), 4),
        ] if train else [],
        "folds": {g: folds[g] for g in sorted(folds)},
    }

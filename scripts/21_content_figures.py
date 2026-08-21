#!/usr/bin/env python3
"""Figures for the content experiment. PNG and PDF, numbering continues at 06.

Reads:  results/metrics/12..20_*.json
Writes: figures/06..10_*.png + .pdf

Bars are annotated with the counts behind them wherever a rate is drawn. The
naming label is close to balanced, so unlike part one there is no need to
protect against an imbalanced-metric illusion -- but a rate with no count next
to it is still a number nobody can check.

Every value plotted comes from a metrics file. Nothing is recomputed here.
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
METRICS = ROOT / "results" / "metrics"
FIGURES = ROOT / "figures"

plt.rcParams.update({
    "figure.dpi": 130,
    "savefig.bbox": "tight",
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
})

NAVY = "#17456B"
OCEAN = "#2E6DA4"
RUST = "#C44E52"
SAND = "#DDA15E"
GREY = "#8C8C8C"

# Every group the feature manifest can emit. A missing key here is a KeyError
# at plot time, not a silently odd colour.
FLAG_COLOUR = {"sequence": OCEAN, "genome": NAVY, "caller": SAND, "db": RUST}


def save(fig, stem):
    FIGURES.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"{stem}.{ext}")
    plt.close(fig)
    print(f"  figures/{stem}.png + .pdf")


# A label sitting on its own reference line reads as a strike-through. Every
# reference-line label goes through this: opaque background, clear of the line.
LABEL_BOX = {"facecolor": "white", "edgecolor": "none", "alpha": 0.85,
             "pad": 1.5}


def load(name):
    path = METRICS / name
    if not path.exists():
        raise SystemExit(f"{path} absent")
    return json.loads(path.read_text())


def fig_cohort(cohort, genomes):
    ct = cohort["cross_tabulation"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 5.2),
                                   gridspec_kw={"width_ratios": [1, 1.25]})

    cells = [("both named", ct["both_named"], OCEAN),
             ("Bakta named only", ct["bakta_named_only"], SAND),
             ("both placeholder", ct["both_placeholder"], GREY),
             ("Prokka named only", ct["prokka_named_only"], RUST)]
    left = 0
    total = ct["total"]
    for label, n, colour in cells:
        ax1.barh([0], [n], left=left, color=colour, edgecolor="white")
        if n / total > 0.03:
            ax1.text(left + n / 2, 0, f"{n:,}\n{n/total:.1%}", ha="center",
                     va="center", fontsize=8,
                     color="white" if colour != GREY else "white")
        left += n
    ax1.set_yticks([])
    ax1.set_xlim(0, total)
    ax1.set_xlabel(f"paired CDS regions (n = {total:,})")
    ax1.set_title("Who named what", loc="left", fontsize=10)
    ax1.grid(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c in cells]
    ax1.legend(handles, [f"{l} ({n:,})" for l, n, _ in cells],
               loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2,
               frameon=False, fontsize=8)

    rows = [r for r in cohort["per_genome"] if r["name_disagreement_rate"]]
    gc = {g["accession"]: g["observed_gc"] for g in genomes["genomes"]}
    rows.sort(key=lambda r: r["name_disagreement_rate"])
    y = np.arange(len(rows))
    ax2.barh(y, [r["name_disagreement_rate"] * 100 for r in rows], color=NAVY)
    ax2.set_yticks(y)
    ax2.set_yticklabels([f"{r['species']}" for r in rows], fontsize=6.5)
    ax2.set_xlabel("regions the two tools name differently (%)")
    ax2.set_title("Name disagreement per genome, both-named regions only",
                  loc="left", fontsize=10)
    overall = cohort["primary_analysis_set"]["name_disagreement"]["primary"]["rate"]
    ax2.axvline(overall * 100, color=RUST, lw=1.2, ls="--")
    # Above the axes, not over a bar: drawn inside the plot it read as a
    # strike-through on whichever species happened to sit at that height.
    ax2.text(overall * 100 + 0.6, 0.02, f"panel {overall:.1%}",
             transform=ax2.get_xaxis_transform(), ha="left", va="bottom",
             fontsize=7.5, color=RUST, bbox=LABEL_BOX)
    for i, r in enumerate(rows):
        ax2.text(1, i, f"n={r['primary_set']:,}", va="center", fontsize=6,
                 color="white")
    fig.suptitle("Both tools call the same regions; they disagree about what "
                 "they are", fontsize=11, x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save(fig, "06_name_disagreement_cohort")


def fig_models(summary):
    base = summary["baselines_test_set"]
    models = summary["models_test_set"]
    labels, mccs, f1s, colours = [], [], [], []
    for n, b in base.items():
        if not b.get("computed"):
            continue
        labels.append(n.replace("_", " "))
        mccs.append(b["test_mcc"])
        f1s.append(b["test_f1"])
        colours.append(GREY)
    for n, m in models.items():
        labels.append(n.replace("_", " "))
        mccs.append(m["test_mcc"])
        f1s.append(m["test_f1"])
        colours.append(NAVY)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6), sharey=True)
    y = np.arange(len(labels))
    ax1.barh(y, mccs, color=colours)
    ax1.set_yticks(y)
    ax1.set_yticklabels(labels, fontsize=8)
    ax1.set_xlabel("Matthews correlation coefficient")
    ax1.set_title("Test-set MCC", loc="left", fontsize=10)
    ax1.axvline(0, color="black", lw=0.8)
    for i, v in enumerate(mccs):
        ax1.text(v + 0.005, i, f"{v:.3f}", va="center", fontsize=7.5)

    ax2.barh(y, f1s, color=colours)
    ax2.set_xlabel("F1")
    ax2.set_title("Test-set F1", loc="left", fontsize=10)
    for i, v in enumerate(f1s):
        ax2.text(v + 0.005, i, f"{v:.3f}", va="center", fontsize=7.5)

    n_rows = summary["split"]["n_test_rows"]
    n_pos = summary["split"]["n_test_positive"]
    fig.suptitle(
        f"Baselines (grey) and models (navy), same {n_rows:,} held-out rows, "
        f"{n_pos:,} positive", fontsize=10.5, x=0.02, ha="left")
    # The majority-class baseline predicts the positive class on every row, so
    # it takes the highest F1 in the panel while its MCC is exactly zero. Read
    # the left panel; the right one is here to show why.
    fig.text(0.02, 0.90,
             "Majority class always predicts the positive class: highest F1 in "
             "the panel, MCC exactly 0. MCC is the headline metric for that "
             "reason.", fontsize=8, ha="left", color=GREY)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    save(fig, "07_content_model_comparison")


def fig_circularity(circ):
    arms = circ["side_by_side"]
    order = ["full", "no_caller", "sequence_only"]
    x = np.arange(len(order))
    cv = [arms[a]["grouped_cv_mean_mcc"] for a in order]
    te = [arms[a]["test_mcc"] for a in order]
    nf = [arms[a]["n_features"] for a in order]

    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    w = 0.36
    ax.bar(x - w / 2, cv, w, label="grouped CV (train genomes)", color=OCEAN)
    ax.bar(x + w / 2, te, w, label="held-out genomes", color=NAVY)
    for i, (a, b) in enumerate(zip(cv, te)):
        ax.text(i - w / 2, a + 0.004, f"{a:.3f}", ha="center", fontsize=7.5)
        ax.text(i + w / 2, b + 0.004, f"{b:.3f}", ha="center", fontsize=7.5)
    ref = circ["reference_points_same_test_rows"]
    ax.axhline(ref["database_coverage_substitute_mcc"], color=RUST, ls="--", lw=1.1)
    # Anchored in axes fraction. Anchored at x=len(order)-0.5 it fell outside
    # the x limits and was clipped away entirely.
    ax.annotate(f"database-coverage baseline "
                f"({ref['database_coverage_substitute_mcc']:.3f})",
                xy=(0.012, ref["database_coverage_substitute_mcc"]),
                xycoords=ax.get_yaxis_transform(), xytext=(0, 4),
                textcoords="offset points", ha="left", va="bottom",
                fontsize=7.5, color=RUST, bbox=LABEL_BOX)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{a.replace('_', ' ')}\n{n} features"
                        for a, n in zip(order, nf)])
    ax.set_ylabel("Matthews correlation coefficient")
    ax.set_title("Circularity audit: what survives when the flagged features go",
                 loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    save(fig, "08_content_circularity")


def fig_importance(imp, top=18):
    ranked = imp["ranked_by_permutation"][:top][::-1]
    y = np.arange(len(ranked))
    vals = [r["permutation_importance_mean"] for r in ranked]
    errs = [r["permutation_importance_std"] for r in ranked]
    unknown = sorted({r["group"] for r in ranked} - set(FLAG_COLOUR))
    if unknown:
        raise SystemExit(
            f"FATAL: no colour defined for feature group(s) {unknown}. Add them "
            "to FLAG_COLOUR rather than letting the figure fall back to a "
            "default that reads as a different group.")
    colours = [FLAG_COLOUR[r["group"]] for r in ranked]

    fig, ax = plt.subplots(figsize=(8.4, 6.4))
    ax.barh(y, vals, xerr=errs, color=colours, error_kw={"lw": 0.8, "ecolor": "#444"})
    ax.set_yticks(y)
    ax.set_yticklabels([r["feature"] for r in ranked], fontsize=7.5)
    ax.set_xlabel("drop in test-set MCC when the column is shuffled")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_title(f"Permutation importance, {imp['model'].replace('_', ' ')}, "
                 "measured on held-out genomes", loc="left", fontsize=10)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in FLAG_COLOUR.values()]
    ax.legend(handles, [f"{k} features" for k in FLAG_COLOUR],
              frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    save(fig, "09_content_importance")


def fig_leakage(final, summary, cohort):
    species = {r["genome"]: r["species"] for r in cohort["per_genome"]}
    per_genome = final["per_genome_test"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4),
                                   gridspec_kw={"width_ratios": [1.15, 1]})

    rows = sorted(per_genome, key=lambda r: (r["mcc"] if r["mcc"] is not None else -1))
    y = np.arange(len(rows))
    ax1.barh(y, [r["mcc"] or 0 for r in rows], color=NAVY)
    ax1.set_yticks(y)
    ax1.set_yticklabels([species.get(r["genome"], r["genome"]) for r in rows],
                        fontsize=7.5)
    ax1.set_xlabel("MCC")
    ax1.axvline(final["test"]["mcc"], color=RUST, ls="--", lw=1.2)
    ax1.text(final["test"]["mcc"] - 0.004, 0.02,
             f"pooled {final['test']['mcc']:.3f}",
             transform=ax1.get_xaxis_transform(), ha="right", va="bottom",
             fontsize=7.5, color=RUST, bbox=LABEL_BOX)
    for i, r in enumerate(rows):
        ax1.text(0.005, i, f"n={r['n_rows']:,}", va="center", fontsize=6.5,
                 color="white")
    ax1.set_title("Per held-out genome", loc="left", fontsize=10)

    labels = ["OOB\n(bootstraps rows)", "grouped CV\n(holds out genomes)",
              "held-out\ngenomes"]
    vals = [final["oob"]["oob_mcc"], final["grouped_cv_mean_mcc"],
            final["test"]["mcc"]]
    ax2.bar(np.arange(3), vals, color=[SAND, OCEAN, NAVY], width=0.6)
    for i, v in enumerate(vals):
        ax2.text(i, v + 0.006, f"{v:.3f}", ha="center", fontsize=8)
    ax2.set_xticks(np.arange(3))
    ax2.set_xticklabels(labels, fontsize=8)
    ax2.set_ylabel("MCC")
    gap = final["oob_vs_grouped_cv"]["gap"]
    ax2.set_title(f"Genome-level leakage: OOB exceeds grouped CV by {gap:+.3f}",
                  loc="left", fontsize=10)
    fig.tight_layout()
    save(fig, "10_content_oob_vs_grouped_cv")


def main():
    cohort = load("12_content_cohort.json")
    genomes = load("01_genomes.json")
    summary = load("20_content_summary.json")
    circ = load("19_content_circularity.json")
    imp = load("18_content_importance_random_forest.json")
    final = load("17_content_final_random_forest.json")

    FIGURES.mkdir(exist_ok=True)
    fig_cohort(cohort, genomes)
    fig_models(summary)
    fig_circularity(circ)
    fig_importance(imp)
    fig_leakage(final, summary, cohort)


if __name__ == "__main__":
    main()

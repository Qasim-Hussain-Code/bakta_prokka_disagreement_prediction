#!/usr/bin/env python3
"""All figures. PNG and PDF.

Reads:  results/metrics/*.json
        data/interim/preds_forest.npz, preds_sequence_only.npz
Writes: figures/*.png, figures/*.pdf

Precision-recall, not ROC, for the headline curve. The label is heavily
imbalanced, and ROC curves look reassuring on imbalanced problems because the
large negative class makes the false-positive rate move slowly. The no-skill
line is drawn on every PR panel so a curve can be read against the floor it
has to clear rather than against the top of the axes.
"""

import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import precision_recall_curve

ROOT = Path(__file__).resolve().parent.parent
METRICS = ROOT / "results" / "metrics"
INTERIM = ROOT / "data" / "interim"
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


def save(fig, stem):
    FIGURES.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"{stem}.{ext}")
    plt.close(fig)
    print(f"  figures/{stem}.png + .pdf")


def load_json(name):
    path = METRICS / name
    return json.loads(path.read_text()) if path.exists() else None


def fig_disagreement_by_genome(dis):
    """Positive rate per genome against genome GC."""
    rows = sorted(dis["per_genome"], key=lambda r: r["unique_to_bakta_rate"])
    genomes = load_json("01_genomes.json")
    gc = {g["accession"]: g["observed_gc"] for g in genomes["genomes"]}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 6),
                                   gridspec_kw={"width_ratios": [1.35, 1]})

    labels = [f"{r['species']}" for r in rows]
    vals = [r["unique_to_bakta_rate"] * 100 for r in rows]
    ax1.barh(range(len(rows)), vals, color="#4C72B0")
    ax1.set_yticks(range(len(rows)))
    ax1.set_yticklabels(labels, fontsize=7)
    ax1.set_xlabel("Bakta CDS with no overlapping Prokka feature (%)")
    ax1.set_title("Disagreement rate per genome", loc="left", fontsize=10)

    x = [gc[r["genome"]] for r in rows]
    ax2.scatter(x, vals, s=28, color="#C44E52", zorder=3)
    ax2.set_xlabel("genome GC (%)")
    ax2.set_ylabel("disagreement rate (%)")
    ax2.set_title("Disagreement vs genome GC", loc="left", fontsize=10)
    if len(x) > 2:
        r = np.corrcoef(x, vals)[0, 1]
        ax2.annotate(f"Pearson r = {r:.2f}", xy=(0.05, 0.93),
                     xycoords="axes fraction", fontsize=8)

    fig.suptitle("Label is tool output: Bakta called, Prokka did not. "
                 "It does not mean a gene is present.",
                 fontsize=8, y=1.0, color="#555555")
    save(fig, "01_disagreement_by_genome")


def fig_pr_curves():
    """Held-out precision-recall for both models against the no-skill floor."""
    panels = []
    for stem, title in (("preds_forest", "All features"),
                        ("preds_sequence_only", "Sequence only")):
        path = INTERIM / f"{stem}.npz"
        if path.exists():
            d = np.load(path, allow_pickle=True)
            panels.append((title, d["y"], d["proba"]))
    if not panels:
        return

    fig, ax = plt.subplots(figsize=(6.5, 5))
    colors = ["#4C72B0", "#DD8452"]
    floor = None
    for (title, y, proba), color in zip(panels, colors):
        precision, recall, _ = precision_recall_curve(y, proba)
        ap = np.trapezoid(precision[::-1], recall[::-1]) if hasattr(np, "trapezoid") \
            else np.trapz(precision[::-1], recall[::-1])
        floor = float(np.mean(y))
        ax.plot(recall, precision, color=color, lw=1.8,
                label=f"{title} (AP≈{abs(ap):.3f})")
    if floor is not None:
        ax.axhline(floor, ls="--", lw=1.2, color="#888888",
                   label=f"no skill ({floor:.3f})")
    ax.set_xlabel("recall")
    ax.set_ylabel("precision")
    ax.set_title("Held-out genomes: precision-recall", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    save(fig, "02_precision_recall")


def fig_importance(imp):
    rows = imp["ranked_by_permutation"][:15][::-1]
    palette = {"sequence": "#4C72B0", "genome": "#55A868", "caller": "#C44E52"}
    fig, ax = plt.subplots(figsize=(7.5, 6))
    ax.barh(range(len(rows)),
            [r["permutation_importance_mean"] for r in rows],
            xerr=[r["permutation_importance_std"] for r in rows],
            color=[palette.get(r["group"], "#888888") for r in rows],
            error_kw={"lw": 0.8, "ecolor": "#444444"})
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r["feature"] for r in rows], fontsize=8)
    ax.set_xlabel("permutation importance (drop in average precision)")
    ax.set_title("What the forest used, measured on held-out genomes",
                 loc="left", fontsize=10)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in palette.values()]
    ax.legend(handles, palette.keys(), frameon=False, fontsize=8,
              title="feature group", title_fontsize=8, loc="lower right")
    save(fig, "03_feature_importance")


def fig_audit(audit):
    comp = audit.get("comparison")
    if not comp:
        return
    fig, ax = plt.subplots(figsize=(6, 4.5))
    names = ["all features", "sequence only"]
    vals = [comp["average_precision_all_features"],
            comp["average_precision_sequence_only"]]
    floor = comp["no_skill_average_precision"]
    ax.bar(names, vals, color=["#4C72B0", "#DD8452"], width=0.55)
    ax.axhline(floor, ls="--", lw=1.2, color="#888888")
    ax.annotate(f"no skill ({floor:.3f})", xy=(1.42, floor), fontsize=8,
                color="#555555", va="bottom", ha="right")
    for i, v in enumerate(vals):
        ax.annotate(f"{v:.3f}", xy=(i, v), ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("average precision (held-out genomes)")
    ax.set_title("Circularity audit: what survives without\n"
                 "gene-caller-derived features", loc="left", fontsize=10)
    save(fig, "04_circularity_audit")


def fig_length(dis):
    path = INTERIM / "features.tsv.gz"
    if not path.exists():
        return
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        pos, neg = [], []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            (pos if f[idx["label"]] == "1" else neg).append(float(f[idx["length_bp"]]))
    if not pos or not neg:
        return
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    bins = np.logspace(np.log10(30), np.log10(max(max(pos), max(neg))), 60)
    ax.hist(neg, bins=bins, color="#4C72B0", alpha=0.65, density=True,
            label=f"matched by Prokka (n={len(neg):,})")
    ax.hist(pos, bins=bins, color="#C44E52", alpha=0.65, density=True,
            label=f"Prokka silent (n={len(pos):,})")
    ax.set_xscale("log")
    ax.set_xlabel("called CDS length (bp, log scale)")
    ax.set_ylabel("density")
    ax.set_title("Length of Bakta calls, by whether Prokka agreed",
                 loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    save(fig, "05_length_by_label")


def main():
    FIGURES.mkdir(exist_ok=True)
    print("writing figures")
    dis = load_json("05_disagreement.json")
    if dis:
        fig_disagreement_by_genome(dis)
        fig_length(dis)
    fig_pr_curves()
    imp = load_json("09_feature_importance.json")
    if imp:
        fig_importance(imp)
    audit = load_json("10_circularity_audit.json")
    if audit:
        fig_audit(audit)
    print("done")


if __name__ == "__main__":
    main()

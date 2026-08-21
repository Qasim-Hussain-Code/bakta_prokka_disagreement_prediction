#!/usr/bin/env python3
"""Placeholder detection and product-name normalisation. Imported by 12-20.

This module is the whole content experiment's definition of its label, so it
lives in one file and every script that needs a rule imports it rather than
restating it. The rules were fixed and approved before any model was fitted;
they are not to be adjusted after seeing a score.

Two separate questions, deliberately kept apart:

  is_placeholder(product)
      Does this tool's product string carry any claim at all? Exact literal
      match against an enumerated list, after case-folding and whitespace
      normalisation. Never a substring test -- 'hypothetical' occurs inside
      12 distinct Bakta products of which 11 are real names
      ('Phage conserved hypothetical protein C-terminal domain-containing
      protein'), and 'conserved' inside 37 of which nearly all are real
      ('Conserved virulence factor B').

  compare(...)
      Given two products that both carry a claim, are they the same claim?

The label is disagreement between two software products about a name. It is
not a claim that either name is correct, and it is not a claim about what the
protein does.
"""

import re
import unicodedata

# --- placeholders ----------------------------------------------------------
#
# Enumerated from the observed data (87,859 paired CDS). Counts at the time of
# writing are in results/metrics/13_name_rules.json, regenerated on every run.
# The last three are declared but not observed in this panel; they are kept so
# the rule is complete rather than fitted to these 25 genomes.
PLACEHOLDERS = frozenset({
    "hypothetical protein",             # bakta 10,930   prokka 34,776
    "putative protein",                 # bakta      0   prokka    530
    "protein",                          # bakta      7   prokka     45
    "conserved protein",                # bakta     12   prokka      0
    "uncharacterized protein",          # bakta      1   prokka      0
    "putative uncharacterized protein", # bakta      0   prokka      1
    "uncharacterised protein",          # not observed
    "unknown protein",                  # not observed
    "unnamed protein product",          # not observed
})

# Kept as NAMED on purpose. Each is a database match carrying a structural
# claim, so treating it as uninformative would be our judgement rather than
# the tool's. They are near-deterministic disagreements and their counts are
# reported as declared subgroups in 12_content_cohort.json so that nobody
# discovers them after the fact.
NAMED_BUT_GENERIC = {
    "uncharacterized_orf_name": re.compile(
        r"(?i)^uncharacteri[sz]ed protein [A-Za-z0-9_.\-]+$"),
    "domain_containing": re.compile(r"(?i)domain-containing protein$"),
    "duf": re.compile(r"\bDUF\d+"),
}

# --- normalisation ---------------------------------------------------------

HEDGES = ("putative", "probable", "possible", "predicted", "conserved")

# DnaA, GyrA, FabG, GroES, TssE1, RpoB. Requires at least four characters with
# an internal capital, so bare words ('Porin', 'Integrase') and three-letter
# tokens are never mistaken for symbols.
GENE_SHAPE = re.compile(r"^[A-Za-z]{3}[A-Z][A-Za-z0-9]{0,2}$")

TRAIL_BRACKET = re.compile(r"\s*[\(\[][^\(\)\[\]]*[\)\]]\s*$")

# Dropped by the 'loose' sensitivity check only.
GENERIC_TOKENS = frozenset({
    "protein", "putative", "family", "domain", "containing", "probable",
    "possible", "predicted", "conserved", "subunit", "type", "like", "related",
})

LEVELS = ("strict", "primary", "loose")


def casefold_ws(s):
    """Case-fold, collapse whitespace, drop a trailing full stop."""
    return re.sub(r"\s+", " ", unicodedata.normalize("NFKC", s).strip()).rstrip(".").lower()


def is_placeholder(product):
    return casefold_ws(product) in PLACEHOLDERS


def strip_trailing_brackets(s):
    prev = None
    while prev != s:
        prev, s = s, TRAIL_BRACKET.sub("", s).strip()
    return s


def symbol_candidates(product_a, gene_a, product_b, gene_b):
    """Trailing tokens treated as gene symbols FOR THIS PAIR, from both tools.

    Built once per pair and applied to both sides. Keying it on each record's
    own gene= attribute instead is wrong and was caught in review: Bakta emits
    gene= on 19,216 of 87,859 CDS against Prokka's 48,201, so a per-record rule
    fires on one side only and drives identical raw strings apart. Before the
    fix, 'GTPase Era' vs 'GTPase Era' was scored as a disagreement 23 times.
    assert_symmetric() below is the regression test for exactly that.
    """
    out = set()
    for g in (gene_a, gene_b):
        if g and g.strip():
            out.add(g.strip().lower())
    for p in (product_a, product_b):
        toks = strip_trailing_brackets(unicodedata.normalize("NFKC", p).strip()).split()
        if len(toks) >= 2:
            last = toks[-1].rstrip(",;")
            if GENE_SHAPE.match(last):
                out.add(last.lower())
    return frozenset(out)


def strip_leading_hedges(s):
    changed = True
    while changed:
        changed = False
        for h in HEDGES:
            m = re.match(rf"(?i)^{h}\b[\s,:-]+", s)
            if m:
                s, changed = s[m.end():], True
    return s


def normalise(product, symbols=frozenset(), level="primary"):
    """One product string -> its comparable form under `level`.

    primary:  NFKC -> strip trailing (...)/[...] -> strip one trailing gene
              symbol -> strip leading hedges -> case-fold -> every run of
              non-alphanumerics becomes one space -> collapse.
    strict:   case-fold and whitespace only. Nothing stripped.
    loose:    primary, then drop generic tokens and compare as an unordered
              set of the remainder.
    """
    if level not in LEVELS:
        raise ValueError(f"unknown level {level!r}; expected one of {LEVELS}")
    s = unicodedata.normalize("NFKC", product).strip()
    if level == "strict":
        return re.sub(r"\s+", " ", s).rstrip(".").lower()
    s = strip_trailing_brackets(s)
    toks = s.split()
    if len(toks) >= 2 and toks[-1].rstrip(",;").lower() in symbols:
        s = " ".join(toks[:-1]).rstrip(" ,;")
    s = strip_leading_hedges(s)
    s = re.sub(r"[^a-z0-9]+", " ", s.lower())
    s = re.sub(r"\s+", " ", s).strip()
    if level == "loose":
        kept = [t for t in s.split() if t not in GENERIC_TOKENS]
        return " ".join(sorted(set(kept))) if kept else s
    return s


def compare(product_bakta, gene_bakta, product_prokka, gene_prokka, level="primary"):
    """-> (agree, normalised_bakta, normalised_prokka)."""
    symbols = (frozenset() if level == "strict"
               else symbol_candidates(product_bakta, gene_bakta,
                                      product_prokka, gene_prokka))
    a = normalise(product_bakta, symbols, level)
    b = normalise(product_prokka, symbols, level)
    return a == b, a, b


def assert_symmetric(pairs):
    """Identical raw product strings must never be scored as a disagreement.

    Fails loudly rather than returning a count. A normalisation rule that can
    split a string from itself is broken, and a broken rule that merely logs a
    warning would still write a label file.
    """
    broken = []
    for p in pairs:
        if p["product_bakta"] != p["product_prokka"]:
            continue
        for level in LEVELS:
            agree, a, b = compare(p["product_bakta"], p["gene_bakta"],
                                  p["product_prokka"], p["gene_prokka"], level)
            if not agree:
                broken.append({"level": level, "product": p["product_bakta"],
                               "normalised_bakta": a, "normalised_prokka": b})
    if broken:
        for b in broken[:10]:
            print(f"  {b['level']}: {b['product']!r} -> {b['normalised_bakta']!r} "
                  f"vs {b['normalised_prokka']!r}")
        raise SystemExit(
            f"FATAL: {len(broken)} pairs with identical raw product strings were "
            "scored as disagreements. The normalisation rule is not symmetric "
            "between the tools. Refusing to write a label built on it."
        )
    return len(pairs)

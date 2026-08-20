#!/usr/bin/env bash
# Resolve the annotation executables and the Bakta database.
#
# Preference order, highest first:
#   1. whatever is on PATH  (a conda env created from environment.yml)
#   2. the locally cached bactopia tool environments on this machine
#
# Whichever is resolved is what gets recorded in provenance/. The published
# reproduction recipe stays environment.yml; this file only lets the pipeline
# run on a machine where those executables already exist elsewhere, without
# spending several gigabytes duplicating them.
#
# Source this file. Do not execute it.

set -euo pipefail

BPDP_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export BPDP_ROOT

_bactopia_bakta="/home/qasim/.bactopia/conda/bioconda--bakta-1.12.0"
_bactopia_prokka="/home/qasim/.bactopia/conda/bioconda--prokka-1.15.6"
_bactopia_amrfinder="/home/qasim/.bactopia/conda/bioconda--ncbi-amrfinderplus-4.2.7"

# --- bakta -----------------------------------------------------------------
if command -v bakta >/dev/null 2>&1; then
  BAKTA_BIN="$(command -v bakta)"
  BAKTA_DB_BIN="$(command -v bakta_db)"
  BAKTA_PREFIX="$(dirname "$(dirname "$BAKTA_BIN")")"
elif [[ -x "$_bactopia_bakta/bin/bakta" ]]; then
  BAKTA_BIN="$_bactopia_bakta/bin/bakta"
  BAKTA_DB_BIN="$_bactopia_bakta/bin/bakta_db"
  BAKTA_PREFIX="$_bactopia_bakta"
else
  echo "FATAL: bakta not found on PATH and no cached environment at $_bactopia_bakta" >&2
  echo "       create it with: conda env create -f environment.yml" >&2
  exit 1
fi

# --- prokka ----------------------------------------------------------------
# prokka is a perl script and breaks if PERL5LIB points at another prefix.
# It also shells out to minced, which is Java. An activated conda base env
# exports JAVA_HOME and JAVA_LD_LIBRARY_PATH; prokka's bundled JVM honours
# them, then looks for its own runtime libraries under the wrong prefix and
# dies with 'Could not resolve "ZIP_Open"' before any annotation happens.
# Both are cleared in run_prokka below.
if command -v prokka >/dev/null 2>&1; then
  PROKKA_BIN="$(command -v prokka)"
  PROKKA_PREFIX="$(dirname "$(dirname "$PROKKA_BIN")")"
elif [[ -x "$_bactopia_prokka/bin/prokka" ]]; then
  PROKKA_BIN="$_bactopia_prokka/bin/prokka"
  PROKKA_PREFIX="$_bactopia_prokka"
else
  echo "FATAL: prokka not found on PATH and no cached environment at $_bactopia_prokka" >&2
  echo "       create it with: conda env create -f environment.yml" >&2
  exit 1
fi

# --- barrnap / nhmmer libgsl -----------------------------------------------
# prokka's cached env ships libgsl.so.27, but the nhmmer binary in it was
# linked against libgsl.so.25. barrnap shells out to nhmmer, so rRNA detection
# fails at runtime and prokka reports "Found 0 rRNAs" for every genome while
# still exiting 0. That is a broken-install artefact, and a damaging one here:
# with no rRNA from Prokka, every rRNA locus Bakta calls looks like a
# tool disagreement when it is really a missing shared library.
#
# Resolve the missing soname from a sibling cached env, appended AFTER prokka's
# own lib dir so prokka's libraries still win every other lookup.
PROKKA_LD_EXTRA=""
if [[ ! -e "$PROKKA_PREFIX/lib/libgsl.so.25" ]]; then
  for _cand in /home/qasim/.bactopia/conda/*/lib; do
    if [[ -e "$_cand/libgsl.so.25" ]]; then PROKKA_LD_EXTRA="$_cand"; break; fi
  done
  if [[ -z "$PROKKA_LD_EXTRA" ]]; then
    echo "WARNING: libgsl.so.25 not found; prokka will report 0 rRNAs" >&2
  fi
fi

# --- amrfinder (bakta_db download calls amrfinder_update) -------------------
if command -v amrfinder >/dev/null 2>&1; then
  AMRFINDER_PREFIX="$(dirname "$(dirname "$(command -v amrfinder)")")"
elif [[ -x "$_bactopia_amrfinder/bin/amrfinder" ]]; then
  AMRFINDER_PREFIX="$_bactopia_amrfinder"
else
  AMRFINDER_PREFIX=""
fi

# --- bakta database --------------------------------------------------------
# db-light, not db-full. 
BAKTA_DB="${BAKTA_DB:-$BPDP_ROOT/db/db-light}"

export BAKTA_BIN BAKTA_DB_BIN BAKTA_PREFIX PROKKA_BIN PROKKA_PREFIX AMRFINDER_PREFIX BAKTA_DB PROKKA_LD_EXTRA

# Run bakta with its own prefix first on PATH, isolated from any leaked PERL5LIB.
run_bakta()    { env -u PERL5LIB PATH="$BAKTA_PREFIX/bin:${AMRFINDER_PREFIX:+$AMRFINDER_PREFIX/bin:}$PATH" "$BAKTA_BIN" "$@"; }
run_bakta_db() { env -u PERL5LIB PATH="$BAKTA_PREFIX/bin:${AMRFINDER_PREFIX:+$AMRFINDER_PREFIX/bin:}$PATH" "$BAKTA_DB_BIN" "$@"; }
run_prokka()   { env -u PERL5LIB -u JAVA_HOME -u JAVA_LD_LIBRARY_PATH LD_LIBRARY_PATH="$PROKKA_PREFIX/lib${PROKKA_LD_EXTRA:+:$PROKKA_LD_EXTRA}" PATH="$PROKKA_PREFIX/bin:$PATH" "$PROKKA_BIN" "$@"; }

#!/usr/bin/env bash
# run_add.sh -- build the combined contextual alignments.
#
# Run from anywhere; paths are resolved relative to this script, so the same
# command works in the container and on the Mac.
#
# The Nextclade REFERENCE is added alongside the Melghat genomes. That is what
# makes the coordinate map exact: reference positions are read directly off the
# reference row of the combined alignment, not inferred from sample genomes.
#
# Requires: mafft on PATH. Set INFORM_DIR to the directory holding the two
# INFORM alignments (defaults to a sibling 'inform' directory).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ="$(cd "$HERE/.." && pwd)"
INFORM_DIR="${INFORM_DIR:-$HERE/inform}"
THREADS="${THREADS:-8}"

command -v mafft >/dev/null || { echo "ERROR: mafft not found on PATH" >&2; exit 1; }

# INFORM_PLUS=1 uses the panel extended with contemporaneous genomes prepared by
# scripts/21_prepare_context.py (see ADD_2021_22_CONTEXT.md). That file is
# UNALIGNED, so mafft realigns the whole panel rather than adding to an existing
# alignment -- slower, and the combined column numbering changes, so
# 15_contextual_prep.py must be rerun afterwards.
for ty in A B; do
  if [ "${INFORM_PLUS:-0}" = "1" ]; then
    INFORM="$INFORM_DIR/RSV.$ty.aln.wGB.noRC.knownDate.plus2021.fasta"
    [ -s "$INFORM" ] || { echo "ERROR: INFORM_PLUS=1 but $INFORM is missing -- run scripts/21_prepare_context.py" >&2; exit 1; }
    echo "[$ty] using EXTENDED panel (unaligned; full realignment)"
  else
    INFORM="$INFORM_DIR/RSV.$ty.aln.wGB.noRC.knownDate.fasta"
  fi
  MELGHAT="$PROJ/data/nextclade_input_RSV-$ty.fasta"
  REF="$HERE/reference/ref_RSV-$ty.fasta"
  for f in "$INFORM" "$MELGHAT" "$REF"; do
    [ -s "$f" ] || { echo "ERROR: missing or empty input: $f" >&2; exit 1; }
  done
  echo "[$ty] building query set (reference + Melghat)"
  cat "$REF" "$MELGHAT" > "$HERE/query_RSV-$ty.fasta"
  echo "[$ty] mafft --add  ($(grep -c '^>' "$HERE/query_RSV-$ty.fasta") queries into $(grep -c '^>' "$INFORM") contextual)"
  if [ "${INFORM_PLUS:-0}" = "1" ]; then
    cat "$INFORM" "$HERE/query_RSV-$ty.fasta" > "$HERE/all_RSV-$ty.fasta"
    mafft --auto --thread "$THREADS" "$HERE/all_RSV-$ty.fasta" \
          > "$HERE/combined_RSV-$ty.fasta" 2> "$HERE/mafft_RSV-$ty.log"
  else
    mafft --add "$HERE/query_RSV-$ty.fasta" --thread "$THREADS" "$INFORM" \
          > "$HERE/combined_RSV-$ty.fasta" 2> "$HERE/mafft_RSV-$ty.log"
  fi
done

echo "[*] checksums"
( cd "$HERE" && sha256sum combined_RSV-A.fasta combined_RSV-B.fasta \
      query_RSV-A.fasta query_RSV-B.fasta \
      reference/ref_RSV-A.fasta reference/ref_RSV-B.fasta \
      reference/annotation_RSV-A.gff3 reference/annotation_RSV-B.gff3 > SHA256SUMS.txt )
# record the INFORM inputs with paths relative to $HERE, so `sha256sum -c`
# run from this directory resolves them (a bare filename does not)
( cd "$HERE" && for f in "$INFORM_DIR"/RSV.[AB].aln.wGB.noRC.knownDate.fasta; do
      rel="${f#$HERE/}"; sha256sum "$rel" 2>/dev/null || sha256sum "$f"
  done >> SHA256SUMS.txt )
echo "[*] done -- see $HERE/SHA256SUMS.txt"

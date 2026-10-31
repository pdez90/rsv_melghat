#!/usr/bin/env bash
# run_tree.sh -- exploratory maximum-likelihood trees from the combined alignments.
#
# SUPPORT VALUES: FastTree's node values are SH-like LOCAL support (Shimodaira-
# Hasegawa-like test on each split, computed by resampling site likelihoods).
# They are NOT whole-tree bootstrap proportions and are not interchangeable with
# them: they test each split against the two alternative nearest-neighbour
# topologies only, holding the rest of the tree fixed. Price, Dehal & Arkin 2010,
# PLoS ONE 5:e9490 (PMC2835736). Treat these trees as exploratory.
set -euo pipefail
# SHA256SUMS: REPLACE the entry for each rebuilt file instead of appending one.
# Appending (the previous behaviour) leaves two different digests under the same
# filename after a rebuild, so the checksum file contradicts itself and cannot be
# used to verify anything.
record_sums() {
  local d="$1"; shift
  ( cd "$d" || exit 1
    touch SHA256SUMS.txt
    for f in "$@"; do
      grep -vF "  ${f}" SHA256SUMS.txt > SHA256SUMS.next || true
      mv SHA256SUMS.next SHA256SUMS.txt
    done
    sha256sum "$@" >> SHA256SUMS.txt )
}

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# PRE-FLIGHT: resolve the binary BEFORE touching any output file.
#
# This block exists because of a real failure on 2026-10-31. The build line was
#     "$FT" ... > tree_RSV-A.nwk
# and FastTree was not installed. The shell performs the redirection FIRST, so
# the tree was truncated to 0 bytes before the missing command was even looked
# up; `set -e` cannot help, because nothing had failed yet at that point. A
# re-freeze then archived the empty file and reported success. So: check the
# binary first, and build to a temporary file that is only moved into place
# once the build has actually succeeded.
FT=""
for c in fasttree FastTree FastTreeDbl fasttreeDbl; do
  command -v "$c" >/dev/null 2>&1 && { FT="$c"; break; }
done
if [ -z "$FT" ]; then
  echo "ERROR: no FastTree binary on PATH (looked for: fasttree, FastTree, FastTreeDbl, fasttreeDbl)." >&2
  echo "       Nothing has been modified. Install it first, e.g.  brew install fasttree" >&2
  exit 127
fi
echo "[*] using $FT -- $("$FT" 2>&1 | head -1 || true)"
for ty in A B; do
  IN="$HERE/combined_RSV-$ty.fasta"
  [ -s "$IN" ] || { echo "ERROR: missing $IN -- run run_add.sh first" >&2; exit 1; }
  echo "[$ty] FastTree -gtr -nt  ($(grep -c '^>' "$IN") sequences)"
  "$FT" -gtr -nt -log "$HERE/fasttree_RSV-$ty.log" "$IN" > "$HERE/tree_RSV-$ty.nwk.partial" 2> "$HERE/fasttree_RSV-$ty.err"
  [ -s "$HERE/tree_RSV-$ty.nwk.partial" ] || { echo "ERROR: $FT produced no tree for $ty; see fasttree_RSV-$ty.err. Existing tree left untouched." >&2; rm -f "$HERE/tree_RSV-$ty.nwk.partial"; exit 1; }
  mv "$HERE/tree_RSV-$ty.nwk.partial" "$HERE/tree_RSV-$ty.nwk"
done
record_sums "$HERE" tree_RSV-A.nwk tree_RSV-B.nwk
# BUILD PROVENANCE: record the digest of the exact input this tree was built from,
# beside the tree. A taxon-set check cannot establish this -- two alignments can
# carry identical sequence names but different nucleotides or columns -- so the
# input digest is written here and verified later by 28_dependency_check.R.
for ty in A B; do
  ( cd "$HERE" && printf 'alignment_sha256=%s\nalignment_file=%s\nbuilt=%s\nbuilder=run_tree.sh\n' \
      "$(sha256sum "combined_RSV-$ty.fasta" | cut -d" " -f1)" \
      "combined_RSV-$ty.fasta" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
      > "tree_RSV-$ty.nwk.provenance" )
done
echo "[*] done"

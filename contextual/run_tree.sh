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
command -v fasttree >/dev/null && FT=fasttree || FT=FastTree
for ty in A B; do
  IN="$HERE/combined_RSV-$ty.fasta"
  [ -s "$IN" ] || { echo "ERROR: missing $IN -- run run_add.sh first" >&2; exit 1; }
  echo "[$ty] FastTree -gtr -nt  ($(grep -c '^>' "$IN") sequences)"
  "$FT" -gtr -nt -log "$HERE/fasttree_RSV-$ty.log" "$IN" > "$HERE/tree_RSV-$ty.nwk" 2> "$HERE/fasttree_RSV-$ty.err"
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

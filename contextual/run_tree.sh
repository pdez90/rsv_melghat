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
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
command -v fasttree >/dev/null && FT=fasttree || FT=FastTree
for ty in A B; do
  IN="$HERE/combined_RSV-$ty.fasta"
  [ -s "$IN" ] || { echo "ERROR: missing $IN -- run run_add.sh first" >&2; exit 1; }
  echo "[$ty] FastTree -gtr -nt  ($(grep -c '^>' "$IN") sequences)"
  "$FT" -gtr -nt -log "$HERE/fasttree_RSV-$ty.log" "$IN" > "$HERE/tree_RSV-$ty.nwk" 2> "$HERE/fasttree_RSV-$ty.err"
done
( cd "$HERE" && sha256sum tree_RSV-A.nwk tree_RSV-B.nwk >> SHA256SUMS.txt )
echo "[*] done"

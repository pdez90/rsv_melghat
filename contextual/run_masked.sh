#!/usr/bin/env bash
# run_masked.sh -- rebuild trees with uncertain columns removed.
# Mask = MASK-HVR2 from 17_repeat_region.py: low-occupancy columns PLUS the
# subgroup-specific G HVR2 window (RSV-A aa 230-322, RSV-B aa 215-311), with every
# boundary derived through the reference-to-alignment map and verified.
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

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; PROJ="$(cd "$HERE/.." && pwd)"
# PRE-FLIGHT: resolve the binary BEFORE touching any output file. See the note
# in run_tree.sh: on 2026-10-31 a missing FastTree truncated a tree to 0 bytes,
# because the shell performs the redirection before looking the command up.
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
for ty in RSV-A RSV-B; do
  python3 -I - "$HERE/combined_$ty.fasta" "$PROJ/output/17_mask_HVR2_$ty.txt" \
      "$HERE/masked_$ty.fasta" <<'PY'
import sys
inp, maskf, outp = sys.argv[1:4]
mask = {int(l) for l in open(maskf) if l.strip()}
def it(p):
    k, buf = None, []
    for line in open(p):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if k: yield k, "".join(buf)
            k, buf = line[1:], []
        elif k is not None: buf.append(line.strip())
    if k: yield k, "".join(buf)
with open(outp, "w") as o:
    n = 0
    for name, s in it(inp):
        keep = "".join(c for i, c in enumerate(s, 1) if i not in mask)
        o.write(f">{name}\n{keep}\n"); n += 1
    print(f"  {outp}: {n} seqs, {len(keep)} cols (dropped {len(mask)})")
PY
  echo "[$ty] FastTree on masked alignment"
  "$FT" -gtr -nt -log "$HERE/fasttree_masked_$ty.log" "$HERE/masked_$ty.fasta" \
      > "$HERE/tree_masked_$ty.nwk.partial" 2> "$HERE/fasttree_masked_$ty.err"
  [ -s "$HERE/tree_masked_$ty.nwk.partial" ] || { echo "ERROR: $FT produced no tree for masked $ty; see fasttree_masked_$ty.err. Existing tree left untouched." >&2; rm -f "$HERE/tree_masked_$ty.nwk.partial"; exit 1; }
  mv "$HERE/tree_masked_$ty.nwk.partial" "$HERE/tree_masked_$ty.nwk"
done
record_sums "$HERE" tree_masked_RSV-A.nwk tree_masked_RSV-B.nwk
# BUILD PROVENANCE: the masked trees depend on the alignment AND the mask, so both
# digests are recorded. Matching tip labels cannot detect a superseded mask.
for ty in RSV-A RSV-B; do
  ( cd "$HERE" && printf 'alignment_sha256=%s\nalignment_file=%s\nmask_sha256=%s\nmask_file=%s\nbuilt=%s\nbuilder=run_masked.sh\n' \
      "$(sha256sum "combined_$ty.fasta" | cut -d" " -f1)" "combined_$ty.fasta" \
      "$(sha256sum "$PROJ/output/17_mask_HVR2_$ty.txt" | cut -d" " -f1)" \
      "output/17_mask_HVR2_$ty.txt" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
      > "tree_masked_$ty.nwk.provenance" )
done
echo "[*] masked trees done"

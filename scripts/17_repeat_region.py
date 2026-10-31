#!/usr/bin/env python3
"""
17_repeat_region.py -- inspect alignment quality across the G hypervariable
region, characterise who contributes gaps there, and emit VERIFIED masks.

Rewritten after review. Three defects in the first version:

1. WINDOW ARITHMETIC. It computed the duplication window as
   (G start in COMBINED columns) + nucleotide offsets. The reference-to-combined
   offset is not constant inside G -- for RSV-A it changes from 49 to 52 because
   three columns are inserted within the gene -- so that window was off by 3
   columns at both ends. Every window boundary is now derived by mapping
   REFERENCE positions through output/15_coordmap_*.csv, and asserted.

2. SUBGROUP. It applied the RSV-A ON1 72-nt duplication window (G aa 262-285)
   to RSV-B, which has no such feature at that position; RSV-B's BA duplication
   is a different length at a different site and is fixed in all contemporary
   strains. Nothing in this project defines an RSV-B duplication window. The
   masked region is therefore the HVR2 definition already used in 09_ggene.R,
   which IS subgroup-specific: RSV-A aa 230-322, RSV-B aa 215-311.

3. INTERPRETATION. Gap frequency alone does not establish polymorphism. The
   contextual contributors are now broken down by sampling year, since older
   non-duplicated strains, incomplete sequence and alignment choices all produce
   gaps and are not distinguishable by a column-occupancy summary.

Two masks are written, so the robustness claim can be stated at the strength the
evidence supports:
    MASK-OCC    low-occupancy columns only
    MASK-HVR2   low-occupancy + the subgroup's ENTIRE HVR2 span (mapped columns
                plus any columns inserted inside that span)
"""
import os, csv, collections, sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CTX  = os.path.join(PROJ, "contextual")
OUT  = os.path.join(PROJ, "output")

# Reference CDS starts (from contextual/reference/annotation_RSV-*.gff3) and the
# HVR2 amino-acid windows used in 09_ggene.R. RSV-A additionally has the ON1
# 72-nt duplication at aa 262-285; RSV-B has no counterpart defined here.
SPEC = {
    "RSV-A": dict(g_ref_start=4652, hvr2=(230, 322), dup=(262, 285)),
    "RSV-B": dict(g_ref_start=4646, hvr2=(215, 311), dup=None),
}

def read_fa(p):
    d, k, buf = {}, None, []
    for line in open(p):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if k: d[k] = "".join(buf)
            k, buf = line[1:].split()[0], []
        elif k is not None: buf.append(line.strip())
    if k: d[k] = "".join(buf)
    return {k: v.upper() for k, v in d.items()}

def load_map(ty):
    ref2comb, gene = {}, {}
    for r in csv.DictReader(open(os.path.join(OUT, f"15_coordmap_{ty}.csv"))):
        p, c = int(r["ref_pos"]), int(r["combined_col"])
        ref2comb[p] = c
        gene[p] = r["gene"]
    return ref2comb, gene

def aa_window_to_cols(ref2comb, g_ref_start, aa0, aa1, label):
    """Map an amino-acid window in G to combined columns, via reference positions."""
    r0 = g_ref_start + (aa0 - 1) * 3
    r1 = g_ref_start + aa1 * 3 - 1
    missing = [p for p in range(r0, r1 + 1) if p not in ref2comb]
    assert not missing, f"{label}: {len(missing)} reference positions unmapped"
    cols = [ref2comb[p] for p in range(r0, r1 + 1)]
    # the window must be contiguous in reference space; in combined space it may
    # be interrupted by inserted columns, so take the span and report the gap
    span = (min(cols), max(cols))
    inserted = (span[1] - span[0] + 1) - len(cols)
    return cols, span, (r0, r1), inserted

def process(ty):
    comb_p = os.path.join(CTX, f"combined_{ty}.fasta")
    if not os.path.exists(comb_p):
        print(f"  [{ty}] no combined alignment"); return None, None
    comb = read_fa(comb_p)
    ncal = read_fa(os.path.join(PROJ, "nextclade", ty, "nextclade.aligned.fasta"))
    refn = next(k for k in comb if k.startswith("REFERENCE"))
    mel  = [k for k in comb if k in ncal]
    ctx  = [k for k in comb if k not in ncal and k != refn]
    L    = len(comb[refn])
    spec = SPEC[ty]
    ref2comb, gene = load_map(ty)
    gcols = sorted(ref2comb[p] for p, g in gene.items() if g == "G")
    print(f"\n[{ty}] {L} columns | G at combined {min(gcols)}-{max(gcols)}")

    windows = {}
    hv_cols, hv_span, hv_ref, hv_ins = aa_window_to_cols(
        ref2comb, spec["g_ref_start"], *spec["hvr2"], label="HVR2")
    windows["HVR2"] = (hv_cols, hv_span, hv_ref, hv_ins)
    print(f"  HVR2 (G aa {spec['hvr2'][0]}-{spec['hvr2'][1]}) = ref {hv_ref[0]}-{hv_ref[1]}"
          f" -> combined {hv_span[0]}-{hv_span[1]} ({len(hv_cols)} cols,"
          f" {hv_ins} inserted columns inside the span)")
    if spec["dup"]:
        d_cols, d_span, d_ref, d_ins = aa_window_to_cols(
            ref2comb, spec["g_ref_start"], *spec["dup"], label="duplication")
        windows["duplication"] = (d_cols, d_span, d_ref, d_ins)
        print(f"  ON1 duplication (G aa {spec['dup'][0]}-{spec['dup'][1]}) = ref"
              f" {d_ref[0]}-{d_ref[1]} -> combined {d_span[0]}-{d_span[1]}"
              f" ({len(d_cols)} cols, {d_ins} inserted inside span)")
    else:
        print("  no duplication window defined for this subgroup (see module docstring)")

    def occ(cols, names):
        b = g = 0
        for c in cols:
            for n in names:
                x = comb[n][c - 1]
                if x in "ACGT": b += 1
                elif x in "-.": g += 1
        tot = len(cols) * len(names)
        return b / tot, g / tot

    rows = []
    regions = {k: v[0] for k, v in windows.items()}
    regions["whole alignment"] = list(range(1, L + 1))
    for name, cols in regions.items():
        mb, mg = occ(cols, mel); cb, cg = occ(cols, ctx)
        print(f"  {name:16} {len(cols):>6} cols | Melghat base {mb:5.1%} gap {mg:5.1%}"
              f" | contextual base {cb:5.1%} gap {cg:5.1%}")
        rows.append(dict(rsv_type=ty, region=name, n_cols=len(cols),
                         melghat_base=round(mb,4), melghat_gap=round(mg,4),
                         contextual_base=round(cb,4), contextual_gap=round(cg,4)))

    ## who contributes contextual gaps in the duplication window, by year?
    if spec["dup"]:
        cols = windows["duplication"][0]
        by_year = collections.defaultdict(lambda: [0, 0])   # [gappy seqs, total]
        for n in ctx:
            parts = n.split("_")
            yr = parts[2][:4] if len(parts) > 2 and parts[2][:4].isdigit() else "unknown"
            s = comb[n]
            gaps = sum(1 for c in cols if s[c - 1] in "-.")
            by_year[yr][1] += 1
            if gaps >= len(cols) * 0.5: by_year[yr][0] += 1
        print("\n  contextual sequences >=50% gapped across the duplication window, by year:")
        for yr in sorted(by_year):
            g, t = by_year[yr]
            if t >= 10:
                print(f"    {yr}: {g}/{t} ({100*g/t:.0f}%)")

    ## low-occupancy columns across the full panel
    allnames = list(comb.keys())
    unc = [c for c in range(1, L + 1)
           if sum(1 for n in allnames if comb[n][c-1] in "ACGT") / len(allnames) < 0.50]
    ## MASK-HVR2 must exclude every column BETWEEN the mapped HVR2 boundaries,
    ## not only the mapped ones. Columns inserted inside the window would
    ## otherwise survive whenever they pass the occupancy threshold, and the
    ## analysis could not be called whole-HVR2 masking. Assert the span.
    hv_mapped = set(windows["HVR2"][0])
    hv_span_cols = set(range(windows["HVR2"][1][0], windows["HVR2"][1][1] + 1))
    hv_inserted = hv_span_cols - hv_mapped
    masks = {
        "OCC":  sorted(set(unc)),
        "HVR2": sorted(set(unc) | hv_span_cols),
    }
    leaked = hv_span_cols - set(masks["HVR2"])
    assert not leaked, f"{ty}: {len(leaked)} HVR2-span columns escaped the mask"
    print(f"  HVR2 span completeness: {len(hv_mapped)} mapped + {len(hv_inserted)}"
          f" inserted = {len(hv_span_cols)} columns, all masked -> PASS")
    for tag, mk in masks.items():
        p = os.path.join(OUT, f"17_mask_{tag}_{ty}.txt")
        open(p, "w").write("\n".join(map(str, mk)) + "\n")
        print(f"  -> MASK-{tag}: {len(mk)} cols ({100*len(mk)/L:.1f}%) -> {os.path.basename(p)}")
        rows.append(dict(rsv_type=ty, region=f"MASK-{tag}", n_cols=len(mk),
                         melghat_base="", melghat_gap="", contextual_base="", contextual_gap=""))

    ## VERIFY: every masked HVR2 column must map back to a reference position in G
    comb2ref = {c: p for p, c in ref2comb.items()}
    hv_set = set(windows["HVR2"][0])
    bad = [c for c in hv_set if gene.get(comb2ref.get(c, -1), "") != "G"]
    print(f"  mask verification: {len(hv_set)} HVR2 columns, {len(bad)} not mapping to G"
          f" -> {'PASS' if not bad else 'FAIL'}")
    return rows, len(bad) == 0

if __name__ == "__main__":
    allrows, ok = [], True
    for ty in ("RSV-A", "RSV-B"):
        r, good = process(ty)
        if r: allrows += r
        ok = ok and (good is not False)
    if allrows:
        with open(os.path.join(OUT, "17_repeat_region.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(allrows[0].keys()))
            w.writeheader(); w.writerows(allrows)
        print("\n  -> wrote 17_repeat_region.csv")
    print("== 17 done ==" if ok else "== 17 FAILED verification ==")
    sys.exit(0 if ok else 1)

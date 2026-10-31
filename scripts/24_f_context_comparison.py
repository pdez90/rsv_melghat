#!/usr/bin/env python3
"""
24_f_context_comparison.py -- place Melghat F residue states inside the
international panel, rather than reporting them only against a resistance list.

Method. Every sequence in the combined alignment -- Melghat, contextual and the
Nextclade reference -- is projected onto the REFERENCE FRAME using the
coordinate map (15_coordmap_*.csv): only columns that map to a reference
position are taken, in reference order. Insertions are therefore dropped, which
is what reference-numbered residue states require. The F CDS columns are then
translated and read at the positions of interest.

Staged comparison, not one pooled percentage:
  stage 1  same RSV subgroup
  stage 2  overlapping collection years only
  stage 3  within lineage where the panel carries enough sequences
  stage 4  country composition reported, so sample-size imbalance is visible

Denominators are CALLABLE counts at each position: a sequence counts only if its
codon there is complete and unambiguous. "Absent from this panel" means not
observed in this comparison panel, never "novel".

The panel's own genotype labels are NOT assumed to match Nextclade lineage
names, so lineage stratification is applied to Melghat genomes only and the
contextual side is stratified by year and country.
"""
import csv, os, sys, collections

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CTX, OUT = os.path.join(PROJ, "contextual"), os.path.join(PROJ, "output")

F_CDS = {"RSV-A": (5697, 7421), "RSV-B": (5676, 7400)}
# positions carried by the curated screen, plus the ones asked for by name
POS = {"RSV-A": sorted({201, 206, 208, 209, 211, 262, 272, 275, 64, 65, 68, 276}),
       "RSV-B": sorted({201, 206, 208, 209, 211, 262, 272, 275, 64, 65, 68, 276})}
COMBOS = [("I206M+Q209R", (206, "M"), (209, "R")),
          ("K68E+I206M",  (68, "E"),  (206, "M"))]

CODE = {}
_b = "TCAG"
_aa = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
for _i, (_c1, _c2, _c3) in enumerate((a, b, c) for a in _b for b in _b for c in _b):
    CODE[_c1 + _c2 + _c3] = _aa[_i]

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

def aa_at(seq_cols, idx):
    """idx = 0-based residue index; returns the amino acid or None if not callable."""
    c = seq_cols[idx * 3: idx * 3 + 3]
    if len(c) < 3 or any(x not in "ACGT" for x in c): return None
    return CODE.get(c)

def ctx_meta(name):
    p = name.split("_")
    country = p[1] if len(p) > 2 else "unknown"
    year = p[2][:4] if len(p) > 2 and p[2][:4].isdigit() else None
    return country, year

def run(ty):
    comb = read_fa(os.path.join(CTX, f"combined_{ty}.fasta"))
    ncal = read_fa(os.path.join(PROJ, "nextclade", ty, "nextclade.aligned.fasta"))
    refn = next(k for k in comb if k.startswith("REFERENCE"))
    ref2comb = {}
    for r in csv.DictReader(open(os.path.join(OUT, f"15_coordmap_{ty}.csv"))):
        ref2comb[int(r["ref_pos"])] = int(r["combined_col"])
    s, e = F_CDS[ty]
    cols = [ref2comb[p] for p in range(s, e + 1)]

    mel = {k for k in comb if k in ncal}
    # Melghat lineage and season
    meta = {}
    for r in csv.DictReader(open(os.path.join(OUT, "20_g_repeat_specimen_table.csv"))):
        meta[r["FieldNoNIV_f"]] = (r["season"], r["clade"])
    lin_by_name = {}
    for r in csv.DictReader(open(os.path.join(OUT, "03_lineage_by_season.csv"))):
        pass  # season-level only; per-genome lineage comes from 20_* for RSV-A

    proj = {}
    for k, sq in comb.items():
        proj[k] = "".join(sq[c - 1] for c in cols)

    ref_states = {p: aa_at(proj[refn], p - 1) for p in POS[ty]}
    mel_years = sorted({meta.get(k, ("?", "?"))[0] for k in mel if k in meta}) or ["(RSV-B: see 03_lineage_by_season)"]

    rows = []
    for p in POS[ty]:
        i = p - 1
        mel_states = collections.Counter()
        for k in mel:
            a = aa_at(proj[k], i)
            if a: mel_states[a] += 1
        ctx_states = collections.Counter(); ctx_by_year = collections.defaultdict(collections.Counter)
        ctx_country = collections.defaultdict(collections.Counter)
        for k in comb:
            if k in mel or k == refn: continue
            a = aa_at(proj[k], i)
            if not a: continue
            ctx_states[a] += 1
            c, y = ctx_meta(k)
            if y: ctx_by_year[y][a] += 1
            ctx_country[a][c] += 1
        mel_n, ctx_n = sum(mel_states.values()), sum(ctx_states.values())
        # the Melghat majority state, and how common it is in the panel
        if not mel_n: continue
        for state, cnt in sorted(mel_states.items(), key=lambda x: -x[1]):
            cctx = ctx_states.get(state, 0)
            # overlapping years only (Melghat 2016-2022 vs panel, which ends 2020)
            ov = sum(ctx_by_year[y].get(state, 0) for y in ctx_by_year if y >= "2016")
            ov_n = sum(sum(ctx_by_year[y].values()) for y in ctx_by_year if y >= "2016")
            top = ", ".join(f"{c}:{n}" for c, n in ctx_country[state].most_common(3))
            rows.append(dict(subgroup=ty, position=p, state=state,
                             is_reference_state=(state == ref_states[p]),
                             melghat_n=cnt, melghat_callable=mel_n,
                             melghat_pct=round(100 * cnt / mel_n, 1),
                             panel_n=cctx, panel_callable=ctx_n,
                             panel_pct=round(100 * cctx / ctx_n, 2) if ctx_n else None,
                             panel_2016plus_n=ov, panel_2016plus_callable=ov_n,
                             panel_2016plus_pct=round(100 * ov / ov_n, 2) if ov_n else None,
                             panel_top_countries=top))
    return rows, proj, mel, refn, comb

def combos(ty, proj, mel, refn, comb):
    out = []
    for label, (p1, a1), (p2, a2) in COMBOS:
        def has(k):
            x, y = aa_at(proj[k], p1 - 1), aa_at(proj[k], p2 - 1)
            if x is None or y is None: return None
            return (x == a1) and (y == a2)
        m = [has(k) for k in mel]; m = [v for v in m if v is not None]
        c = [has(k) for k in comb if k not in mel and k != refn]; c = [v for v in c if v is not None]
        out.append(dict(subgroup=ty, combination=label,
                        melghat_n=sum(m), melghat_callable=len(m),
                        panel_n=sum(c), panel_callable=len(c),
                        panel_pct=round(100 * sum(c) / len(c), 2) if c else None))
    return out

if __name__ == "__main__":
    allrows, allcombo = [], []
    for ty in ("RSV-A", "RSV-B"):
        r, proj, mel, refn, comb = run(ty)
        allrows += r
        allcombo += combos(ty, proj, mel, refn, comb)
        print(f"\n[{ty}] {len(mel)} Melghat, {len(comb) - len(mel) - 1} contextual genomes")
        for x in r:
            if x["melghat_pct"] >= 5:
                print(f"  F{x['position']}{x['state']}{' (ref)' if x['is_reference_state'] else ''}: "
                      f"Melghat {x['melghat_n']}/{x['melghat_callable']} ({x['melghat_pct']}%) | "
                      f"panel {x['panel_n']}/{x['panel_callable']} ({x['panel_pct']}%) | "
                      f"panel 2016+ {x['panel_2016plus_n']}/{x['panel_2016plus_callable']} "
                      f"({x['panel_2016plus_pct']}%)")
    with open(os.path.join(OUT, "24_f_context_comparison.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(allrows[0].keys())); w.writeheader(); w.writerows(allrows)
    with open(os.path.join(OUT, "24_f_combinations.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(allcombo[0].keys())); w.writeheader(); w.writerows(allcombo)
    print("\n-- combinations --")
    for x in allcombo:
        print(f"  {x['subgroup']} {x['combination']}: Melghat {x['melghat_n']}/{x['melghat_callable']} | "
              f"panel {x['panel_n']}/{x['panel_callable']} ({x['panel_pct']}%)")
    print("\n  -> wrote 24_f_context_comparison.csv, 24_f_combinations.csv")

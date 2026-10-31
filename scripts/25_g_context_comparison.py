#!/usr/bin/env python3
"""
25_g_context_comparison.py -- apply the SAME alignment-based G definitions to
the contextual panel as to Melghat, restricted to overlapping collection years.

Definitions, identical on both sides and both kept separate:
  apparent repeat absence     >=50% of the 72 columns spanning G aa 262-285 are gaps
  apparent 3' shortening      >=60 gapped columns in G 3' of aa 285

Only RSV-A is analysed: the 72-nt ON1 duplication is an RSV-A feature and
nothing in this project defines an RSV-B counterpart.

Staging. Melghat seasons 2016/17 and 2019/20 overlap the panel, which ends in
April 2020; 2021/22 does not, so for that season the panel gives historical
context only and the comparison is reported but marked as non-overlapping.
Gap-run boundaries are compared so that "the same pattern" means the same
columns, not merely the same summary.
"""
import csv, os, collections

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CTX, OUT = os.path.join(PROJ, "contextual"), os.path.join(PROJ, "output")
TY, G0 = "RSV-A", 4652

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

comb = read_fa(os.path.join(CTX, f"combined_{TY}.fasta"))
ncal = read_fa(os.path.join(PROJ, "nextclade", TY, "nextclade.aligned.fasta"))
refn = next(k for k in comb if k.startswith("REFERENCE"))
ref2comb, gene = {}, {}
for r in csv.DictReader(open(os.path.join(OUT, f"15_coordmap_{TY}.csv"))):
    ref2comb[int(r["ref_pos"])] = int(r["combined_col"]); gene[int(r["ref_pos"])] = r["gene"]
gend = max(p for p, g in gene.items() if g == "G")
dup   = [ref2comb[p] for p in range(G0 + (262 - 1) * 3, G0 + 285 * 3)]
three = [ref2comb[p] for p in range(G0 + 285 * 3, gend + 1)]
mel = {k for k in comb if k in ncal}

season = {r["FieldNoNIV_f"]: r["season"] for r in csv.DictReader(
    open(os.path.join(OUT, "20_g_repeat_specimen_table.csv")))}

def feats(k):
    s = comb[k]
    dg = sum(1 for c in dup if s[c - 1] in "-.")
    tg = sum(1 for c in three if s[c - 1] in "-.")
    return dg / len(dup), tg

def run_boundaries(k):
    """first and last gapped column of the duplication-window gap run, in reference positions"""
    s = comb[k]
    comb2ref = {c: p for p, c in ref2comb.items()}
    gapped = [c for c in dup if s[c - 1] in "-."]
    if not gapped: return None
    return (comb2ref.get(min(gapped)), comb2ref.get(max(gapped)), len(gapped))

def year_of(k):
    p = k.split("_")
    return p[2][:4] if len(p) > 2 and p[2][:4].isdigit() else None

rows = []
# Melghat, by season
by_season = collections.defaultdict(lambda: [0, 0, 0])
for k in mel:
    dgf, tg = feats(k)
    s = season.get(k, "?")
    by_season[s][0] += 1
    if dgf >= 0.50: by_season[s][1] += 1
    if tg >= 60:    by_season[s][2] += 1
# contextual, by year
by_year = collections.defaultdict(lambda: [0, 0, 0])
for k in comb:
    if k in mel or k == refn: continue
    y = year_of(k)
    if not y: continue
    dgf, tg = feats(k)
    by_year[y][0] += 1
    if dgf >= 0.50: by_year[y][1] += 1
    if tg >= 60:    by_year[y][2] += 1

print(f"\n== {TY} G features, same definitions both sides ==")
print(f"  duplication window = reference {G0+(262-1)*3}-{G0+285*3-1} ({len(dup)} cols);"
      f" 3' region = {len(three)} cols\n")
print("  Melghat, by season:")
for s in sorted(by_season):
    n, a, b = by_season[s]
    print(f"    {s}: n={n:3}  repeat absence {a:3} ({100*a/n:5.1f}%)  3' shortening {b:2} ({100*b/n:4.1f}%)")
print("\n  Contextual panel, by collection year (overlapping years marked *):")
for y in sorted(by_year):
    n, a, b = by_year[y]
    if n < 10: continue
    star = "*" if y in ("2016", "2017", "2018", "2019", "2020") else " "
    print(f"   {star}{y}: n={n:4}  repeat absence {a:4} ({100*a/n:5.1f}%)  3' shortening {b:3} ({100*b/n:4.1f}%)")
    rows.append(dict(side="panel", stratum=y, overlapping=(star == "*"), n=n,
                     repeat_absence=a, repeat_absence_pct=round(100*a/n, 2),
                     three_prime=b, three_prime_pct=round(100*b/n, 2)))
for s in sorted(by_season):
    n, a, b = by_season[s]
    rows.append(dict(side="Melghat", stratum=s, overlapping=(s in ("2016/17", "2019/20")), n=n,
                     repeat_absence=a, repeat_absence_pct=round(100*a/n, 2),
                     three_prime=b, three_prime_pct=round(100*b/n, 2)))

# overlapping-years head-to-head
mel_ov = [k for k in mel if season.get(k) in ("2016/17", "2019/20")]
ctx_ov = [k for k in comb if k not in mel and k != refn and (year_of(k) or "") in ("2016","2017","2018","2019","2020")]
def rate(ks, i):
    n = len(ks); c = sum(1 for k in ks if (feats(k)[0] >= 0.50 if i == 0 else feats(k)[1] >= 60))
    return c, n, (100*c/n if n else None)
print("\n  Overlapping years only (Melghat 2016/17 + 2019/20 vs panel 2016-2020):")
for i, lab in ((0, "apparent repeat absence"), (1, "apparent 3' shortening")):
    c1, n1, p1 = rate(mel_ov, i); c2, n2, p2 = rate(ctx_ov, i)
    print(f"    {lab:26} Melghat {c1}/{n1} ({p1:.1f}%)   panel {c2}/{n2} ({p2:.1f}%)")
    rows.append(dict(side="overlapping-years comparison", stratum=lab, overlapping=True,
                     n=n1, repeat_absence=c1, repeat_absence_pct=round(p1, 2),
                     three_prime=c2, three_prime_pct=round(p2, 2)))

# do gap boundaries agree?
print("\n  Gap-run boundaries in the duplication window (reference positions):")
mb = collections.Counter(run_boundaries(k) for k in mel if feats(k)[0] >= 0.50)
cb = collections.Counter(run_boundaries(k) for k in comb
                         if k not in mel and k != refn and feats(k)[0] >= 0.50)
for lab, cnt in (("Melghat", mb), ("panel", cb)):
    print(f"    {lab}:")
    for b, n in cnt.most_common(3):
        print(f"      {b[0]}-{b[1]} ({b[2]} gapped cols): {n}")

with open(os.path.join(OUT, "25_g_context_comparison.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("\n  -> wrote 25_g_context_comparison.csv")

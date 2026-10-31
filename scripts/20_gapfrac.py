#!/usr/bin/env python3
"""Helper for 20_*.R: per-genome gap fraction in the ON1 duplication window and
gap length 3' of it, both via the reference->combined coordinate map."""
import os, csv
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT, CTX = os.path.join(PROJ, "output"), os.path.join(PROJ, "contextual")
def read_fa(p):
    d,k,buf={},None,[]
    for l in open(p):
        l=l.rstrip("\n")
        if l.startswith(">"):
            if k: d[k]="".join(buf)
            k,buf=l[1:].split()[0],[]
        elif k is not None: buf.append(l.strip())
    if k: d[k]="".join(buf)
    return {k:v.upper() for k,v in d.items()}
ty="RSV-A"; G0=4652
ref2comb, gene = {}, {}
for r in csv.DictReader(open(os.path.join(OUT, f"15_coordmap_{ty}.csv"))):
    ref2comb[int(r["ref_pos"])]=int(r["combined_col"]); gene[int(r["ref_pos"])]=r["gene"]
dup=[ref2comb[p] for p in range(G0+(262-1)*3, G0+285*3)]
gend=max(p for p,g in gene.items() if g=="G")
three=[ref2comb[p] for p in range(G0+285*3, gend+1)]
comb=read_fa(os.path.join(CTX,f"combined_{ty}.fasta"))
ncal=read_fa(os.path.join(PROJ,"nextclade",ty,"nextclade.aligned.fasta"))
rows=[]
for k in comb:
    if k not in ncal: continue
    s=comb[k]
    gd=sum(1 for c in dup if s[c-1] in "-.")
    g3=sum(1 for c in three if s[c-1] in "-.")
    rows.append(dict(seqName=k, dup_gap_cols=gd, dup_gap_frac=round(gd/len(dup),4),
                     three_prime_gap=g3))
with open(os.path.join(OUT,"20_gapfrac_RSV-A.csv"),"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(f"  -> 20_gapfrac_RSV-A.csv ({len(rows)} genomes; dup window {len(dup)} cols, 3' {len(three)} cols)")

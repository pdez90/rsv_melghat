#!/usr/bin/env python3
"""
21_prepare_context.py -- normalise NCBI Virus downloads into the INFORM panel's
header convention and merge them with the existing contextual sequences.

Input : contextual/incoming/ncbi_RSV-{A,B}_*.fasta   (FASTA from NCBI Virus,
        defline "Accession |Geo Location |Collection Date")
Output: contextual/inform/RSV.{A,B}.aln.wGB.noRC.knownDate.plus2021.fasta
        (UNALIGNED -- run_add.sh realigns; do not concatenate with the aligned panel)

Quality rules, applied and reported rather than silent:
  * length >= 14000 nt after removing gaps
  * <= 5% N
  * a parseable collection date; year-only is kept but flagged
  * accessions already present in the INFORM panel are dropped
"""
import os, re, sys, glob, collections

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CTX  = os.path.join(PROJ, "contextual")
INC  = os.path.join(CTX, "incoming")
INF  = os.path.join(CTX, "inform")
MIN_LEN, MAX_N = 14000, 0.05

def read_fa(p):
    d, k, buf = {}, None, []
    for line in open(p):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if k: d[k] = "".join(buf)
            k, buf = line[1:], []
        elif k is not None: buf.append(line.strip())
    if k: d[k] = "".join(buf)
    return d

def norm_country(s):
    s = (s or "").split(":")[0].strip()
    return re.sub(r"[^A-Za-z]", "", s) or "unknown"

def norm_date(s):
    s = (s or "").strip()
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", s): return s, "full"
    if re.fullmatch(r"\d{4}-\d{2}", s):       return s, "year-month"
    if re.fullmatch(r"\d{4}", s):             return s, "year-only"
    return None, "unparsed"

def main():
    if not os.path.isdir(INC):
        print(f"no {INC} -- see contextual/ADD_2021_22_CONTEXT.md"); return 1
    for sub in ("A", "B"):
        files = sorted(glob.glob(os.path.join(INC, f"ncbi_RSV-{sub}_*.fasta")))
        if not files:
            print(f"[RSV-{sub}] no incoming files, skipped"); continue
        panel_p = os.path.join(INF, f"RSV.{sub}.aln.wGB.noRC.knownDate.fasta")
        have = set()
        if os.path.exists(panel_p):
            for h in read_fa(panel_p):
                have.add(h.split("_")[0].split(".")[0])
        kept, drop = {}, collections.Counter()
        prec = collections.Counter(); byyear = collections.Counter()
        byctry = collections.Counter()
        for f in files:
            for hdr, seq in read_fa(f).items():
                parts = [x.strip() for x in hdr.split("|")]
                acc = parts[0].split()[0].split(".")[0]
                ctry = norm_country(parts[1] if len(parts) > 1 else "")
                date, p = norm_date(parts[2] if len(parts) > 2 else "")
                s = seq.replace("-", "").replace(".", "").upper()
                if acc in have:            drop["already in panel"] += 1; continue
                if date is None:           drop["unparseable date"] += 1; continue
                if len(s) < MIN_LEN:       drop["too short"] += 1; continue
                if s.count("N")/len(s) > MAX_N: drop["too many Ns"] += 1; continue
                name = f"{acc}_{ctry}_{date}_{sub}"
                if name in kept:           drop["duplicate"] += 1; continue
                kept[name] = s; prec[p] += 1
                byyear[date[:4]] += 1; byctry[ctry] += 1
        print(f"\n[RSV-{sub}] {len(kept)} sequences kept from {len(files)} file(s)")
        for r, n in drop.most_common(): print(f"    dropped, {r}: {n}")
        print(f"    date precision: {dict(prec)}")
        print(f"    by year: {dict(sorted(byyear.items()))}")
        print(f"    top countries: {byctry.most_common(8)}")
        if not kept: continue
        out = os.path.join(INF, f"RSV.{sub}.aln.wGB.noRC.knownDate.plus2021.fasta")
        # the existing panel is ALIGNED; write it ungapped so mafft realigns all of it
        with open(out, "w") as o:
            if os.path.exists(panel_p):
                for h, s in read_fa(panel_p).items():
                    o.write(f">{h}\n{s.replace('-','').replace('.','')}\n")
            for h, s in kept.items():
                o.write(f">{h}\n{s}\n")
        print(f"    -> {os.path.basename(out)} (UNALIGNED; run_add.sh with INFORM_PLUS=1)")
    return 0

if __name__ == "__main__":
    sys.exit(main())

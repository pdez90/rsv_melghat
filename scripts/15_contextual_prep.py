#!/usr/bin/env python3
"""
15_contextual_prep.py -- harmonise Melghat genomes with the INFORM contextual
alignment and emit a validated coordinate map, an insertion accounting, a
date-precision audit and a provenance record.

Rewritten after review. What changed and why:

1. Gene boundaries come from the Nextclade dataset's own reference annotation
   (contextual/reference/annotation_RSV-*.gff3), NOT from the per-sequence CDS
   records in nextclade.gff. Those per-sequence records are in each sample's own
   unaligned coordinates, so taking their modal value and applying it to a
   reference-alignment position mixes two coordinate systems. It gave G =
   4583-5476 for RSV-A where the reference annotation says 4652-5617.

2. The reference sequence is added to the combined alignment by run_add.sh, so
   ref_pos -> combined_col is READ OFF THE REFERENCE ROW rather than inferred
   from sample genomes. There is no voting and no ambiguity to resolve.

3. Every Melghat genome is checked by BASE IDENTITY, not by length agreement.
   Equal ungapped length does not establish correspondence. Genomes that fail
   are counted and named, not silently skipped.

4. Insertions stripped by Nextclade are accounted for per sequence from
   nextclade.csv, rather than asserted in aggregate. Measured result:
   totalInsertions is ZERO for all 241 genomes. An earlier claim that Nextclade
   removed ~270 nt of G-duplication sequence was wrong -- those 270 characters
   are DASHES carried in data/nextclade_input_RSV-A.fasta, which is not strictly
   unaligned. The RSV-A reference (EPI_ISL_412866, England/397/2017) is itself an
   ON1-era genome and already carries the 72-nt G duplication, so Melghat
   genomes have no insertion relative to it.
"""
import os, re, csv, sys, hashlib, collections

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CTX  = os.path.join(PROJ, "contextual")
OUT  = os.path.join(PROJ, "output")
REF_ACC = {"RSV-A": "EPI_ISL_412866", "RSV-B": "EPI_ISL_1653999"}

def read_fa(path, upper=True):
    d, k, buf = {}, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if k is not None: d[k] = "".join(buf)
            k, buf = line[1:].split()[0], []
        elif k is not None:
            buf.append(line.strip())
    if k is not None: d[k] = "".join(buf)
    return {k: (v.upper() if upper else v) for k, v in d.items()}

def sha256(path, n=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(n), b""): h.update(b)
    return h.hexdigest()

def ref_genes(ty):
    """CDS boundaries in reference coordinates, from the dataset annotation."""
    g = {}
    for line in open(os.path.join(CTX, "reference", f"annotation_RSV-{ty[-1]}.gff3")):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        # attribute key differs between the two datasets: RSV-A uses
        # gene_name=, RSV-B uses locus_tag=/Name=
        m = (re.search(r"gene_name=([^;]+)", f[8]) or
             re.search(r"locus_tag=([^;]+)", f[8]) or
             re.search(r"Name=([^;]+)", f[8]))
        if m: g[m.group(1)] = (int(f[3]), int(f[4]))
    return g

def base_cols(seq):
    return [i for i, c in enumerate(seq) if c not in "-."]

def nextclade_insertions(ty):
    """seqName -> (totalInsertions, [(after_ref_pos, inserted_seq), ...])"""
    # Nextclade writes .csv (semicolon) or .tsv (tab) depending on how it was
    # run; accept either rather than assuming one.
    cand = [os.path.join(PROJ, "nextclade", ty, n)
            for n in ("nextclade.csv", "nextclade.tsv")]
    p = next((c for c in cand if os.path.exists(c)), None)
    if p is None:
        print(f"  [{ty}] no nextclade.csv or nextclade.tsv -- "
              f"insertion accounting skipped")
        return {}
    with open(p, newline="") as f:
        # detect on the HEADER LINE only: data rows carry comma-separated
        # substitution lists that swamp a byte-sample heuristic
        header = f.readline(); f.seek(0)
        delim = max((";", ",", "\t"), key=header.count)
        out = {}
        for row in csv.DictReader(f, delimiter=delim):
            ins = []
            raw = (row.get("insertions") or "").strip()
            if raw:
                for item in raw.split(","):
                    if ":" in item:
                        pos, s = item.split(":", 1)
                        try: ins.append((int(pos), s))
                        except ValueError: pass
            try: tot = int(row.get("totalInsertions") or 0)
            except ValueError: tot = 0
            out[row["seqName"]] = (tot, ins)
    return out

def process(ty):
    comb_p = os.path.join(CTX, f"combined_{ty}.fasta")
    if not os.path.exists(comb_p) or os.path.getsize(comb_p) == 0:
        print(f"  [{ty}] combined alignment missing -- run contextual/run_add.sh"); return None
    comb  = read_fa(comb_p)
    ncal  = read_fa(os.path.join(PROJ, "nextclade", ty, "nextclade.aligned.fasta"))
    refun = read_fa(os.path.join(CTX, "reference", f"ref_RSV-{ty[-1]}.fasta"))
    genes = ref_genes(ty)
    refname = next(iter(refun))
    reflen  = len(next(iter(refun.values())))
    if refname not in comb:
        print(f"  [{ty}] reference {refname} not in combined alignment -- rerun run_add.sh"); return None

    melghat = sorted(k for k in comb if k in ncal)
    ctx     = sorted(k for k in comb if k not in ncal and k != refname)
    Lc      = len(comb[refname])
    print(f"\n[{ty}] combined {Lc} cols | {len(ctx)} contextual + {len(melghat)} Melghat + 1 reference")
    print(f"  reference {refname}: {reflen} nt (nextclade aligned length "
          f"{len(next(iter(ncal.values())))})")

    # ---- exact map, read off the reference row ------------------------------
    rc = base_cols(comb[refname])
    if len(rc) != reflen:
        print(f"  [{ty}] FATAL: reference lost bases in the combined alignment "
              f"({len(rc)} of {reflen})"); return None
    ref2comb = {i + 1: c + 1 for i, c in enumerate(rc)}       # 1-based both sides
    print(f"  mapped {len(ref2comb)}/{reflen} reference positions exactly (no voting)")

    def gene_of(p):
        hits = [g for g, (s, e) in genes.items() if s <= p <= e]
        return "|".join(sorted(hits)) if hits else "intergenic/UTR"

    with open(os.path.join(OUT, f"15_coordmap_{ty}.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["ref_pos", "combined_col", "gene", "ref_base"])
        refseq = next(iter(refun.values()))
        for p in range(1, reflen + 1):
            w.writerow([p, ref2comb[p], gene_of(p), refseq[p - 1]])

    for g in ("G", "F"):
        if g in genes:
            s, e = genes[g]
            print(f"  {g}: reference {s}-{e} ({e-s+1} nt, {(e-s+1)//3} codons) "
                  f"-> combined {ref2comb[s]}-{ref2comb[e]}")

    # ---- validate EVERY Melghat genome by base identity ---------------------
    ok, len_mismatch, base_mismatch = [], [], []
    for k in melghat:
        a, b = base_cols(ncal[k]), base_cols(comb[k])
        sa = "".join(ncal[k][i] for i in a)
        sb = "".join(comb[k][i] for i in b)
        if len(sa) != len(sb):
            len_mismatch.append((k, len(sa), len(sb))); continue
        if sa != sb:
            nd = sum(1 for x, y in zip(sa, sb) if x != y)
            base_mismatch.append((k, nd)); continue
        ok.append(k)
    print(f"  base-identity validation over all {len(melghat)} Melghat genomes: "
          f"{len(ok)} identical, {len(len_mismatch)} length mismatch, "
          f"{len(base_mismatch)} base mismatch")
    for k, la, lb in len_mismatch[:5]:
        print(f"    LENGTH {k}: nextclade {la} vs combined {lb}")
    for k, nd in base_mismatch[:5]:
        print(f"    BASES  {k}: {nd} differing positions")

    # ---- insertions Nextclade stripped, per sequence -------------------------
    ins = nextclade_insertions(ty)
    if not ins:
        print("  insertion accounting unavailable for this subgroup")
    tot = {k: ins.get(k, (0, []))[0] for k in melghat}
    n_with = sum(1 for v in tot.values() if v > 0)
    allins = sorted(sum(tot.values()) for _ in [0])[0]
    print(f"  Nextclade stripped insertions: {allins} nt across {n_with}/{len(melghat)} genomes"
          f" (max {max(tot.values()) if tot else 0} nt in one genome)")
    byloc = collections.Counter()
    for k in melghat:
        for pos, s in ins.get(k, (0, []))[1]:
            byloc[(gene_of(pos), len(s))] += 1
    with open(os.path.join(OUT, f"15_insertions_{ty}.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["seqName", "total_inserted_nt", "after_ref_pos", "gene", "inserted_len"])
        for k in melghat:
            t, lst = ins.get(k, (0, []))
            if not lst: w.writerow([k, t, "", "", ""])
            for pos, s in lst: w.writerow([k, t, pos, gene_of(pos), len(s)])
    print("  insertion hotspots (gene, length -> count):")
    for (g, L), n in byloc.most_common(6):
        print(f"    {g:28} {L:>4} nt   x{n}")

    # ---- contextual date precision ------------------------------------------
    prec = collections.Counter()
    for k in ctx:
        parts = k.split("_"); d = parts[2] if len(parts) > 2 else ""
        if   re.fullmatch(r"\d{4}-\d{2}-\d{2}", d): prec["full"] += 1
        elif re.fullmatch(r"\d{4}-\d{2}", d):       prec["year-month"] += 1
        elif re.fullmatch(r"\d{4}", d):             prec["year-only"] += 1
        else:                                       prec["unparsed"] += 1
    print(f"  contextual date precision: {dict(prec)}")

    return dict(subgroup=ty, reference=refname, ref_accession=REF_ACC[ty],
                ref_len=reflen, combined_cols=Lc, n_contextual=len(ctx),
                n_melghat=len(melghat), mapped_positions=len(ref2comb),
                validated_identical=len(ok), length_mismatch=len(len_mismatch),
                base_mismatch=len(base_mismatch),
                stripped_insertion_nt=allins, genomes_with_insertions=n_with,
                **{f"date_{k}": v for k, v in prec.items()})

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    rows = [r for r in (process(t) for t in ("RSV-A", "RSV-B")) if r]
    if rows:
        with open(os.path.join(OUT, "15_contextual_prep_summary.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=sorted({k for r in rows for k in r}))
            w.writeheader(); w.writerows(rows)
    # provenance
    prov = os.path.join(OUT, "15_provenance.txt")
    with open(prov, "w") as f:
        f.write("15_contextual_prep.py provenance\n")
        f.write(f"python {sys.version.split()[0]}\n")
        f.write("Nextclade dataset: nextstrain/rsv {a:EPI_ISL_412866, b:EPI_ISL_1653999}, "
                "reference tree updated 2025-08-23 (dataset builds 2025-08-25 and 2025-09-09 "
                "carry byte-identical annotations)\n")
        f.write("Contextual panel: Langedijk et al. 2024 Nat Commun 15:3083, "
                "github.com/bramvrancken/RSV_INFORM\n\n")
        for p in ("contextual/SHA256SUMS.txt",):
            fp = os.path.join(PROJ, p)
            if os.path.exists(fp): f.write(open(fp).read())
        for t in ("RSV-A", "RSV-B"):
            for rel in (f"output/15_coordmap_{t}.csv", f"output/15_insertions_{t}.csv"):
                fp = os.path.join(PROJ, rel)
                if os.path.exists(fp): f.write(f"{sha256(fp)}  {rel}\n")
    print(f"\n  -> wrote 15_coordmap_*.csv, 15_insertions_*.csv, "
          f"15_contextual_prep_summary.csv, 15_provenance.txt")
    print("== 15 done ==")

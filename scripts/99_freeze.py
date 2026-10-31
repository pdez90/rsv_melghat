#!/usr/bin/env python3
"""
99_freeze.py -- freeze the analysis as a DATED SNAPSHOT OF THE ACTUAL FILES plus
a manifest of real SHA-256 digests, and verify a snapshot on demand.

  python3 scripts/99_freeze.py            create freeze/<date>/ and its manifest
  python3 scripts/99_freeze.py --verify   re-hash every entry against the live
                                          files; exit 1 if any is missing/changed

Why this replaced the R version: a manifest RECORDS a state, it does not PRESERVE
one. The first version wrote md5 digests under a column named sha256, listed only
8 inputs (omitting the external cohort files, the Nextclade outputs and the
contextual alignments), and was run in one environment while the files it
described lived in another -- so five entries did not verify. This is Python so
it runs in any environment that can see the project, which is what let the
mismatch be caught at all.

Large files with a .gz sibling are snapshotted gzipped; the manifest records the
digest of the UNCOMPRESSED file, so gunzip reproduces a verifying copy.
"""
import csv, hashlib, os, shutil, subprocess, sys, datetime, platform

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
OUT  = os.path.join(PROJ, "output")
MANIFEST = "99_freeze_manifest.csv"

# mirrors DATA_DIRS in 00_setup.R
DATA_DIRS = [os.path.join(PROJ, "data"),
             os.path.join(PROJ, "..", "Final", "Data"),
             os.path.join(PROJ, "..", "Final", "LRTI data files"),
             os.path.join(PROJ, "..", "GISproject", "GISProject", "capstone project Mattie"),
             os.path.join(PROJ, "..")]
EXTERNAL = ["ari.csv", "demographics.csv", "lrtideaths.csv", "critical_data_gpsdata.csv",
            "93 village list.csv", "NP_results_with_clade_sequencing.csv",
            "village_characteristics.dta"]
AUTHORITATIVE = [
 "03_cohort_rsv_by_season.csv","03_lineage_by_season.csv","03_sequencing_representativeness.csv",
 "04_F_evidence_table.csv","12_O5_variant_screen.csv",
 "06_severity_model_sensitivity.csv","12_O2_severity.csv",
 "07_sample_disposition.csv","07_disposition_summary.csv","07_excluded_samples.csv",
 "12_O1_coverage.csv","12_O4_spatial.csv","13_comparison_table.csv","14_village_influence.csv",
 "15_contextual_prep_summary.csv","16_support_by_relevance.csv","16_melghat_clades.csv",
 "17_repeat_region.csv","19_mask_comparison.csv",
 "20_g_repeat_specimen_table.csv","20_priority_specimens_for_reads.csv",
 "22_spatial_conditioned.csv","23_lineage_by_setting.csv","23_rarefied_richness.csv",
 "23_diversity_by_setting.csv",
 "24_f_context_comparison.csv","24_f_combinations.csv","25_g_context_comparison.csv",
 "26_ascertainment_testing_step.csv","26_ascertainment_sequencing_step.csv","26_village_coverage.csv",
 "26_sequencing_step_within_season.csv",
 "27_discordant_cohort_status.csv","27_coverage_declared_definition.csv",
 "06_convergence_diagnosis.csv","28_dependency_check.csv",
 "25_g_boundary_agreement.csv",
 "29_household_key_diagnostic.csv","29_household_key_partition_identity.csv",
 "29_household_key_mother_field.csv"]

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    d = h.hexdigest()
    assert len(d) == 64, f"not a sha256 digest for {p}"
    return d

def data_path(name):
    for d in DATA_DIRS:
        p = os.path.join(d, name)
        if os.path.exists(p): return os.path.realpath(p)
    return None

def listdir(sub, exts=None, recursive=False):
    root = os.path.join(PROJ, sub)
    if not os.path.isdir(root): return []
    out = []
    if recursive:
        for dp, _, fs in os.walk(root):
            for f in fs: out.append(os.path.join(dp, f))
    else:
        out = [os.path.join(root, f) for f in os.listdir(root)
               if os.path.isfile(os.path.join(root, f))]
    if exts: out = [p for p in out if os.path.splitext(p)[1] in exts]
    return sorted(out)

def collect():
    items = []
    def add(paths, group):
        for p in paths:
            if p: items.append((group, p))
    add(listdir("scripts", {".R", ".py", ".sh"}), "script")
    add(listdir("data"), "input: project data")
    add([data_path(f) for f in EXTERNAL], "input: external cohort")
    add(listdir("nextclade", recursive=True), "input: nextclade")
    add(listdir(os.path.join("contextual", "reference")), "input: reference")
    add([p for p in listdir(os.path.join("contextual", "inform")) if p.endswith(".fasta")],
        "input: contextual panel")
    add([p for p in listdir("contextual") if os.path.basename(p).startswith("combined_RSV-")
         and p.endswith(".fasta")], "intermediate: alignment")
    add([p for p in listdir("contextual")
         if p.endswith(".nwk") or p.endswith(".sh")], "intermediate: tree and driver")
    # build-provenance sidecars written by run_tree.sh / run_masked.sh: these ARE
    # part of the reproducibility package, since they are what lets a later run
    # establish that a tree was built from the alignment and mask on disk. They
    # were previously collected by no group at all, so a recorded provenance would
    # not have been archived or verified.
    add([p for p in listdir("contextual") if p.endswith(".provenance")],
        "intermediate: tree build provenance")
    add([p for p in listdir("output") if os.path.basename(p).startswith("17_mask_")], "intermediate: mask")
    add([os.path.join(OUT, f) for f in AUTHORITATIVE], "authoritative result")
    add(listdir("manuscript", {".docx", ".pdf"}), "manuscript")
    # a file can be reached by two routes (an external cohort file copied into
    # data/, say); keep the first listing so each path appears exactly once
    seen, uniq = set(), []
    for g, q in items:
        key = os.path.realpath(q)
        if key in seen: continue
        seen.add(key); uniq.append((g, q))
    return uniq

def relpath(p):
    rp = os.path.relpath(os.path.realpath(p), PROJ)
    return os.path.join("external", os.path.basename(p)) if rp.startswith("..") else rp

def main():
    verify = "--verify" in sys.argv
    stamp  = datetime.date.today().isoformat()
    freeze = os.path.join(PROJ, "freeze", stamp)

    if verify:
        mp = os.path.join(freeze, MANIFEST)
        if not os.path.exists(mp): mp = os.path.join(OUT, MANIFEST)
        if not os.path.exists(mp): sys.exit("no manifest found; run the freeze first")
        rows = list(csv.DictReader(open(mp)))
        print(f"\n== 99 verify ==\n  manifest: {mp}")
        ok = changed = missing = skipped = 0
        # SNAPSHOT INTEGRITY is a SEPARATE check from live-file agreement, and
        # was previously not performed at all: the loop below only re-hashed the
        # LIVE files, so a corrupted or truncated snapshot copy would pass. Each
        # archived copy is now hashed too, gzipped copies after decompression,
        # since the manifest records the digest of the UNCOMPRESSED bytes.
        snap_ok = snap_bad = snap_absent = 0
        def snap_digest(rec):
            # the manifest's `snapshot` field ALREADY carries the .gz suffix when
            # the copy is compressed -- appending another one looks for a file
            # that never existed and reports every compressed entry as absent.
            sp = os.path.join(freeze, rec["snapshot"])
            gz = str(rec.get("snapshot_gzipped", "")).strip().lower() == "true"
            if not os.path.exists(sp): return None
            if not gz: return sha256(sp)
            import gzip
            h = hashlib.sha256()
            with gzip.open(sp, "rb") as fh:
                for b in iter(lambda: fh.read(1 << 20), b""): h.update(b)
            return h.hexdigest()
        for r in rows:
            if r["path"].startswith("external/"):
                # verify these too: by the absolute path recorded at freeze time,
                # or failing that against the snapshot copy, so nothing is
                # silently unchecked
                p = r.get("abs_path") or ""
                if not p or not os.path.exists(p):
                    p = os.path.join(freeze, r["snapshot"])
                if not os.path.exists(p):
                    print(f"  UNCHECKED {r['path']} (source gone, no snapshot copy)")
                    skipped += 1; continue
                if sha256(p) != r["sha256"]:
                    print(f"  CHANGED   {r['path']}"); changed += 1
                else: ok += 1
                continue
            p = os.path.join(PROJ, r["path"])
            if not os.path.exists(p): print(f"  MISSING   {r['path']}"); missing += 1
            elif sha256(p) != r["sha256"]: print(f"  CHANGED   {r['path']}"); changed += 1
            else: ok += 1
        for r in rows:
            d = snap_digest(r)
            if d is None:
                print(f"  SNAPSHOT ABSENT  {r['snapshot']}"); snap_absent += 1
            elif d != r["sha256"]:
                print(f"  SNAPSHOT CORRUPT {r['snapshot']}"); snap_bad += 1
            else: snap_ok += 1
        print(f"  live files:  {ok} ok, {changed} changed, {missing} missing, "
              f"{skipped} unchecked")
        print(f"  snapshots:   {snap_ok} ok, {snap_bad} corrupt, {snap_absent} absent")
        # An UNCHECKED entry is a FAILURE, not a pass with a note: it means an
        # external input is absent from its recorded location AND absent from the
        # snapshot, so the manifest no longer evidences anything about that file.
        # Treating it as success let a verification pass while part of the
        # analysis inputs could not be confirmed at all.
        if changed or missing or skipped or snap_bad or snap_absent:
            print("== VERIFICATION FAILED ==")
            if snap_bad or snap_absent:
                print("  reason: the archived snapshot copy of one or more entries "
                      "is absent or does not match its recorded digest, so the "
                      "snapshot is not a faithful archive even where the live "
                      "files agree. These are separate guarantees.")
            if skipped and not (changed or missing):
                print("  reason: unchecked entries only -- an external input is "
                      "missing from both its recorded path and the snapshot, so "
                      "it cannot be verified. Restore the file or re-freeze.")
            sys.exit(1)
        print("== verification passed ==")
        return

    print("\n== 99 freeze ==")
    items = collect()
    gone = [p for g, p in items if not os.path.exists(p)]
    if gone:
        print("  MISSING, nothing frozen:")
        for p in gone: print("    -", relpath(p))
        sys.exit(1)
    os.makedirs(freeze, exist_ok=True)
    man, nbytes = [], 0
    for group, p in items:
        rel = relpath(p)
        gz  = p + ".gz"
        use_gz = os.path.exists(gz) and os.path.getsize(p) > 5_000_000
        src  = gz if use_gz else p
        dest = os.path.join(freeze, rel + (".gz" if use_gz else ""))
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copy2(src, dest); nbytes += os.path.getsize(src)
        man.append(dict(group=group, path=rel, bytes=os.path.getsize(p), sha256=sha256(p),
                        snapshot=os.path.relpath(dest, freeze), snapshot_gzipped=use_gz,
                        abs_path=os.path.realpath(p) if rel.startswith("external/") else ""))
    for target in (os.path.join(freeze, MANIFEST), os.path.join(OUT, MANIFEST)):
        with open(target, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(man[0].keys())); w.writeheader(); w.writerows(man)
    info = os.path.join(freeze, "99_freeze_sessionInfo.txt")
    with open(info, "w") as f:
        f.write("Melghat RSV genomics -- frozen analysis\n")
        f.write(f"frozen: {datetime.datetime.now().isoformat(timespec='seconds')} on {platform.node()}\n")
        f.write("seed: 20261030 throughout\nsnapshot: " + freeze + "\n\n")
        f.write("Nextclade dataset: nextstrain/rsv a/EPI_ISL_412866, b/EPI_ISL_1653999;\n")
        f.write("  reference tree updated 2025-08-23\n")
        f.write("Contextual panel: Langedijk et al. 2024 Nat Commun 15:3083;\n")
        f.write("  github.com/bramvrancken/RSV_INFORM (ends April 2020)\n\n")
        f.write(f"python {platform.python_version()} on {platform.platform()}\n")
        for tool in ("R", "mafft", "fasttree"):
            try:
                v = subprocess.run([tool, "--version"], capture_output=True, text=True, timeout=10)
                f.write(f"{tool}: {(v.stdout or v.stderr).splitlines()[0] if (v.stdout or v.stderr) else 'unknown'}\n")
            except Exception: f.write(f"{tool}: not found in this environment\n")
    shutil.copy2(info, os.path.join(OUT, "99_freeze_sessionInfo.txt"))
    counts = {}
    for m in man: counts[m["group"]] = counts.get(m["group"], 0) + 1
    for g in sorted(counts): print(f"  {g:28} {counts[g]:>4} files")
    print(f"\n  {len(man)} files snapshotted ({nbytes/1e6:.0f} MB) -> {freeze}")
    print("  verify with:  python3 scripts/99_freeze.py --verify\n== 99 done ==")

if __name__ == "__main__":
    main()

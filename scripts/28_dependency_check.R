# 28_dependency_check.R -- do the derived contextual artefacts actually
# correspond to the inputs currently on disk?
#
# WHY: run_all.R previously checked only that prerequisite FILES EXIST. Existence
# is not currency: a tree left over from an earlier alignment, or a masked tree
# built from a superseded mask, satisfies an existence test while invalidating
# every downstream analysis that reads it. The driver's claim to refuse stale
# inputs was therefore inaccurate.
#
# Timestamps are not a sound test either. In this project contextual/inform/*.fasta
# is 4 seconds NEWER than the combined alignment built from it, because run_add.sh
# touches the decompressed input during the same build; a naive mtime rule calls
# that stale when it is not. So mtimes are reported as ADVISORY only.
#
# What this script checks, stated at its real strength:
#
#  1. TAXON-SET CONSISTENCY. The set of tip labels in each tree equals the set of
#     sequence names in the alignment it is supposed to describe. This is NOT a
#     provenance test: two alignments can carry identical sequence names but
#     different nucleotides, different columns or a different mask, and both would
#     pass. An earlier version of this header claimed such a check "cannot pass for
#     a tree built from a different alignment" -- that was false and is withdrawn.
#     Failing it proves the tree does not belong to this alignment; passing it
#     leaves build provenance open.
#
#  2. BUILD PROVENANCE, where it was recorded. run_tree.sh and run_masked.sh now
#     write <tree>.provenance beside each tree, holding the SHA-256 of the exact
#     alignment (and, for masked trees, the mask) used. Where that file exists the
#     digest is re-computed from the inputs now on disk and compared, which does
#     establish correspondence. Where it does not exist -- any tree built before
#     provenance recording was added -- this is reported as UNVERIFIED rather than
#     passed, and the only ways to close it are to rebuild from the current inputs
#     or to accept the taxon-set check as the weaker guarantee it is.
#
# Neither check verifies branch lengths or substitution model settings.
source(file.path("scripts", "00_setup.R"))
cat("\n== 28 dependency check (derived artefacts vs current inputs) ==\n")

digest_sha256 <- function(p) {
  # no hard dependency on a hashing package: shell out to the platform tool
  out <- tryCatch(system2(if (nzchar(Sys.which("sha256sum"))) "sha256sum" else "shasum",
                          c(if (nzchar(Sys.which("sha256sum"))) character(0) else c("-a","256"),
                            shQuote(p)), stdout = TRUE, stderr = FALSE),
                  error = function(e) NA_character_)
  if (length(out) == 0 || is.na(out[1])) return(NA_character_)
  sub("\\s.*$", "", out[1])
}

fasta_names <- function(p) {
  h <- grep("^>", readLines(p, warn = FALSE), value = TRUE)
  sub("^>", "", sub("[ \t].*$", "", h))
}

rows <- list(); fail <- character(0)
for (ty in c("RSV-A", "RSV-B")) {
  s   <- substr(ty, 5, 5)
  aln <- file.path(PROJ, "contextual", sprintf("combined_%s.fasta", ty))
  if (!file.exists(aln)) { cat(sprintf("  %s: alignment absent, skipped\n", ty)); next }
  an  <- fasta_names(aln)
  for (lbl in c("tree", "tree_masked")) {
    tf <- file.path(PROJ, "contextual", sprintf("%s_%s.nwk", lbl, ty))
    if (!file.exists(tf)) { cat(sprintf("  %s %s: absent\n", ty, lbl)); next }
    tips <- ape::read.tree(tf)$tip.label
    ## build provenance, if it was recorded at build time
    pf <- paste0(tf, ".provenance")
    prov <- "UNVERIFIED: no provenance file (tree predates provenance recording)"
    if (file.exists(pf)) {
      kv <- strsplit(readLines(pf, warn = FALSE), "=", fixed = TRUE)
      kv <- setNames(vapply(kv, function(x) paste(x[-1], collapse = "="), ""),
                     vapply(kv, `[`, "", 1))
      want_aln <- unname(kv["alignment_sha256"])
      got_aln  <- digest_sha256(aln)
      ok_aln   <- !is.na(want_aln) && identical(want_aln, got_aln)
      ok_mask  <- TRUE
      if ("mask_sha256" %in% names(kv)) {
        mf <- file.path(PROJ, unname(kv["mask_file"]))
        ok_mask <- file.exists(mf) && identical(unname(kv["mask_sha256"]),
                                                digest_sha256(mf))
      }
      prov <- if (ok_aln && ok_mask) "VERIFIED against recorded input digests"
              else paste0("MISMATCH: recorded inputs differ from those on disk",
                          if (!ok_aln) " (alignment)" else "",
                          if (!ok_mask) " (mask)" else "")
    }
    only_tree <- setdiff(tips, an); only_aln <- setdiff(an, tips)
    agree <- length(only_tree) == 0 && length(only_aln) == 0
    rows[[paste(ty, lbl)]] <- tibble::tibble(
      subgroup = ty, artefact = basename(tf), alignment = basename(aln),
      n_alignment = length(an), n_tree_tips = length(tips),
      only_in_tree = length(only_tree), only_in_alignment = length(only_aln),
      labels_agree = agree, build_provenance = prov,
      tree_mtime = format(file.mtime(tf)), alignment_mtime = format(file.mtime(aln)),
      tree_newer_than_alignment = file.mtime(tf) >= file.mtime(aln))
    cat(sprintf("  %-22s tips %5d vs alignment %5d  %s\n", basename(tf),
                length(tips), length(an),
                if (agree) "labels agree" else
                  sprintf("MISMATCH: %d tree-only, %d alignment-only",
                          length(only_tree), length(only_aln))))
    cat(sprintf("                         provenance: %s\n", prov))
    if (!agree) fail <- c(fail, basename(tf))
  }
}

if (length(rows)) {
  out <- dplyr::bind_rows(rows)
  write_out(out, "28_dependency_check.csv")
  cat("\n  labels_agree is a TAXON-SET check: failing it disproves correspondence,\n")
  cat("  passing it does not establish that the tree was built from these inputs.\n")
  cat("  build_provenance is the stronger test and applies only where the builder\n")
  cat("  recorded input digests. mtime ordering is advisory -- see the header.\n")
  if (any(grepl("^UNVERIFIED", out$build_provenance)))
    cat("  NOTE: at least one tree has no recorded provenance. Rebuild with\n",
        "       contextual/run_tree.sh and run_masked.sh to establish it.\n", sep = "")
} else {
  write_out(tibble::tibble(note = "no contextual artefacts present to check"),
            "28_dependency_check.csv")
}

if (length(fail)) {
  stop(sprintf(paste("tip labels do not match the current alignment for: %s.",
                     "These trees were built from different input; rebuild with",
                     "contextual/run_tree.sh and contextual/run_masked.sh before",
                     "running any tree-dependent analysis."),
               paste(fail, collapse = ", ")))
}
cat("== 28 done ==\n")

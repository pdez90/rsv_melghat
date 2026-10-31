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
# The decisive check is CONTENT CORRESPONDENCE: the set of tip labels in each tree
# must equal the set of sequence names in the alignment it is supposed to describe.
# That cannot pass for a tree built from a different alignment, and it needs no
# build-time bookkeeping. It does NOT verify branch lengths or that the same
# masking was applied -- for the masked trees it establishes the same taxon set,
# and the mask itself is checked by 17_repeat_region.py's own assertions.
source(file.path("scripts", "00_setup.R"))
cat("\n== 28 dependency check (derived artefacts vs current inputs) ==\n")

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
    only_tree <- setdiff(tips, an); only_aln <- setdiff(an, tips)
    agree <- length(only_tree) == 0 && length(only_aln) == 0
    rows[[paste(ty, lbl)]] <- tibble::tibble(
      subgroup = ty, artefact = basename(tf), alignment = basename(aln),
      n_alignment = length(an), n_tree_tips = length(tips),
      only_in_tree = length(only_tree), only_in_alignment = length(only_aln),
      labels_agree = agree,
      tree_mtime = format(file.mtime(tf)), alignment_mtime = format(file.mtime(aln)),
      tree_newer_than_alignment = file.mtime(tf) >= file.mtime(aln))
    cat(sprintf("  %-22s tips %5d vs alignment %5d  %s\n", basename(tf),
                length(tips), length(an),
                if (agree) "labels agree" else
                  sprintf("MISMATCH: %d tree-only, %d alignment-only",
                          length(only_tree), length(only_aln))))
    if (!agree) fail <- c(fail, basename(tf))
  }
}

if (length(rows)) {
  out <- dplyr::bind_rows(rows)
  write_out(out, "28_dependency_check.csv")
  cat("\n  mtime ordering is ADVISORY only -- see the header. The decisive\n")
  cat("  column is labels_agree.\n")
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

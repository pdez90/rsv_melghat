# run_all.R -- reproduce the Melghat RSV genomics analysis.
# Run from the project folder:  Rscript scripts/run_all.R
#
#   --core          stages A and C only (skips the slow resampling panel)
#   --allow-skips   exit 0 even if a stage was skipped for missing prerequisites
#                   (WITHOUT this flag, any skipped stage makes the run exit 1,
#                    so a partial run can never be mistaken for a full one)
#
# Requires R: readr dplyr tidyr stringr purrr ggplot2 jsonlite ape haven lme4
#             broom.mixed
# Requires python3 (standard library only) for stages D, E and F.
# Nextclade output must be present in nextclade/RSV-A and nextclade/RSV-B.
#
# WHY STAGES: an earlier version of this driver stopped after step 14. Scripts
# 15-27 were run by hand, so a "successful" run followed by a freeze could
# preserve stale results from the later analyses -- the conditioned spatial test,
# the community-facility comparison, the ascertainment tables, the F/G contextual
# comparisons and the cohort-status audit. Every script that produces a result
# quoted in the manuscript is now driven from here.
#
# NOT driven from here, by design -- they need external binaries and rebuild the
# contextual alignments and trees themselves:
#   contextual/run_add.sh     combined alignments  (needs mafft)
#   contextual/run_tree.sh    exploratory trees    (needs FastTree)
#   contextual/run_masked.sh  masked trees         (needs FastTree)
#   scripts/21_prepare_context.py  only when ADDING new contextual genomes
# Run those three when the alignments or trees must be regenerated; the delivered
# copies carry checksums in contextual/SHA256SUMS.txt. Stages D and E below check
# for their outputs and refuse to run on absent or stale inputs rather than
# silently producing nothing.

args        <- commandArgs(trailingOnly = TRUE)
core        <- "--core"        %in% args
allow_skips <- "--allow-skips" %in% args

PROJ <- normalizePath(".")
OUT  <- file.path(PROJ, "output")
skipped <- character(0)
ran     <- character(0)

have <- function(...) all(file.exists(file.path(PROJ, c(...))))

run_r <- function(s) {
  cat("\n", strrep("=", 70), "\n", s, "\n", strrep("=", 70), "\n", sep = "")
  source(file.path("scripts", s), echo = FALSE)
  ran <<- c(ran, s)
}

run_cmd <- function(label, cmd, cmd_args) {
  cat("\n", strrep("=", 70), "\n", label, "\n", strrep("=", 70), "\n", sep = "")
  st <- system2(cmd, cmd_args)
  if (!identical(st, 0L)) stop(sprintf("%s exited with status %s", label, st))
  ran <<- c(ran, label)
}

stage <- function(name, needs, body) {
  missing <- needs[!file.exists(file.path(PROJ, needs))]
  if (length(missing)) {
    cat(sprintf("\n-- SKIPPED stage %s: missing %s\n", name,
                paste(missing, collapse = ", ")))
    skipped <<- c(skipped, sprintf("%s (missing %s)", name,
                                   paste(missing, collapse = ", ")))
    return(invisible(FALSE))
  }
  body(); invisible(TRUE)
}

t0 <- Sys.time()

## ---- stage A: cohort linkage, typing, primary analyses ---------------------
for (s in c("00_verify_data.R", "01_link_sequencing.R", "01b_export_fasta.R",
            "02_nextclade_results.R", "03_lineage_turnover.R",
            "04_fgene_antigenic.R", "05_phylogeography.R",
            "06_severity_reinfection.R", "07_data_gaps.R", "08_robustness.R",
            "09_ggene.R", "10_validation.R", "11_provenance.R")) run_r(s)

## ---- stage B: the slow resampling panel ------------------------------------
if (core) {
  cat("\n-- SKIPPED stage B (sensitivity panel): --core requested\n")
  skipped <- c(skipped, "B sensitivity panel (--core requested)")
} else {
  for (s in c("12_sensitivity.R", "13_sensitivity_summary.R",
              "14_village_influence.R")) run_r(s)
}

## ---- stage C: analyses needing only stage A output -------------------------
## These were previously run by hand and are the ones most at risk of going
## stale: every figure they produce is quoted in the manuscript.
for (s in c("22_spatial_conditioned.R", "23_community_vs_facility.R",
            "26_ascertainment.R", "27_cohort_status_audit.R")) run_r(s)

## ---- stage D: contextual alignment harmonisation and masks -----------------
## Must precede stage E: 15_contextual_prep.py writes output/15_coordmap_*.csv,
## which 20_gapfrac.py reads, and 17_repeat_region.py maps every mask window
## through that same coordinate map.
stage("D contextual prep + masks",
      c("scripts/15_contextual_prep.py", "contextual/combined_RSV-A.fasta",
        "contextual/combined_RSV-B.fasta"),
      function() {
        run_cmd("15_contextual_prep.py", "python3", c("scripts/15_contextual_prep.py"))
        run_cmd("17_repeat_region.py",   "python3", c("scripts/17_repeat_region.py"))
      })

## ---- stage E: G-repeat specimen table (needs stage D's coordinate map) -----
stage("E G-repeat table",
      c("scripts/20_gapfrac.py", "output/15_coordmap_RSV-A.csv"),
      function() {
        run_cmd("20_gapfrac.py", "python3", c("scripts/20_gapfrac.py"))
        run_r("20_g_repeat_specimen_table.R")
      })

stage("E tree-dependent analyses",
      c("contextual/tree_RSV-A.nwk", "contextual/tree_RSV-B.nwk",
        "contextual/tree_masked_RSV-A.nwk", "contextual/tree_masked_RSV-B.nwk"),
      function() {
        for (s in c("16_tree_annotation.R", "18_annotated_tree.R",
                    "19_mask_comparison.R")) run_r(s)
      })

## ---- stage F: F and G states against the contextual panel ------------------
stage("F contextual comparisons",
      c("scripts/24_f_context_comparison.py", "scripts/25_g_context_comparison.py",
        "contextual/combined_RSV-A.fasta"),
      function() {
        run_cmd("24_f_context_comparison.py", "python3",
                c("scripts/24_f_context_comparison.py"))
        run_cmd("25_g_context_comparison.py", "python3",
                c("scripts/25_g_context_comparison.py"))
      })

## ---- report ----------------------------------------------------------------
cat(sprintf("\n%s\nRan %d steps in %.1f min. Results are in output/.\n",
            strrep("=", 70), length(ran),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
if (length(skipped)) {
  cat("\nSKIPPED STAGES -- the analysis is INCOMPLETE and any freeze taken now\n")
  cat("will carry stale results for these:\n")
  for (k in skipped) cat("  - ", k, "\n", sep = "")
  if (!allow_skips) {
    cat("\nExiting non-zero so this does not read as a full run.\n")
    cat("Pass --allow-skips if a partial run is what you intended.\n")
    quit(status = 1)
  }
} else {
  cat("\nAll stages ran. Freeze with:  python3 scripts/99_freeze.py\n")
}

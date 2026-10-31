# run_all.R -- reproduce the whole Melghat RSV genomics analysis.
# Run from the project folder:  Rscript scripts/run_all.R
#
# Requires: readr dplyr tidyr stringr purrr ggplot2 jsonlite ape haven lme4
#           broom.mixed
# Nextclade output must already be present in nextclade/RSV-A and nextclade/RSV-B.
#
# Steps 12-14 are the sensitivity panel and the village-influence analysis.
# They are slower than the rest (permutation and bootstrap resampling; a few
# minutes). Skip them with:  Rscript scripts/run_all.R --core
#
# NOT run from here:
#   scripts/15_contextual_prep.py   contextual alignment prep -- needs python3 and
#                                   contextual/combined_RSV-*.fasta
#   contextual/run_add.sh           rebuilds those alignments -- needs mafft
#   contextual/run_tree.sh          exploratory trees -- needs FastTree
# Those three are only needed if the alignments themselves must be regenerated;
# the delivered copies carry checksums in contextual/SHA256SUMS.txt.

args  <- commandArgs(trailingOnly = TRUE)
core  <- "--core" %in% args

core_steps <- c("00_verify_data.R", "01_link_sequencing.R", "01b_export_fasta.R",
                "02_nextclade_results.R", "03_lineage_turnover.R",
                "04_fgene_antigenic.R", "05_phylogeography.R",
                "06_severity_reinfection.R", "07_data_gaps.R", "08_robustness.R",
                "09_ggene.R", "10_validation.R", "11_provenance.R")

slow_steps <- c("12_sensitivity.R", "13_sensitivity_summary.R",
                "14_village_influence.R")

steps <- if (core) core_steps else c(core_steps, slow_steps)

t0 <- Sys.time()
for (s in steps) {
  cat("\n", strrep("=", 70), "\n", s, "\n", strrep("=", 70), "\n", sep = "")
  source(file.path("scripts", s), echo = FALSE)
}
cat(sprintf("\nAll steps complete in %.1f min. Results are in output/.\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
if (core) cat("Skipped 12-14 (--core). Drop the flag to run them.\n")

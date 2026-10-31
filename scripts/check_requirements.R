# check_requirements.R -- run this FIRST. It reports anything missing and
# prints the one command that installs it. Changes nothing on your machine.
#
#   Rscript scripts/check_requirements.R

cat("\n=============== MELGHAT RSV PIPELINE: REQUIREMENTS CHECK ===============\n\n")
ok_all <- TRUE
say <- function(status, what, detail = "") {
  cat(sprintf("  [%-4s] %-28s %s\n", status, what, detail))
  if (status == "MISS") ok_all <<- FALSE
}

## ---- R version ---------------------------------------------------------------
cat("R\n")
rv <- getRversion()
say(if (rv >= "4.1.0") "OK" else "MISS", R.version.string,
    if (rv >= "4.1.0") "" else "need R >= 4.1.0")

## ---- packages ----------------------------------------------------------------
cat("\nR packages\n")
REQ <- c(readr = "1.4.0", dplyr = "1.1.0", tidyr = "1.1.0", stringr = "1.4.0",
         purrr = "0.3.0", ggplot2 = "3.3.0", jsonlite = "1.7.0", ape = "5.5",
         haven = "2.4.0", lme4 = "1.1-27", broom.mixed = "0.2.9")
missing <- character(0); tooold <- character(0)
for (p in names(REQ)) {
  if (!requireNamespace(p, quietly = TRUE)) { say("MISS", p, "not installed"); missing <- c(missing, p) }
  else {
    v <- as.character(utils::packageVersion(p))
    if (utils::compareVersion(v, REQ[[p]]) < 0) {
      say("OLD", p, sprintf("%s installed, need >= %s", v, REQ[[p]])); tooold <- c(tooold, p)
    } else say("OK", p, v)
  }
}
# dplyr >= 1.1.0 is a hard requirement: the linkage uses relationship = "many-to-many"
if (requireNamespace("dplyr", quietly = TRUE) &&
    utils::compareVersion(as.character(utils::packageVersion("dplyr")), "1.1.0") < 0)
  cat("\n  !! dplyr must be >= 1.1.0 -- 01_link_sequencing.R uses relationship=\"many-to-many\"\n")

## ---- input data --------------------------------------------------------------
cat("\nInput data\n")
PROJ <- getwd()
DIRS <- c(file.path(PROJ, "data"),
          file.path(PROJ, "..", "Final", "Data"),
          file.path(PROJ, "..", "GISproject", "GISProject", "capstone project Mattie"),
          file.path(PROJ, ".."))
find1 <- function(f) { h <- file.path(DIRS, f); h <- h[file.exists(h)]
                       if (length(h)) normalizePath(h[1]) else NA_character_ }
NEEDED <- c("NP_results_with_clade_sequencing.csv", "ari.csv", "lrtideaths.csv",
            "demographics.csv", "93 village list.csv", "critical_data_gpsdata.csv",
            "village_characteristics.dta")
for (f in NEEDED) {
  p <- find1(f)
  if (is.na(p)) say("MISS", f, "not found in data/, ../Final/Data/, ../GISproject/..., ..")
  else say("OK", f, sub(paste0(normalizePath(file.path(PROJ, "..")), "/"), "../", p, fixed = TRUE))
}

## ---- Nextclade output --------------------------------------------------------
cat("\nNextclade output (produced in the browser, nothing to install)\n")
for (ty in c("RSV-A", "RSV-B")) {
  need <- c("nextclade.tsv", "nextclade.aligned.fasta",
            "nextclade.cds_translation.F.fasta", "nextclade.cds_translation.G.fasta")
  have <- file.exists(file.path(PROJ, "nextclade", ty, need))
  say(if (all(have)) "OK" else "MISS", paste("nextclade/", ty, sep = ""),
      sprintf("%d of %d files", sum(have), length(need)))
}

## ---- optional ----------------------------------------------------------------
cat("\nOptional (only for 02b_clade_assign_crosscheck.R, not run by run_all.R)\n")
m <- suppressWarnings(nchar(Sys.which("mafft")))
say(if (m > 0) "OK" else "note", "mafft", if (m > 0) Sys.which("mafft") else "absent; the cross-check script will not run")

## ---- verdict -----------------------------------------------------------------
cat("\n-----------------------------------------------------------------------\n")
if (length(missing) || length(tooold)) {
  cat("Install what is missing with:\n\n")
  cat(sprintf("  install.packages(c(%s))\n\n",
              paste(sprintf('"%s"', unique(c(missing, tooold))), collapse = ", ")))
}
if (ok_all) {
  cat("All requirements satisfied. Run:  Rscript scripts/run_all.R\n")
} else {
  cat("Resolve the MISS items above, then re-run this check.\n")
}
cat("\n")

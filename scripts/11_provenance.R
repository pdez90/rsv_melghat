# 11_provenance.R -- record exactly what produced the results.
#
# Captures: checksums of every input and script, R and package versions, and
# whatever Nextclade metadata is recoverable from its own output. Anything not
# recoverable is written as "UNKNOWN" rather than omitted, so the gap is visible.
#
# Run this LAST. It writes output/11_provenance.txt and 11_provenance_files.csv.

source(file.path("scripts", "00_setup.R"))
suppressPackageStartupMessages(library(tools))
cat("\n== 11 provenance ==\n")

con <- file(file.path(OUT, "11_provenance.txt"), open = "wt")
say <- function(...) { line <- sprintf(...); cat(line, "\n"); writeLines(line, con) }

say("MELGHAT RSV GENOMICS -- PROVENANCE RECORD")
say("generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
say("project:   %s", normalizePath(PROJ))
say("")

## ---- inputs -----------------------------------------------------------------
say("================= INPUT FILES =================")
INPUTS <- c("NP_results_with_clade_sequencing.csv", "nextclade_input_RSV-A.fasta",
            "nextclade_input_RSV-B.fasta", "melghat_rsv_genomes.fasta",
            "ari.csv", "lrtideaths.csv", "demographics.csv",
            "93 village list.csv", "critical_data_gpsdata.csv",
            "village_characteristics.dta")
rows <- list()
for (f in INPUTS) {
  pth <- data_path(f, required = FALSE)
  if (is.na(pth)) { say("%-42s MISSING", f); next }
  md5 <- unname(md5sum(pth))
  say("%-42s %10d B  %s", f, file.size(pth), md5)
  rows[[f]] <- tibble(kind = "input", file = f, path = pth,
                      bytes = file.size(pth), md5 = md5,
                      mtime = format(file.mtime(pth), "%Y-%m-%d %H:%M:%S"))
}

## ---- Nextclade --------------------------------------------------------------
say("")
say("================= NEXTCLADE =================")
say("Nextclade was run in the browser (clades.nextstrain.org). run_all.R does not")
say("regenerate this stage, so the output below is taken on trust from the files.")
for (ty in c("RSV-A", "RSV-B")) {
  aus <- file.path(PROJ, "nextclade", ty, "nextclade.auspice.json")
  say("")
  say("-- %s --", ty)
  if (file.exists(aus)) {
    j <- fromJSON(aus, simplifyVector = FALSE)
    say("  reference tree title : %s", j$meta$title %||% "UNKNOWN")
    say("  tree built           : %s  (metadata only -- this is the reference TREE",
        j$meta$updated %||% "UNKNOWN")
    say("                           build date, NOT an identifier for the Nextclade")
    say("                           dataset release. The dataset tag stays UNKNOWN.)")
    say("  data provenance      : %s",
        paste(vapply(j$meta$data_provenance %||% list(), function(x) x$name %||% "?", character(1)),
              collapse = ", "))
    say("  auspice schema       : %s", j$version %||% "UNKNOWN")
  } else say("  auspice json         : MISSING")
  for (f in c("nextclade.tsv", "nextclade.aligned.fasta", "nextclade.auspice.json",
              "nextclade.cds_translation.F.fasta", "nextclade.cds_translation.G.fasta")) {
    pth <- file.path(PROJ, "nextclade", ty, f)
    if (file.exists(pth)) {
      md5 <- unname(md5sum(pth))
      say("  %-34s %10d B  %s", f, file.size(pth), md5)
      rows[[paste(ty, f)]] <- tibble(kind = paste("nextclade", ty), file = f, path = pth,
                                     bytes = file.size(pth), md5 = md5,
                                     mtime = format(file.mtime(pth), "%Y-%m-%d %H:%M:%S"))
    } else say("  %-34s MISSING", f)
  }
}
say("")
say("  Nextclade application version : UNKNOWN (v3.28.0 was reported in the browser")
say("    at the time of the run, but is not recorded in any output file)")
say("  dataset tag / release         : UNKNOWN. Not emitted by Nextclade Web output.")
say("    The reference accessions are RSV-A EPI_ISL_412866 and RSV-B EPI_ISL_1653999,")
say("    and the tree build date is recorded above, but neither identifies the dataset")
save_note <- "    release; two releases can share a reference and a tree build date."
say("%s", save_note)
say("  => To close these gaps, re-run via the Nextclade CLI, which records both.")

## ---- scripts ----------------------------------------------------------------
say("")
say("================= SCRIPTS =================")
for (f in sort(list.files(file.path(PROJ, "scripts"), pattern = "[.]R$"))) {
  pth <- file.path(PROJ, "scripts", f)
  md5 <- unname(md5sum(pth))
  say("%-36s %8d B  %s", f, file.size(pth), md5)
  rows[[paste("script", f)]] <- tibble(kind = "script", file = f, path = pth,
                                       bytes = file.size(pth), md5 = md5,
                                       mtime = format(file.mtime(pth), "%Y-%m-%d %H:%M:%S"))
}

## ---- planning and reporting documents ---------------------------------------
say("")
say("================= DOCUMENTS =================")
for (f in c("SENSITIVITY_PLAN.md", "README.md", "Melghat_RSV_report.docx")) {
  pth <- file.path(PROJ, f)
  if (file.exists(pth)) {
    md5 <- unname(md5sum(pth))
    say("%-36s %8d B  %s", f, file.size(pth), md5)
    rows[[paste("doc", f)]] <- tibble(kind = "document", file = f, path = pth,
                                      bytes = file.size(pth), md5 = md5,
                                      mtime = format(file.mtime(pth), "%Y-%m-%d %H:%M:%S"))
  } else say("%-36s MISSING", f)
}

## ---- software ---------------------------------------------------------------
say("")
say("================= SOFTWARE =================")
say("%s", R.version.string)
say("platform: %s", R.version$platform)
say("RNG kind: %s / %s / %s", RNGkind()[1], RNGkind()[2], RNGkind()[3])
say("seed used in scripts 05, 08, 09, 10: 20261030")
say("")
PKGS <- c("readr", "dplyr", "tidyr", "stringr", "purrr", "ggplot2", "jsonlite",
          "ape", "haven", "lme4", "broom.mixed")
for (p in PKGS)
  say("  %-14s %s", p,
      if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "NOT INSTALLED")

say("")
say("================= FULL sessionInfo() =================")
si <- capture.output(utils::sessionInfo())
for (l in si) say("%s", l)

say("")
say("================= NOT CAPTURED =================")
say("- Read-level sequencing data (depth, base quality, strand balance). Only")
say("  consensus genomes are available, so variant calls cannot be checked at")
say("  read level and the G repeat length cannot be resolved.")
say("- Laboratory metadata: extraction batch, sequencing run, assembly pipeline")
say("  and version used to produce the consensus genomes.")
say("- Nextclade application version and dataset tag (see above).")

close(con)
write_out(bind_rows(rows), "11_provenance_files.csv")
cat(sprintf("  -> wrote 11_provenance.txt\n"))
cat("== 11 done ==\n")

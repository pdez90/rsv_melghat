# 00_setup.R -- shared paths, libraries, helpers for the Melghat RSV genomics pipeline
# Priyanka deSouza / MAHAN Melghat cohort.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(ggplot2); library(jsonlite); library(ape); library(haven)
})

## Project root: run these scripts from the project folder (the one holding
## data/, scripts/ and output/), e.g.  Rscript scripts/01_link_sequencing.R
PROJ <- Sys.getenv("MELGHAT_PROJ", unset = "")
if (!nzchar(PROJ)) {
  cand <- c(getwd(), dirname(getwd()))
  hit  <- cand[dir.exists(file.path(cand, "data"))]
  PROJ <- if (length(hit)) hit[1] else getwd()
}

DATA <- file.path(PROJ, "data")

## Input files live in several folders of the RSV_India tree, not one:
##   data/                      sequencing inputs kept with this analysis
##   ../Final/Data/             cohort data (ari, lrtideaths, demographics)
##   ../GISproject/...          household GPS
##   ../                        village_characteristics.dta
## data_path() finds a file by name across those, so the scripts run whether
## the inputs are spread across the tree or copied into one folder.
DATA_DIRS <- c(
  DATA,
  file.path(PROJ, "..", "Final", "Data"),
  file.path(PROJ, "..", "Final", "LRTI data files"),
  file.path(PROJ, "..", "GISproject", "GISProject", "capstone project Mattie"),
  file.path(PROJ, ".."))

data_path <- function(file, required = TRUE) {
  hits <- file.path(DATA_DIRS, file)
  hit  <- hits[file.exists(hits)]
  if (!length(hit)) {
    if (!required) return(NA_character_)
    stop(sprintf("Could not find '%s'. Looked in:\n  %s",
                 file, paste(normalizePath(DATA_DIRS, mustWork = FALSE), collapse = "\n  ")),
         call. = FALSE)
  }
  normalizePath(hit[1])
}
OUT  <- file.path(PROJ, "output")
REF  <- file.path(PROJ, "ncdata_untrusted")   # downloaded Nextclade reference data (untrusted input)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

## Stata-style date strings ("21jul2021") parsed locale-independently.
parse_dmy <- function(x) {
  x <- tolower(trimws(as.character(x)))
  mon <- c(jan=1,feb=2,mar=3,apr=4,may=5,jun=6,jul=7,aug=8,sep=9,oct=10,nov=11,dec=12)
  d <- suppressWarnings(as.integer(str_sub(x, 1, 2)))
  m <- mon[str_sub(x, 3, 5)]
  y <- suppressWarnings(as.integer(str_sub(x, 6, 9)))
  out <- as.Date(NA)
  ok <- !is.na(d) & !is.na(m) & !is.na(y)
  out <- rep(as.Date(NA), length(x))
  out[ok] <- as.Date(sprintf("%04d-%02d-%02d", y[ok], m[ok], d[ok]))
  out
}

## RSV season label: Melghat RSV activity runs with the monsoon, so seasons are
## defined Jun-May and named by their starting year.
rsv_season <- function(date) {
  y <- as.integer(format(date, "%Y")); m <- as.integer(format(date, "%m"))
  ifelse(is.na(date), NA_character_, paste0(ifelse(m >= 6, y, y - 1), "/", substr(ifelse(m >= 6, y + 1, y), 3, 4)))
}

## Translate a nucleotide vector/string to amino acids using the standard
## genetic code. Avoids a dependency on seqinr for one function.
GENETIC_CODE <- local({
  b <- c("T", "C", "A", "G")
  codons <- paste0(rep(b, each = 16), rep(rep(b, each = 4), 4), rep(b, 16))
  aa <- strsplit(paste0("FFLLSSSSYY**CC*W", "LLLLPPPPHHQQRRRR",
                        "IIIMTTTTNNKKSSRR", "VVVVAAAADDEEGGGG"), "")[[1]]
  setNames(aa, codons)
})

translate_nt <- function(x) {
  v <- toupper(if (length(x) == 1 && nchar(x) > 1) strsplit(x, "")[[1]] else x)
  n <- length(v) %/% 3
  if (n == 0) return(character(0))
  vapply(seq_len(n), function(i) {
    c3 <- v[(3*i - 2):(3*i)]
    if (any(!c3 %in% c("A", "C", "G", "T"))) "X"
    else unname(GENETIC_CODE[paste0(c3, collapse = "")])
  }, character(1))
}

## Read an aligned FASTA into a character matrix (rows = sequences, cols = sites).
read_fasta_matrix <- function(path) {
  l <- readLines(path)
  h <- grep("^>", l)
  ids <- sub("^>\\s*", "", l[h])
  ends <- c(h[-1] - 1, length(l))
  sq <- mapply(function(a, b) paste(l[(a + 1):b], collapse = ""), h, ends)
  m <- do.call(rbind, strsplit(toupper(sq), ""))
  rownames(m) <- ids
  m
}

write_out <- function(x, name) {
  readr::write_csv(x, file.path(OUT, name))
  cat(sprintf("  -> wrote %s (%d rows)\n", name, nrow(x)))
  invisible(x)
}

save_plot <- function(p, name, w = 9, h = 6, dpi = 200) {
  ggplot2::ggsave(file.path(OUT, name), p, width = w, height = h, dpi = dpi)
  cat(sprintf("  -> wrote %s\n", name))
}

theme_melghat <- function() {
  theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(colour = "grey35", size = 9),
          strip.text = element_text(face = "bold"))
}

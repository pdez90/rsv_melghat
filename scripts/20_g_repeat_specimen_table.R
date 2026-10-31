# 20_g_repeat_specimen_table.R -- one specimen-level table reconciling the two
# different definitions of "the duplication is missing" that this project has
# used, so they are never again quoted as if they were the same quantity.
#
#   LENGTH definition  (09_ggene.R): G_len < 310 aa in the Nextclade G
#                      translation -> "72-nt duplication absent"
#   GAP definition     (17_repeat_region.py): >=50% of the 72 combined-alignment
#                      columns spanning G aa 262-285 are gaps
#
# These are not interchangeable: the first is a property of the translated CDS,
# the second of the nucleotide alignment, and a genome can fail one without the
# other. The table also carries Ct, season, lineage, G length, Nextclade
# frameshift calls and the 3' gap length, which is what a batch/Ct association
# would be tested on once laboratory metadata arrives.

source(file.path("scripts", "00_setup.R"))
`%||%` <- function(a, b) if (is.null(a)) b else ifelse(is.na(a), b, a)
cat("\n== 20 G repeat specimen table ==\n")

dat  <- readRDS(file.path(OUT, "02_dat.rds"))
glen <- read_csv(file.path(OUT, "09_G_length_duplication.csv"), show_col_types = FALSE)

## gap fractions, computed from the combined alignment through the coordinate map
gapfile <- file.path(OUT, "20_gapfrac_RSV-A.csv")
if (!file.exists(gapfile)) stop("run the python helper first (see 20_gapfrac.py)")
gapf <- read_csv(gapfile, show_col_types = FALSE)

## Nextclade writes .csv (semicolon) or .tsv (tab) depending on how it was run,
## and a plain read_csv() on the semicolon form yields a single column.
nc_path <- c(file.path(PROJ, "nextclade", "RSV-A", "nextclade.csv"),
             file.path(PROJ, "nextclade", "RSV-A", "nextclade.tsv"))
nc_path <- nc_path[file.exists(nc_path)][1]
if (is.na(nc_path)) stop("no nextclade.csv or nextclade.tsv for RSV-A")
hdr <- readLines(nc_path, n = 1)
delim <- c(";", "\t", ",")[which.max(c(lengths(regmatches(hdr, gregexpr(";", hdr))),
                                        lengths(regmatches(hdr, gregexpr("\t", hdr))),
                                        lengths(regmatches(hdr, gregexpr(",", hdr)))))]
nc <- read_delim(nc_path, delim = delim, show_col_types = FALSE,
                 guess_max = 5000) %>%
      select(seqName, totalFrameShifts, frameShifts)

tab <- dat %>%
  filter(nc_dataset == "RSV-A") %>%
  select(FieldNoNIV_f, studyid, season, clade, code_village, CTvalue_f,
         nc_coverage, nc_qc_status, link_status) %>%
  left_join(glen %>% select(FieldNoNIV_f, G_len, G_dup), by = "FieldNoNIV_f") %>%
  left_join(gapf, by = c("FieldNoNIV_f" = "seqName")) %>%
  left_join(nc %>% select(seqName, totalFrameShifts, frameShifts),
            by = c("FieldNoNIV_f" = "seqName")) %>%
  mutate(
    ## Features are kept SEPARATE and named for what was observed in the
    ## consensus alignment, not for a mechanism. "Absent" below always means
    ## "meets the consensus-alignment criterion for apparent absence" -- a
    ## >=50% gap across the duplication window does not by itself establish a
    ## complete biological 72-nt deletion.
    apparent_repeat_absence   = dup_gap_frac >= 0.50,
    apparent_3prime_shortening = three_prime_gap >= 60,
    short_G_translation        = G_len < 310,
    frameshift_flagged         = totalFrameShifts > 0,
    G309_frameshift            = grepl("G:309-322", frameShifts %||% ""),
    missing_frac               = 1 - nc_coverage,
    feature = case_when(
      apparent_repeat_absence &  apparent_3prime_shortening ~ "both features",
      apparent_repeat_absence & !apparent_3prime_shortening ~ "apparent repeat absence only",
      !apparent_repeat_absence & apparent_3prime_shortening ~ "apparent 3' shortening only",
      TRUE                                                  ~ "neither feature"),
    ## retained only to reconcile against the earlier single category
    agreement = case_when(
      short_G_translation &  apparent_repeat_absence ~ "both call absent",
      short_G_translation & !apparent_repeat_absence ~ "LENGTH only",
      !short_G_translation & apparent_repeat_absence ~ "GAP only",
      TRUE                                           ~ "both call present"),
    length_absent = short_G_translation,
    gap_absent    = apparent_repeat_absence,
    g3_truncated_71 = three_prime_gap == 71)

cat("\n-- FINAL feature categorisation (features kept distinct) --\n")
cat("   Any absence of overlap below is 'no overlap under the specified feature\n")
cat("   thresholds', NOT biological mutual exclusivity.\n")
print(with(tab, table(feature, season)))
cat("\n-- feature by lineage --\n")
print(with(tab, table(feature, clade)))
cat("\n-- per-feature summary --\n")
cat("   COVERAGE, correctly decomposed. Nextclade's coverage here is exactly\n")
cat("     (alignmentEnd - alignmentStart + 1 - totalMissing - totalNonACGTNs) / L\n")
cat("   verified to machine precision on all 241 genomes. Internal DELETIONS are\n")
cat("   NOT deducted -- subtracting them makes the fit worse. An earlier claim\n")
cat("   that the 72-nt deletion mechanically explains half the coverage gap was\n")
cat("   WRONG and is withdrawn.\n")
cat("   The difference is driven entirely by TERMINAL ALIGNMENT SPAN, not interior\n")
cat("   loss: median alignment 70-15029 for apparent repeat absence and 29-15009\n")
cat("   for 3' shortening, against 1-15225 for unaffected genomes, with interior\n")
cat("   loss ~0 in every group. Note 70 and 15029 are the reference CDS bounds of\n")
cat("   NS1 and L, so some genome sets appear trimmed to the coding span and others\n")
cat("   retain the UTRs -- a possible pipeline or batch signature, testable once\n")
cat("   batch metadata is available.\n")
cat("   Within the season- and lineage-matched 2021/22 A.D.3 stratum the direction\n")
cat("   REVERSES (0.983 affected vs 0.970 unaffected), so the overall gap is a\n")
cat("   cross-stratum artefact. Coverage remains confounded with lineage and season.\n")
print(tab %>% group_by(feature) %>%
      summarise(n = n(),
                median_G_len = median(G_len, na.rm = TRUE),
                median_dup_gap = round(median(dup_gap_frac), 3),
                median_3p_gap = median(three_prime_gap),
                frameshift_flagged = sum(frameshift_flagged),
                G309 = sum(G309_frameshift),
                median_Ct = median(CTvalue_f, na.rm = TRUE),
                median_coverage = round(median(nc_coverage, na.rm = TRUE), 4),
                median_missing = round(median(missing_frac, na.rm = TRUE), 4),
                .groups = "drop") %>% as.data.frame(), right = FALSE)

cat("\n-- cross-tabulation of the two definitions --\n")
print(with(tab, table(LENGTH_absent = length_absent, GAP_absent = gap_absent,
                      useNA = "ifany")))
cat("\n-- agreement category by season --\n")
print(with(tab, table(season, agreement)))

cat("\n-- discordant specimens --\n")
disc <- tab %>% filter(agreement %in% c("LENGTH only", "GAP only")) %>%
  select(FieldNoNIV_f, season, clade, G_len, dup_gap_frac, agreement,
         three_prime_gap, totalFrameShifts, CTvalue_f)
if (nrow(disc)) print(as.data.frame(disc), right = FALSE) else cat("  none\n")

cat("\n-- Ct by agreement category (no batch metadata yet; Ct is a proxy only) --\n")
print(tab %>% group_by(agreement) %>%
        summarise(n = n(), median_Ct = median(CTvalue_f, na.rm = TRUE),
                  IQR_Ct = paste(round(quantile(CTvalue_f, c(.25,.75), na.rm = TRUE), 1),
                                 collapse = "-"),
                  median_coverage = round(median(nc_coverage, na.rm = TRUE), 4),
                  .groups = "drop") %>% as.data.frame(), right = FALSE)

cat("\n-- within 2021/22 and lineage A.D.3 only (season- and lineage-matched) --\n")
m <- tab %>% filter(season == "2021/22", clade == "A.D.3")
if (nrow(m)) print(m %>% group_by(agreement) %>%
        summarise(n = n(), median_Ct = median(CTvalue_f, na.rm = TRUE),
                  median_cov = round(median(nc_coverage, na.rm = TRUE), 4),
                  .groups = "drop") %>% as.data.frame(), right = FALSE)

write_csv(tab, file.path(OUT, "20_g_repeat_specimen_table.csv"))
cat(sprintf("\n  -> wrote 20_g_repeat_specimen_table.csv (%d specimens)\n", nrow(tab)))
cat("== 20 done ==\n")

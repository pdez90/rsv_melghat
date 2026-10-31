# 02_nextclade_results.R -- ingest Nextclade Web output (RSV-A and RSV-B runs).
#
# Nextclade v3.28.0 was run in the browser against the Nextstrain RSV-A
# (EPI_ISL_412866) and RSV-B (EPI_ISL_1653999) reference datasets. These calls
# are authoritative; 02_clade_assign.R is an independent R re-implementation
# kept only as a cross-check.

source(file.path("scripts", "00_setup.R"))
cat("\n== 02 Nextclade results ==\n")

NCA <- file.path(PROJ, "nextclade", "RSV-A")
NCB <- file.path(PROJ, "nextclade", "RSV-B")
linked <- readRDS(file.path(OUT, "01_linked.rds"))

read_nc <- function(dir, type) {
  read_tsv(file.path(dir, "nextclade.tsv"), show_col_types = FALSE,
           guess_max = 500, progress = FALSE) %>%
    transmute(
      FieldNoNIV_f        = seqName,
      nc_clade            = clade,
      nc_escape           = suppressWarnings(as.numeric(escape)),
      nc_qc_status        = qc.overallStatus,
      nc_qc_score         = qc.overallScore,
      nc_coverage         = coverage,
      nc_subs             = totalSubstitutions,
      nc_missing          = totalMissing,
      nc_aa_subs_n        = totalAminoacidSubstitutions,
      nc_aa_subs          = aaSubstitutions,
      nc_private_subs     = `privateNucMutations.totalPrivateSubstitutions`,
      nc_clade_founder    = `founderMuts['clade'].nodeName`,
      nc_dataset          = type)
}

nc <- bind_rows(read_nc(NCA, "RSV-A"), read_nc(NCB, "RSV-B"))
cat(sprintf("Nextclade rows: %d (A %d / B %d)\n", nrow(nc),
            sum(nc$nc_dataset == "RSV-A"), sum(nc$nc_dataset == "RSV-B")))
stopifnot(nrow(nc) == 241, !any(duplicated(nc$FieldNoNIV_f)))

## ---- QC ---------------------------------------------------------------------
cat("\n-- Nextclade QC --\n"); print(table(nc$nc_qc_status, nc$nc_dataset))
cat(sprintf("coverage: median %.3f, min %.3f; %d genomes <90%% covered\n",
            median(nc$nc_coverage), min(nc$nc_coverage), sum(nc$nc_coverage < 0.90)))

## ---- merge onto the cohort-linked dataset -----------------------------------
dat <- linked %>%
  left_join(nc, by = "FieldNoNIV_f") %>%
  mutate(supplied_clade = str_trim(CLADE_f),
         clade = nc_clade,
         clade_agrees = case_when(
           is.na(clade) ~ NA,
           supplied_clade == clade ~ TRUE,
           supplied_clade %in% c("A.D.", "A.D") & str_starts(clade, "A.D") ~ TRUE,
           TRUE ~ FALSE),
         newly_resolved = supplied_clade %in% c("A.D.", "A.D") & !clade %in% c("A.D.", "A.D"))

cat("\n-- supplied CLADE_f vs Nextclade clade --\n")
print(table(supplied = dat$supplied_clade, nextclade = dat$clade, useNA = "ifany"))
cat(sprintf("\nexact agreement: %d/%d (%.1f%%)\n",
            sum(dat$supplied_clade == dat$clade, na.rm = TRUE), nrow(dat),
            100 * mean(dat$supplied_clade == dat$clade, na.rm = TRUE)))
cat(sprintf("compatible (treating bare 'A.D.' as unresolved): %.1f%%\n",
            100 * mean(dat$clade_agrees, na.rm = TRUE)))
cat(sprintf("genomes whose bare 'A.D.' is now resolved to a finer lineage: %d\n",
            sum(dat$newly_resolved, na.rm = TRUE)))
if (any(dat$newly_resolved, na.rm = TRUE))
  print(dat %>% filter(newly_resolved) %>% count(season, clade))

## ---- subgroup check against the cohort's own A/B typing ---------------------
cat("\n-- subgroup vs ari.csv typing --\n")
print(with(dat %>% filter(linked), table(nextclade = nc_dataset,
                                         ari = ifelse(rsva_episode == 1, "A", "B"), useNA = "ifany")))

## ---- immune escape score (Nextclade RSV F escape model) ---------------------
if (any(!is.na(dat$nc_escape))) {
  cat("\n-- Nextclade escape score --\n")
  print(dat %>% group_by(nc_dataset) %>%
          summarise(n = sum(!is.na(nc_escape)), median = median(nc_escape, na.rm = TRUE),
                    min = min(nc_escape, na.rm = TRUE), max = max(nc_escape, na.rm = TRUE),
                    .groups = "drop"))
}

write_out(dat %>% select(-FASTA_f), "02_nextclade_linked.csv")
saveRDS(dat, file.path(OUT, "02_dat.rds"))
cat("== 02 done ==\n")

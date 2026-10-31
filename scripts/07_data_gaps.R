# 07_data_gaps.R -- a complete, auditable account of what was lost from the
# genotypic analysis and why. Produces the tables behind the data-gaps document.

source(file.path("scripts", "00_setup.R"))
cat("\n== 07 data gaps ==\n")

dat   <- readRDS(file.path(OUT, "02_dat.rds"))
ari   <- read_csv(data_path("ari.csv"), show_col_types = FALSE, guess_max = 20000)
demog <- read_csv(data_path("demographics.csv"), show_col_types = FALSE,
                  guess_max = 20000) %>% select(studyid, birth_date_new)

## ---- per-sample disposition --------------------------------------------------
disp <- dat %>%
  mutate(
    in_demographics = studyid %in% demog$studyid,
    episodes_any    = map_int(studyid, ~ sum(ari$studyid == .x)),
    reason = case_when(
      link_status == "exact" ~ "Linked - exact episode-date match",
      link_status == "nearest_prior_episode_same_RSV_type" ~
        "Linked - nearest prior episode, same subgroup (flagged)",
      link_status == "unlinked_age_ge_2y_outside_LRTI_cohort" ~
        "EXCLUDED - child aged 2y or over at swab; outside the under-2 cohort",
      link_status == "unlinked_not_in_demographics" ~
        "EXCLUDED - studyid absent from demographics.csv",
      link_status == "unlinked_no_matching_LRTI_episode" ~
        "EXCLUDED - child enrolled but no LRTI episode matches the swab date",
      TRUE ~ link_status)) %>%
  select(FieldNoNIV_f, studyid, swab_date, season, nc_dataset, clade, CTvalue_f,
         age_days_swab, in_demographics, episodes_any, nc_qc_status, nc_coverage,
         link_status, reason, linked)

cat("\n-- disposition of all 241 sequenced samples --\n")
summ <- disp %>% count(reason, name = "samples") %>% arrange(desc(samples))
print(as.data.frame(summ))

excluded <- disp %>% filter(!linked) %>% arrange(reason, swab_date)
cat(sprintf("\n-- the %d excluded samples, in full --\n", nrow(excluded)))
print(as.data.frame(excluded %>% select(FieldNoNIV_f, studyid, swab_date, season,
                                        nc_dataset, clade, age_days_swab,
                                        in_demographics, episodes_any, reason)))

## do the 6 "enrolled but no matching episode" children have ANY episodes?
cat("\n-- children enrolled but with no episode matching the swab --\n")
noep <- disp %>% filter(link_status == "unlinked_no_matching_LRTI_episode")
for (i in seq_len(nrow(noep))) {
  eps <- ari %>% filter(studyid == noep$studyid[i]) %>%
    mutate(d = parse_dmy(ai_date_new)) %>% arrange(d)
  cat(sprintf("%s (%s, swab %s): %d episode(s) on record%s\n",
              noep$FieldNoNIV_f[i], noep$studyid[i], noep$swab_date[i], nrow(eps),
              if (nrow(eps)) paste0(" - ", paste(format(eps$d, "%d %b %Y"), collapse = ", ")) else ""))
}

## ---- the 28 over-age samples -------------------------------------------------
over <- disp %>% filter(link_status == "unlinked_age_ge_2y_outside_LRTI_cohort")
cat(sprintf("\n-- over-age samples: n=%d, all season %s --\n", nrow(over),
            paste(unique(over$season), collapse = ", ")))
cat(sprintf("age at swab: median %.1f y (range %.1f-%.1f y)\n",
            median(over$age_days_swab)/365.25, min(over$age_days_swab)/365.25,
            max(over$age_days_swab)/365.25))
cat(sprintf("these are %d of the %d genomes sequenced from 2021/22 (%.0f%%)\n",
            nrow(over), sum(dat$season == "2021/22", na.rm = TRUE),
            100*nrow(over)/sum(dat$season == "2021/22", na.rm = TRUE)))

## ---- sequence-quality losses --------------------------------------------------
cat("\n-- sequence quality --\n")
print(dat %>% count(nc_dataset, nc_qc_status) %>%
        pivot_wider(names_from = nc_qc_status, values_from = n, values_fill = 0) %>% as.data.frame())
lowcov <- dat %>% filter(nc_coverage < 0.90) %>%
  select(FieldNoNIV_f, nc_dataset, clade, nc_coverage, nc_missing, nc_qc_status)
cat("\ngenomes under 90% coverage:\n"); print(as.data.frame(lowcov))

## frameshifts (affect protein-level analysis only)
fs <- bind_rows(
  read_tsv(file.path(PROJ, "nextclade", "RSV-A", "nextclade.tsv"),
           show_col_types = FALSE, guess_max = 500, progress = FALSE) %>% mutate(type = "RSV-A"),
  read_tsv(file.path(PROJ, "nextclade", "RSV-B", "nextclade.tsv"),
           show_col_types = FALSE, guess_max = 500, progress = FALSE) %>% mutate(type = "RSV-B")) %>%
  filter(totalFrameShifts > 0) %>%
  transmute(FieldNoNIV_f = seqName, type, clade, totalFrameShifts, frameShifts, coverage)
cat(sprintf("\ngenomes with frameshifts: %d (all %s)\n", nrow(fs),
            paste(unique(fs$type), collapse = "/")))
cat("affected CDS:\n")
print(table(unlist(str_extract_all(fs$frameShifts, "^[^:]+|(?<=,)[^:]+"))))
f_fs <- fs %>% filter(str_detect(frameShifts, "(^|,)F:"))
cat(sprintf("of these, %d has a frameshift in F and is unusable for the F-protein analysis: %s\n",
            nrow(f_fs), paste(f_fs$FieldNoNIV_f, collapse = ", ")))

## ---- clade disagreements -------------------------------------------------------
disagree <- dat %>%
  filter(str_remove(supplied_clade, "\\.$") != str_remove(clade, "\\.$")) %>%
  select(FieldNoNIV_f, studyid, season, supplied_clade, clade, nc_qc_status, nc_coverage)
cat(sprintf("\n-- clade disagreements, supplied vs Nextclade: %d --\n", nrow(disagree)))
print(as.data.frame(disagree))

## ---- coverage of the cohort ----------------------------------------------------
rsv_eps <- ari %>% filter(rsva_episode == 1 | rsvb_episode == 1)
cat(sprintf("\n-- denominator --\ncohort RSV+ episodes: %d; sequenced and linked: %d (%.1f%%)\n",
            nrow(rsv_eps), sum(dat$linked), 100*sum(dat$linked)/nrow(rsv_eps)))
cat(sprintf("Ct of sequenced samples: median %.0f, range %.0f-%.0f (selection on low Ct)\n",
            median(dat$CTvalue_f, na.rm = TRUE), min(dat$CTvalue_f, na.rm = TRUE),
            max(dat$CTvalue_f, na.rm = TRUE)))
cat(sprintf("LRTI deaths among sequenced episodes: %d\n",
            sum(dat$lrti_death == 1, na.rm = TRUE)))

write_out(disp,     "07_sample_disposition.csv")
write_out(excluded, "07_excluded_samples.csv")
write_out(summ,     "07_disposition_summary.csv")
write_out(fs,       "07_frameshift_genomes.csv")
write_out(disagree, "07_clade_disagreements.csv")
cat("== 07 done ==\n")

# 27_cohort_status_audit.R -- reconcile two different RSV-positive definitions
# before any coverage or ascertainment number is quoted.
#
# THE DISCREPANCY. A sequenced specimen carries a laboratory result
# (RSVResult_f, which agrees with the Nextclade subgroup for all 241 genomes).
# A cohort episode carries its own flags (rsva_episode, rsvb_episode). For three
# linked episodes these disagree. Linkage success (203 of 241) is unaffected;
# COVERAGE OF COHORT-POSITIVE EPISODES is, because the denominator and the
# numerator must use the same definition.
#
# DECLARED DEFINITION, used for all coverage and ascertainment reporting:
#   an episode is RSV-positive if the COHORT flags say so.
# This keeps numerator and denominator on one definition and leaves the cohort
# data as recorded. Nothing is silently recoded. The discordant genomes remain
# in the genomic analyses, where the laboratory result is the relevant one, and
# are flagged there.

source(file.path("scripts", "00_setup.R"))
cat("\n== 27 cohort RSV-status audit ==\n")

dat <- readRDS(file.path(OUT, "02_dat.rds"))
ari <- read_csv(data_path("ari.csv"), show_col_types = FALSE, guess_max = 20000) %>%
  mutate(ep_date = parse_dmy(ai_date_new),
         season = rsv_season(ep_date),
         swabbed = !is.na(nswab_episode) & nswab_episode == 1,
         cohort_a = !is.na(rsva_episode) & rsva_episode == 1,
         cohort_b = !is.na(rsvb_episode) & rsvb_episode == 1,
         cohort_pos = cohort_a | cohort_b)

L <- dat %>% filter(link_status %in% c("exact", "nearest_prior_episode_same_RSV_type")) %>%
  left_join(ari %>% select(studyid, ep_date, cohort_a, cohort_b, cohort_pos, swabbed),
            by = c("studyid", "ep_date")) %>%
  mutate(cohort_call = case_when(cohort_a & cohort_b ~ "both flags set",
                                 cohort_a ~ "RSV-A", cohort_b ~ "RSV-B",
                                 TRUE ~ "neither flag set"),
         concordant = cohort_call == nc_dataset)

cat(sprintf("\n  linked genomes                       %3d\n", nrow(L)))
cat(sprintf("  laboratory subgroup == cohort flag   %3d\n", sum(L$concordant)))
cat(sprintf("  DISCORDANT                           %3d\n", sum(!L$concordant)))

disc <- L %>% filter(!concordant) %>%
  select(FieldNoNIV_f, studyid, ep_date, season, lab_subgroup = nc_dataset,
         cohort_call, link_status, swabbed, lri_cat)
cat("\n-- discordant records (reported, not recoded) --\n")
print(as.data.frame(disc), right = FALSE)
write_csv(disc, file.path(OUT, "27_discordant_cohort_status.csv"))

## coverage under the DECLARED (cohort-flag) definition
cov <- ari %>% filter(swabbed, cohort_pos) %>% count(season, name = "cohort_positive") %>%
  left_join(L %>% filter(cohort_pos) %>% count(season, name = "sequenced_cohort_positive"),
            by = "season") %>%
  left_join(L %>% count(season, name = "sequenced_linked_any"), by = "season") %>%
  mutate(across(where(is.integer), ~coalesce(.x, 0L)),
         pct_declared = round(100 * sequenced_cohort_positive / cohort_positive, 1))
tot <- cov %>% summarise(season = "TOTAL", across(where(is.numeric), sum)) %>%
  mutate(pct_declared = round(100 * sequenced_cohort_positive / cohort_positive, 1))
cov2 <- bind_rows(cov, tot)
cat("\n-- coverage under the declared cohort-flag definition --\n")
print(as.data.frame(cov2), right = FALSE)
write_csv(cov2, file.path(OUT, "27_coverage_declared_definition.csv"))

cat(sprintf("\n  Linkage success stays %d of %d sequenced specimens (%.1f%%).\n",
            nrow(L), nrow(dat), 100 * nrow(L) / nrow(dat)))
cat(sprintf("  Coverage of cohort-positive episodes is %d/%d = %.1f%%,\n",
            sum(L$cohort_pos), sum(ari$swabbed & ari$cohort_pos),
            100 * sum(L$cohort_pos) / sum(ari$swabbed & ari$cohort_pos)))
cat("  NOT 203/644: those two counts use different definitions of RSV-positive.\n")

## ascertainment populations must use the same definition
cat("\n-- effect on the ascertainment populations --\n")
cat(sprintf("  RSV-positive tested episodes (denominator) %4d\n", sum(ari$swabbed & ari$cohort_pos)))
cat(sprintf("  of which sequenced and linked              %4d\n", sum(L$cohort_pos)))
cat(sprintf("  linked but NOT cohort-positive             %4d  (kept in genomic analyses only)\n",
            sum(!L$cohort_pos)))
cat("== 27 done ==\n")

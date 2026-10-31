# 01_link_sequencing.R -- link the 241 sequenced RSV swabs to MAHAN cohort LRTI episodes.
#
# Keys (per Rowena Crow, 3 Nov 2025): link on *studyid*; the migration uid columns
# (uid_original / uid_current / uid_after*move) were never updated correctly and
# are not usable. Episode identity = studyid + ai_date_new.

source(file.path("scripts", "00_setup.R"))
cat("\n== 01 link sequencing to cohort ==\n")

seqdat <- read_csv(data_path("NP_results_with_clade_sequencing.csv"),
                   show_col_types = FALSE, guess_max = 1000)
ari    <- read_csv(data_path("ari.csv"), show_col_types = FALSE, guess_max = 20000)
deaths <- read_csv(data_path("lrtideaths.csv"), show_col_types = FALSE, guess_max = 5000)
demog  <- read_csv(data_path("demographics.csv"), show_col_types = FALSE, guess_max = 20000) %>%
  select(studyid, birth_date_new, code_village_93, tahsil, sexr_dem = sexr,
         wealthscore_rural, quintile_rural, mother_edu_categorized, malnutritionr)

stopifnot(nrow(seqdat) == 241, nrow(ari) == 9236, nrow(deaths) == 205)
cat(sprintf("inputs: %d sequenced swabs | %d AI episodes | %d LRTI deaths | %d enrolled children\n",
            nrow(seqdat), nrow(ari), nrow(deaths), nrow(demog)))

seqdat <- seqdat %>%
  mutate(swab_date = parse_dmy(swab_date_new),
         ai_date   = parse_dmy(ai_date_new),
         rsv_type  = case_when(rsva == 1 ~ "RSV-A", rsvb == 1 ~ "RSV-B", TRUE ~ NA_character_))
ari <- ari %>% mutate(ep_date = parse_dmy(ai_date_new), ep_row = row_number())

# ---- exact link: same child, same episode date -------------------------------
exact <- seqdat %>%
  select(FieldNoNIV_f, studyid, swab_date) %>%
  # many-to-many is expected here: a child can have several episodes
  inner_join(ari %>% select(studyid, ep_date, ep_row), by = "studyid",
             relationship = "many-to-many") %>%
  filter(ep_date == swab_date) %>%
  mutate(link = "exact", days_episode_to_swab = 0L)

# ---- near link: swab taken during an episode that began up to 45 d earlier ----
# Restricted to episodes whose own RSV typing matches the sequenced type, so we
# never attach a swab to an unrelated illness episode.
near <- seqdat %>%
  select(FieldNoNIV_f, studyid, swab_date, rsva, rsvb) %>%
  inner_join(ari %>% select(studyid, ep_date, ep_row, rsva_episode, rsvb_episode),
             by = "studyid", relationship = "many-to-many") %>%
  mutate(gap = as.integer(swab_date - ep_date),
         type_ok = (rsva == 1 & rsva_episode == 1) | (rsvb == 1 & rsvb_episode == 1)) %>%
  filter(gap > 0, gap <= 45, type_ok, !FieldNoNIV_f %in% exact$FieldNoNIV_f) %>%
  group_by(FieldNoNIV_f) %>% slice_min(gap, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(FieldNoNIV_f, studyid, swab_date, ep_row,
            link = "nearest_prior_episode_same_RSV_type", days_episode_to_swab = gap)

links <- bind_rows(exact %>% select(FieldNoNIV_f, ep_row, link, days_episode_to_swab),
                   near  %>% select(FieldNoNIV_f, ep_row, link, days_episode_to_swab))

ep_keep <- ari %>%
  select(ep_row, episode_start_date = ai_date_new, ep_date, lri_cat, community, hosp_lri,
         referred_hosp_episode, age_days_episode, sexr, ai_year, ai_month, code_village,
         subset_33_villages, timeline, lrti_death, died_lt_1d, recover_fu1_episode,
         any_pos_viral_resp, rsva_episode, rsvb_episode, rhino_episode, hmpv_episode,
         piv_episode, adv_episode, flua_episode, flub_episode, sars_cov2_episode,
         weight_kg, muac, nswab_episode)

linked <- seqdat %>%
  left_join(links, by = "FieldNoNIV_f") %>%
  left_join(ep_keep, by = "ep_row") %>%
  left_join(demog, by = "studyid") %>%
  mutate(
    birth_date    = parse_dmy(birth_date_new),
    age_days_swab = as.integer(swab_date - birth_date),
    season        = rsv_season(swab_date),
    link_status = case_when(
      !is.na(link)                 ~ link,
      is.na(birth_date)            ~ "unlinked_not_in_demographics",
      age_days_swab >= 730         ~ "unlinked_age_ge_2y_outside_LRTI_cohort",
      TRUE                         ~ "unlinked_no_matching_LRTI_episode"),
    linked = link_status %in% c("exact", "nearest_prior_episode_same_RSV_type"),
    severity = factor(lri_cat, levels = c(1, 2, 3),
                      labels = c("very severe", "severe", "non-severe")),
    reached_facility = if_else(community == 0, 1L, 0L))

# ---- death linkage -----------------------------------------------------------
linked <- linked %>%
  left_join(deaths %>% select(studyid, death_location, death_date_new, cod,
                              days_to_death_episode, death_age_days), by = "studyid")

cat("\nlink status:\n"); print(table(linked$link_status))
cat(sprintf("\nlinked to an episode: %d / %d (%.1f%%)\n",
            sum(linked$linked), nrow(linked), 100 * mean(linked$linked)))

# A/B typing concordance between the sequencing file and the cohort episode record
conc <- linked %>% filter(linked) %>%
  summarise(a_agree = mean(rsva == rsva_episode, na.rm = TRUE),
            b_agree = mean(rsvb == rsvb_episode, na.rm = TRUE))
cat(sprintf("A/B concordance with ari.csv: RSV-A %.1f%%  RSV-B %.1f%%\n",
            100 * conc$a_agree, 100 * conc$b_agree))

# sequencing coverage of all RSV-positive episodes in the cohort
rsv_eps <- ari %>% filter(rsva_episode == 1 | rsvb_episode == 1)
cat(sprintf("cohort RSV+ episodes: %d (A %d / B %d); sequenced & linked: %d (%.0f%%)\n",
            nrow(rsv_eps), sum(rsv_eps$rsva_episode == 1, na.rm = TRUE),
            sum(rsv_eps$rsvb_episode == 1, na.rm = TRUE),
            sum(linked$linked), 100 * sum(linked$linked) / nrow(rsv_eps)))

# reinfections: same child sequenced more than once
reinf <- linked %>% group_by(studyid) %>% filter(n() > 1) %>%
  arrange(studyid, swab_date) %>%
  select(studyid, FieldNoNIV_f, swab_date, season, rsv_type, CLADE_f,
         severity, age_days_swab, link_status) %>% ungroup()
cat(sprintf("\nchildren sequenced more than once: %d (%d swabs)\n",
            n_distinct(reinf$studyid), nrow(reinf)))
print(as.data.frame(reinf))

write_out(linked, "01_sequencing_linked.csv")
write_out(reinf,  "01_repeat_sequenced_children.csv")
saveRDS(linked, file.path(OUT, "01_linked.rds"))
cat("== 01 done ==\n")

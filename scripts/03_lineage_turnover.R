# 03_lineage_turnover.R -- RSV subgroup and lineage turnover across seasons,
# the timing of RSV activity, and how representative the sequenced subset is.

source(file.path("scripts", "00_setup.R"))
cat("\n== 03 lineage turnover and seasonality ==\n")

dat <- readRDS(file.path(OUT, "02_dat.rds"))
ari <- read_csv(data_path("ari.csv"), show_col_types = FALSE, guess_max = 20000) %>%
  mutate(ep_date = parse_dmy(ai_date_new), season = rsv_season(ep_date),
         # NB rsva_episode/rsvb_episode are coded 0 for episodes that were never
         # swabbed, so positivity must be computed over SWABBED episodes only.
         swabbed = !is.na(nswab_episode) & nswab_episode == 1,
         rsv_pos = rsva_episode == 1 | rsvb_episode == 1,
         subgroup = case_when(rsva_episode == 1 ~ "RSV-A", rsvb_episode == 1 ~ "RSV-B"),
         # row number of the UNFILTERED read, matching 01_link_sequencing.R,
         # so dat$ep_row refers to the same episodes
         ep_row = row_number())

# Fixed categorical order (validated palette, adjacent pairlist)
LIN <- c("A.D", "A.D.2.1", "A.D.2.2", "A.D.3", "A.D.5", "B.D.4.1", "B.D.4.1.1")
PAL <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"), LIN)

d <- dat %>% filter(!is.na(clade)) %>%
  mutate(lineage = factor(clade, levels = LIN),
         season  = factor(season, levels = sort(unique(season))))

## ---- 1. lineage composition by season ---------------------------------------
by_season <- d %>% count(season, nc_dataset, lineage, name = "n") %>%
  group_by(season) %>% mutate(pct_of_season = round(100 * n / sum(n), 1)) %>% ungroup()
cat("\n-- sequenced genomes by season and lineage --\n")
print(as.data.frame(by_season))

# Coverage of the clinical cohort must count LINKED EPISODES, not all genomes:
# 38 sequenced specimens have no cohort episode, and they are not spread evenly
# across seasons (28 of them fall in 2021/22 alone).
season_tot <- dat %>% filter(linked) %>%
  group_by(season) %>% summarise(linked_episodes = n_distinct(ep_row), .groups = "drop")
# DECLARED DEFINITION (see 27_cohort_status_audit.R): coverage is reported on the
# cohort's own rule -- an episode is RSV-positive if the cohort flags say so --
# so the numerator and the denominator below share one definition. Two linked
# episodes are cohort-negative against a positive laboratory result; they are
# excluded from the coverage numerator and retained in the genomic analyses,
# where the laboratory result is the relevant one. linked_episodes is kept
# alongside so the two numerators can be read against each other.
cohort_pos_rows <- ari %>% filter(rsv_pos) %>% pull(ep_row)
season_declared <- dat %>% filter(linked, ep_row %in% cohort_pos_rows) %>%
  group_by(season) %>%
  summarise(linked_cohort_positive = n_distinct(ep_row), .groups = "drop")
all_specimens <- d %>% count(season, name = "all_sequenced")
cat("\n-- subgroup dominance by season --\n")
dom <- d %>% count(season, nc_dataset) %>% group_by(season) %>%
  mutate(pct = round(100 * n / sum(n))) %>% ungroup()
print(as.data.frame(dom))

## ---- 2. all cohort RSV+ episodes by season (the full denominator) -----------
cohort_season <- ari %>% filter(rsv_pos) %>% count(season, subgroup, name = "episodes") %>%
  pivot_wider(names_from = subgroup, values_from = episodes, values_fill = 0) %>%
  mutate(total = `RSV-A` + `RSV-B`) %>%
  left_join(season_tot, by = "season") %>%
  left_join(season_declared, by = "season") %>%
  left_join(all_specimens, by = "season") %>%
  mutate(linked_episodes = replace_na(linked_episodes, 0L),
         linked_cohort_positive = replace_na(linked_cohort_positive, 0L),
         all_sequenced   = replace_na(all_sequenced, 0L),
         pct_cohort_covered = round(100 * linked_cohort_positive / total, 1))
cat("\n-- cohort RSV+ episodes vs LINKED sequenced episodes, by season --\n")
print(as.data.frame(cohort_season))
cat("\npct_cohort_covered = linked AND cohort-positive episodes / cohort RSV+ episodes,\n")
cat("the declared definition. 'linked_episodes' is the wider all-linked numerator and\n")
cat("differs by the two cohort-negative episodes. 'all_sequenced'\n")
cat("counts every specimen from that season including those with no cohort episode;\n")
cat("dividing it by cohort episodes overstates coverage and is NOT done here.\n")
cat("\nNOTE: the 2021/22 season covers the 33-village subset only; earlier seasons cover 93 villages.\n")

## ---- 3. timing within the year ----------------------------------------------
monthly <- ari %>% filter(rsv_pos, !is.na(ep_date)) %>%
  mutate(month = as.integer(format(ep_date, "%m"))) %>%
  count(month, subgroup, name = "episodes")
cat("\n-- RSV+ episodes by calendar month (all cohort) --\n")
print(monthly %>% pivot_wider(names_from = subgroup, values_from = episodes, values_fill = 0) %>%
        mutate(total = `RSV-A` + `RSV-B`) %>% as.data.frame())
monsoon <- ari %>% filter(!is.na(ep_date)) %>%
  mutate(month = as.integer(format(ep_date, "%m")), mons = month %in% 6:9) %>%
  group_by(mons) %>%
  summarise(episodes = n(), swabbed = sum(swabbed), rsv = sum(rsv_pos, na.rm = TRUE),
            pct_swabbed = round(100 * sum(swabbed) / n(), 1),
            pct_rsv_of_swabbed = round(100 * sum(rsv_pos, na.rm = TRUE) / sum(swabbed), 1),
            .groups = "drop")
cat("\n-- RSV positivity, monsoon (Jun-Sep) vs rest (denominator = SWABBED) --\n")
print(as.data.frame(monsoon))
cat("Swabbing is not uniform: 52.3% of monsoon episodes vs 39.9% otherwise, and is\n")
cat("strongly selected on severity, so positivity is among tested episodes only.\n")

## ---- 4. is the sequenced subset representative? -----------------------------
rsv_eps <- ari %>% filter(rsv_pos) %>%
  mutate(key = paste(studyid, ai_date_new),
         sequenced = key %in% paste(dat$studyid[dat$linked], dat$episode_start_date[dat$linked]))
cat("\n-- sequenced vs unsequenced RSV+ episodes --\n")
rep_tab <- rsv_eps %>% group_by(sequenced) %>%
  summarise(n = n(),
            pct_very_severe = round(100 * mean(lri_cat == 1, na.rm = TRUE), 1),
            pct_severe      = round(100 * mean(lri_cat == 2, na.rm = TRUE), 1),
            pct_reached_facility = round(100 * mean(community == 0, na.rm = TRUE), 1),
            median_age_days = median(age_days_episode, na.rm = TRUE),
            deaths = sum(lrti_death == 1, na.rm = TRUE), .groups = "drop")
print(as.data.frame(rep_tab))
sev_test <- suppressWarnings(chisq.test(table(rsv_eps$sequenced, rsv_eps$lri_cat)))
cat(sprintf("severity distribution, sequenced vs not: chi-sq p = %.3f\n", sev_test$p.value))

## ---- 5. figure: lineage composition by season -------------------------------
p <- ggplot(by_season, aes(season, n, fill = lineage)) +
  geom_col(width = 0.72, colour = "#fcfcfb", linewidth = 0.6) +
  geom_text(aes(label = ifelse(n >= 3, n, "")), position = position_stack(vjust = 0.5),
            size = 3.1, colour = "white", fontface = "bold") +
  scale_fill_manual(values = PAL, drop = FALSE, name = "Lineage") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(title = "RSV lineage turnover across seasons, Melghat 2016-2022",
       subtitle = paste("241 sequenced genomes. Each season is dominated by a single subgroup;",
                        "\n2021/22 covers the 33-village subset only."),
       x = "RSV season (Jun-May)", y = "Sequenced genomes") +
  theme_melghat() + theme(legend.position = "right")
save_plot(p, "03_lineage_by_season.png", w = 9, h = 5.5)

write_out(by_season,     "03_lineage_by_season.csv")
write_out(cohort_season, "03_cohort_rsv_by_season.csv")
write_out(monthly,       "03_rsv_by_month.csv")
write_out(rep_tab,       "03_sequencing_representativeness.csv")
cat("== 03 done ==\n")

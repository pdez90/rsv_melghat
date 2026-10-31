## 00_verify_data.R -- Melghat: confirm we are on the correct / latest data
## Checks file identity, row counts, and every headline number in the two manuscripts.
source(file.path("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(tools)})

say <- function(...) cat(sprintf(...), "\n")
chk <- function(label, got, want) {
  ok <- isTRUE(all.equal(got, want))
  say("%-56s %-18s %s", label, paste(got, collapse = "/"),
      if (ok) "OK" else paste0("MISMATCH (expected ", paste(want, collapse = "/"), ")"))
  invisible(ok)
}

say("=============== FILE IDENTITY ===============")
WANT <- c("ari.csv", "lrtideaths.csv", "demographics.csv",
          "NP_results_with_clade_sequencing.csv")
paths <- vapply(WANT, data_path, character(1))
id <- data.frame(file = WANT, bytes = file.size(paths),
                 md5 = substr(md5sum(paths), 1, 12),
                 folder = basename(dirname(paths)), row.names = NULL)
print(id)

## The zip Rowena sent on 27 Oct 2025 shipped a mislabelled deaths file: a byte-identical
## copy of the AI-episodes file. She re-sent the correct 41K export on 3 Nov 2025.
## That corrected export has a known checksum, so confirm against it directly.
MD5_CORRECTED_DEATHS <- "8bb5f66f4bda4f75a3583dae45e7511a"
MD5_AI_EPISODES      <- "041dc79e82079851d27e08a65a9fb2e3"
dmd5 <- unname(md5sum(data_path("lrtideaths.csv")))
chk("deaths file == Rowena's 3 Nov 2025 corrected export", dmd5, MD5_CORRECTED_DEATHS)
if (identical(dmd5, MD5_AI_EPISODES))
  say("  !! This is the mislabelled copy of the AI-episodes file. Use Final/Data/lrtideaths.csv.")

## The stale copy still sitting in Final/LRTI data files/ should not be used.
stale <- data_path("LRTI deaths ALL YEARS.csv", required = FALSE)
if (!is.na(stale) && identical(unname(md5sum(stale)), MD5_AI_EPISODES))
  say("NOTE: %s is the mislabelled duplicate (= AI episodes file). Do not use it.", stale)

ari <- read_csv(data_path("ari.csv"), show_col_types = FALSE, progress = FALSE)
dth <- read_csv(data_path("lrtideaths.csv"), show_col_types = FALSE, progress = FALSE)
dem <- read_csv(data_path("demographics.csv"), show_col_types = FALSE, progress = FALSE)
seq <- read_csv(data_path("NP_results_with_clade_sequencing.csv"), show_col_types = FALSE,
                progress = FALSE)

say("")
say("=============== COUNTS vs MANUSCRIPTS ===============")
chk("AI episodes (Rowena email: 9,236)",            nrow(ari), 9236L)
chk("Enrolled children (LRI_Mortality abstract)",   nrow(dem), 14649L)
chk("LRTI deaths (Rowena email: 205)",              nrow(dth), 205L)
chk("Children with >=1 LRTI episode",               n_distinct(ari$studyid), 5750L)

sev <- table(factor(ari$lri_cat, 1:3))
chk("Episodes v.severe/severe/non-severe (Table 1)",
    as.integer(sev[c(1, 2, 3)]), c(1874L, 2621L, 4741L))

## Child-level severity: ari.csv carries TWO severity definitions that do not agree.
##   (a) lri_cat            1=very severe, 2=severe, 3=non-severe
##   (b) the WHO flags      v_severe_lri / severe_lri_who / non_severe_lri_who
## They disagree on 276 episodes (lri_cat==1 but v_severe_lri==0) and 377 episodes
## (lri_cat==2 but severe_lri_who==0). The published child-level Table 1 numbers
## (1664 / 2189 / 3274) are only reproducible by MIXING the two: 1664 comes from
## lri_cat, while 2189 and 3274 come from the any_* WHO flags.
cm <- ari %>% group_by(studyid) %>% summarise(m = min(lri_cat, na.rm = TRUE), .groups = "drop")
chk("Children, >=1 very severe episode (lri_cat)",
    n_distinct(ari$studyid[ari$lri_cat == 1]), 1664L)
chk("Children, >=1 severe episode (any_severe_lri)",
    n_distinct(ari$studyid[ari$any_severe_lri == 1]), 2189L)
chk("Children, >=1 non-severe episode (any_non_severe_lri)",
    n_distinct(ari$studyid[ari$any_non_severe_lri == 1]), 3274L)
say("NOTE two severity definitions disagree in ari.csv:")
say("  lri_cat==1 but v_severe_lri==0:      %d episodes",
    sum(ari$lri_cat == 1 & ari$v_severe_lri == 0, na.rm = TRUE))
say("  lri_cat==2 but severe_lri_who==0:    %d episodes",
    sum(ari$lri_cat == 2 & ari$severe_lri_who == 0, na.rm = TRUE))
say("  children very severe: %d by lri_cat vs %d by any_v_severe_lri (differ by %d)",
    n_distinct(ari$studyid[ari$lri_cat == 1]),
    n_distinct(ari$studyid[ari$any_v_severe_lri == 1]),
    n_distinct(ari$studyid[ari$lri_cat == 1]) - n_distinct(ari$studyid[ari$any_v_severe_lri == 1]))
say("  -> this analysis uses lri_cat throughout. Confirm the intended definition with Rowena.")

chk("Death location village/hospital/on-way",
    as.integer(table(dth$death_location)[c("village", "hospital", "on way")]), c(99L, 92L, 14L))
chk("Deaths by severity v.sev/sev/non-sev",
    as.integer(table(factor(dth$lri_cat, 1:3))), c(180L, 22L, 3L))
chk("Deaths <=1 day of first assessment",           sum(dth$died_lt_1d == 1, na.rm = TRUE), 149L)

say("")
say("=============== SURVEILLANCE PERIODS ===============")
print(table(ari$timeline, useNA = "ifany"))
say("Episode dates: %s to %s",
    format(min(as.Date(ari$ai_date_new, "%d%b%Y"), na.rm = TRUE)),
    format(max(as.Date(ari$ai_date_new, "%d%b%Y"), na.rm = TRUE)))
say("Death dates:   %s to %s",
    format(min(as.Date(dth$death_date_new, "%d%b%Y"), na.rm = TRUE)),
    format(max(as.Date(dth$death_date_new, "%d%b%Y"), na.rm = TRUE)))

say("")
say("=============== VILLAGE COVERAGE ===============")
v93 <- read_csv(data_path("93 village list.csv"), show_col_types = FALSE, progress = FALSE)
vc  <- read_dta(data_path("village_characteristics.dta"))
vnames <- names(v93)[1]
say("Villages in the 93-village list file: %d", nrow(v93))
say("Distinct villages appearing in ari.csv: %d", n_distinct(ari$code_village))
say("Villages in village_characteristics.dta: %d", nrow(vc))

## village_characteristics.dta keys on a bare 3-digit code ("006"); ari.csv uses "006-Tarubanda".
epi_v <- ari %>% mutate(vcode = substr(code_village, 1, 3)) %>% count(code_village, vcode, name = "episodes")
miss  <- epi_v$code_village[!epi_v$vcode %in% vc$village_code]
say("Episode-villages with NO characteristics row: %d", length(miss))
chk("Every village in ari.csv has a characteristics row", length(miss), 0L)

## The 12 villages on Rowena's "villages_missing_characteristics" list are NOT a gap
## inside this analysis: they sit OUTSIDE the 93-village set. Their children appear in
## demographics.csv but they contribute no episodes to ari.csv and no deaths to
## lrtideaths.csv, and they have no row in village_characteristics.dta.
dem_v <- dem %>% filter(!is.na(code_village_93)) %>%
  mutate(vcode = substr(as.character(code_village_93), 1, 3)) %>%
  count(code_village_93, vcode, name = "children")
outside <- dem_v %>% filter(!vcode %in% epi_v$vcode) %>% arrange(desc(children))
say("Villages in demographics.csv: %d; of these OUTSIDE the 93 analysis villages: %d",
    nrow(dem_v), nrow(outside))
if (nrow(outside)) {
  print(as.data.frame(outside %>% select(code_village_93, children)))
  say("Enrolled children in those villages: %d. Episodes contributed to ari.csv: %d. Deaths in lrtideaths.csv: %d",
      sum(outside$children),
      sum(ari$code_village %in% outside$code_village_93),
      sum(dth$code_village_93 %in% outside$code_village_93))
  say("  -> Rowena's spreadsheet lists 4 deaths in 050-Bhondilava; those are NOT among the 205.")
  say("  -> Ask whether these villages are correctly out of scope.")
}
write_csv(outside, file.path(OUT, "00_villages_outside_93.csv"))

say("")
say("=============== SEQUENCING FILE ===============")
say("Sequenced genomes: %d   RSV-A: %d   RSV-B: %d",
    nrow(seq), sum(seq$rsva == 1), sum(seq$rsvb == 1))
say("Distinct children: %d (2 sequenced twice = reinfections)", n_distinct(seq$studyid))
print(table(seq$year, seq$CLADE_f))
say("Ct range: %d-%d (sequencing was selected on low Ct)",
    min(seq$CTvalue_f, na.rm = TRUE), max(seq$CTvalue_f, na.rm = TRUE))

## RSV+ episodes in the cohort, to show what fraction was sequenced
rsvpos <- sum(ari$rsva_episode == 1 | ari$rsvb_episode == 1, na.rm = TRUE)
say("RSV-positive episodes in ari.csv: %d  ->  %.0f%% of them sequenced",
    rsvpos, 100 * nrow(seq) / rsvpos)

say("")
say("=============== KNOWN-BAD VARIABLES (do not use) ===============")
ga <- dem$gest_age_wk[!is.na(dem$gest_age_wk)]
say("demographics gest_age_wk: %.0f%% <37wk, %.0f%% exactly 36 -> parental recall, unusable",
    100 * mean(ga < 37), 100 * mean(ga == 36))
say("deaths gest_age_lmp_wk present for %d/%d deaths -> use this one instead",
    sum(!is.na(dth$gest_age_lmp_wk)), nrow(dth))
say("referred_hosp_episode flagged on %d episodes but only %.1f%% of them reached a facility",
    sum(ari$referred_hosp_episode == 1, na.rm = TRUE),
    100 * mean(ari$community[ari$referred_hosp_episode == 1] == 0, na.rm = TRUE))
say("  -> NOT a referral-completion variable. Confirm definition with Rowena.")

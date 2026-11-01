# 26_ascertainment.R -- who does the genomic dataset represent, and who is
# under-represented?
#
# Three nested populations:
#   A  all LRTI episodes in the cohort
#   B  tested (swabbed) episodes, a subset of A
#   B+ RSV-positive tested episodes, a subset of B
#   C  sequenced AND linked episodes, a subset of B+
#
# Two selection steps are therefore measurable from these data, and each is
# compared between DISJOINT groups. An earlier version compared the selected
# group against its own parent population, which contains it and therefore
# dilutes the contrast:
#   tested      vs untested                            (disjoint)
#   sequenced   vs unsequenced RSV-positive tested      (disjoint)
#
# Reported as STANDARDISED DIFFERENCES, not significance tests. With ~9,000
# episodes a trivial difference reaches p < 0.05, and with 203 sequenced
# genomes a real one may not; neither tells you about representativeness.
# |SMD| >= 0.10 is flagged as a difference worth noting. Missingness is reported
# per variable per population, because a variable absent from one group cannot
# be balanced on.
#
# This measures MEASURED selection only. No weighting repairs unmeasured
# selection, and none is attempted here.

source(file.path("scripts", "00_setup.R"))
cat("\n== 26 ascertainment and sequencing bias ==\n")

ari <- read_csv(data_path("ari.csv"), show_col_types = FALSE, guess_max = 20000) %>%
  mutate(ep_date = parse_dmy(ai_date_new), season = rsv_season(ep_date),
         swabbed = !is.na(nswab_episode) & nswab_episode == 1,
         rsv_pos = (!is.na(rsva_episode) & rsva_episode == 1) |
                   (!is.na(rsvb_episode) & rsvb_episode == 1),
         very_severe = as.integer(lri_cat == 1),
         female = as.integer(sexr == 2),
         facility = as.integer(hosp_lri == 1),
         age_months = age_days_episode / 30.4)

dat <- readRDS(file.path(OUT, "02_dat.rds"))
seq_keys <- dat %>% filter(link_status %in% c("exact", "nearest_prior_episode_same_RSV_type")) %>%
  distinct(studyid, ep_date = as.Date(ep_date))
ari <- ari %>% mutate(sequenced = paste(studyid, ep_date) %in%
                        paste(seq_keys$studyid, seq_keys$ep_date))

A  <- ari
B  <- ari %>% filter(swabbed)
Bp <- ari %>% filter(swabbed, rsv_pos)
## DECLARED DEFINITION (see 27_cohort_status_audit.R): RSV-positive means the
## COHORT flags say so. Two linked episodes are sequenced but cohort-negative;
## they are excluded here so numerator and denominator share one definition, and
## retained in the genomic analyses where the laboratory result is the relevant
## one. Using all 203 linked episodes here would compare a 203-episode numerator
## against a 644-episode denominator built on a different rule.
Cc <- ari %>% filter(swabbed, rsv_pos, sequenced)
## disjoint complements
Bn <- ari %>% filter(!swabbed)                       # untested
Cn <- Bp %>% filter(!sequenced)                      # RSV-positive tested, not sequenced

cat(sprintf("\n  A  all LRTI episodes          %5d\n", nrow(A)))
cat(sprintf("  B  tested (swabbed)           %5d  (%.1f%% of A)\n", nrow(B), 100*nrow(B)/nrow(A)))
cat(sprintf("  B+ RSV-positive and tested    %5d  (%.1f%% of B)\n", nrow(Bp), 100*nrow(Bp)/nrow(B)))
cat(sprintf("  C  sequenced, linked AND cohort-positive %5d  (%.1f%% of B+)\n",
            nrow(Cc), 100*nrow(Cc)/nrow(Bp)))
cat(sprintf("     (%d further linked episodes are cohort-negative; see 27_cohort_status_audit.R)\n",
            sum(ari$sequenced) - nrow(Cc)))

smd_bin <- function(x, y) {
  p1 <- mean(x, na.rm = TRUE); p2 <- mean(y, na.rm = TRUE)
  s <- sqrt((p1*(1-p1) + p2*(1-p2))/2)
  if (is.na(s) || s == 0) return(NA_real_); (p1 - p2)/s
}
smd_num <- function(x, y) {
  s <- sqrt((var(x, na.rm = TRUE) + var(y, na.rm = TRUE))/2)
  if (is.na(s) || s == 0) return(NA_real_)
  (mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE))/s
}

compare <- function(g1, g2, l1, l2) {
  rows <- list()
  num <- c(age_months = "age_months")
  for (v in names(num)) rows[[v]] <- tibble(
    variable = v, type = "mean",
    val1 = round(mean(g1[[v]], na.rm = TRUE), 2), val2 = round(mean(g2[[v]], na.rm = TRUE), 2),
    smd = round(smd_num(g1[[v]], g2[[v]]), 3),
    miss1 = round(100*mean(is.na(g1[[v]])), 1), miss2 = round(100*mean(is.na(g2[[v]])), 1))
  bin <- c(very_severe = "very_severe", female = "female", facility = "facility")
  for (v in names(bin)) rows[[v]] <- tibble(
    variable = v, type = "%",
    val1 = round(100*mean(g1[[v]], na.rm = TRUE), 1), val2 = round(100*mean(g2[[v]], na.rm = TRUE), 1),
    smd = round(smd_bin(g1[[v]], g2[[v]]), 3),
    miss1 = round(100*mean(is.na(g1[[v]])), 1), miss2 = round(100*mean(is.na(g2[[v]])), 1))
  for (s in sort(unique(na.omit(A$season)))) rows[[paste0("season ", s)]] <- tibble(
    variable = paste0("season ", s), type = "%",
    val1 = round(100*mean(g1$season == s, na.rm = TRUE), 1),
    val2 = round(100*mean(g2$season == s, na.rm = TRUE), 1),
    smd = round(smd_bin(g1$season == s, g2$season == s), 3),
    miss1 = round(100*mean(is.na(g1$season)), 1), miss2 = round(100*mean(is.na(g2$season)), 1))
  out <- bind_rows(rows) %>% rename(!!l1 := val1, !!l2 := val2,
                                    !!paste0("miss_", l1) := miss1, !!paste0("miss_", l2) := miss2)
  out
}

cat(sprintf("\n  disjoint groups: untested %d | RSV-positive tested but not sequenced %d\n",
            nrow(Bn), nrow(Cn)))

cat("\n-- Step 1: tested vs UNTESTED (disjoint) --\n")
S1 <- compare(B, Bn, "tested", "untested")
print(as.data.frame(S1), right = FALSE)

cat("\n-- Step 2: sequenced vs UNSEQUENCED RSV-positive tested (disjoint) --\n")
S2 <- compare(Cc, Cn, "sequenced", "not_sequenced")
print(as.data.frame(S2), right = FALSE)

## within-season, for the sequencing step, where both groups have >= 20
cat("\n-- Step 2 within season (both groups >= 20 episodes) --\n")
ws <- list()
for (s in sort(unique(na.omit(Bp$season)))) {
  g1 <- Cc %>% filter(season == s); g2 <- Cn %>% filter(season == s)
  if (nrow(g1) < 20 || nrow(g2) < 20) {
    cat(sprintf("  %s: sequenced %3d / not sequenced %3d -> not compared (needs >=20 in both)\n",
                s, nrow(g1), nrow(g2))); next }
  r <- tibble(season = s, n_seq = nrow(g1), n_unseq = nrow(g2),
              smd_very_severe = round(smd_bin(g1$very_severe, g2$very_severe), 3),
              smd_facility    = round(smd_bin(g1$facility, g2$facility), 3),
              smd_age         = round(smd_num(g1$age_months, g2$age_months), 3))
  ws[[s]] <- r
  cat(sprintf("  %s: sequenced %3d / not %3d | SMD severe %+0.3f, facility %+0.3f, age %+0.3f\n",
              s, nrow(g1), nrow(g2), r$smd_very_severe, r$smd_facility, r$smd_age))
}
WS <- bind_rows(ws)
if (nrow(WS)) write_csv(WS, file.path(OUT, "26_sequencing_step_within_season.csv"))

## Season indicators are levels of ONE categorical variable, not independent
## characteristics, so they are named rather than counted alongside the others.
flag <- function(S, lab) {
  b  <- S %>% filter(!is.na(smd), abs(smd) >= 0.10)
  bc <- b %>% filter(!grepl("^season ", variable))
  bs <- b %>% filter(grepl("^season ", variable))
  cat(sprintf("\n  %s, imbalances with |SMD| >= 0.10:\n", lab))
  if (nrow(bc)) for (i in seq_len(nrow(bc))) cat(sprintf("    %-18s SMD %+0.3f\n", bc$variable[i], bc$smd[i]))
  else cat("    none among age, sex, severity, facility attendance\n")
  if (nrow(bs)) cat(sprintf("    season composition (one variable): %s\n",
                            paste(sprintf("%s %+0.2f", sub("season ", "", bs$variable), bs$smd), collapse = ", ")))
  else cat("    season composition: no level reaches 0.10\n")
}
flag(S1, "testing step"); flag(S2, "sequencing step")

## village-level coverage: which villages are under-represented?
cat("\n-- village coverage (episodes sequenced / RSV-positive tested) --\n")
vc <- Bp %>% count(code_village, name = "rsv_pos") %>%
  left_join(Cc %>% count(code_village, name = "sequenced"), by = "code_village") %>%
  mutate(sequenced = coalesce(sequenced, 0L), pct = round(100*sequenced/rsv_pos, 1)) %>%
  arrange(desc(rsv_pos))
cat(sprintf("  %d villages contribute RSV-positive episodes; %d contribute none sequenced\n",
            nrow(vc), sum(vc$sequenced == 0)))
cat(sprintf("  coverage across villages: median %.1f%%, IQR %.1f-%.1f%%, range %.1f-%.1f%%\n",
            median(vc$pct), quantile(vc$pct,.25), quantile(vc$pct,.75), min(vc$pct), max(vc$pct)))
print(as.data.frame(head(vc, 10)), right = FALSE)

## ---- adherence to the stated swabbing rule -------------------------------
## The cohort protocol (companion access-to-care paper, S1.4) directs a swab for
## any child assessed with severe or very severe LRTI. That is a DESIGN, and it
## explains the direction of the testing-step imbalance -- but it is not a
## complete account of who was tested, in either direction:
##   - not every severe episode was swabbed (incomplete implementation), and
##   - some non-severe episodes were swabbed anyway.
## Both are quantified here so the manuscript can quote them from a table rather
## than assert that the protocol fully explains the sample.
rule <- tibble::tibble(
  quantity = c("severe or very severe, among tested",
               "tested, among severe or very severe",
               "non-severe among tested",
               "non-severe among tested, recorded community-managed",
               "non-severe among tested, recorded facility"),
  numerator = c(
    sum(ari$swabbed & ari$lri_cat %in% c(1, 2), na.rm = TRUE),
    sum(ari$lri_cat %in% c(1, 2) & ari$swabbed, na.rm = TRUE),
    sum(ari$swabbed & !(ari$lri_cat %in% c(1, 2)), na.rm = TRUE),
    sum(ari$swabbed & !(ari$lri_cat %in% c(1, 2)) & ari$community == 1, na.rm = TRUE),
    sum(ari$swabbed & !(ari$lri_cat %in% c(1, 2)) & ari$hosp_lri == 1, na.rm = TRUE)),
  denominator = c(
    sum(ari$swabbed, na.rm = TRUE),
    sum(ari$lri_cat %in% c(1, 2), na.rm = TRUE),
    sum(ari$swabbed, na.rm = TRUE),
    sum(ari$swabbed & !(ari$lri_cat %in% c(1, 2)), na.rm = TRUE),
    sum(ari$swabbed & !(ari$lri_cat %in% c(1, 2)), na.rm = TRUE))) %>%
  mutate(pct = round(100 * numerator / denominator, 1))
cat("\n-- adherence to the stated swabbing rule --\n")
print(as.data.frame(rule))
cat("  The rule is severity-triggered, so it explains the DIRECTION of the\n")
cat("  testing-step imbalance. It does not fully explain WHO was tested:\n")
cat("  severe episodes were under-sampled and some non-severe were sampled.\n")
write_csv(rule, file.path(OUT, "26_testing_rule_adherence.csv"))

write_csv(S1, file.path(OUT, "26_ascertainment_testing_step.csv"))
write_csv(S2, file.path(OUT, "26_ascertainment_sequencing_step.csv"))
write_csv(vc, file.path(OUT, "26_village_coverage.csv"))
cat("\n  -> wrote 26_ascertainment_testing_step.csv, 26_ascertainment_sequencing_step.csv,\n")
cat("     26_village_coverage.csv, 26_testing_rule_adherence.csv\n")
cat("== 26 done ==\n")

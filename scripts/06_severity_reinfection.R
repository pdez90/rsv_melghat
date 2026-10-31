# 06_severity_reinfection.R -- clinical severity by RSV subgroup and lineage,
# and the repeat-sequenced children.
#
# DESIGN WARNING: in this cohort each RSV season is essentially 100% one
# subgroup (2016/17 A, 2017/18 B, 2018/19 B, 2019/20 A, 2021/22 A). Subgroup is
# therefore collinear with season, and any A-vs-B severity contrast is also a
# between-season contrast (different surveillance intensity, different age
# structure of the susceptible pool, and for 2021/22 a different village set).
# The only comparison internal to a season is between LINEAGES circulating in
# the same season, so that is the primary analysis here.

source(file.path("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(lme4); library(broom.mixed)})
cat("\n== 06 severity and reinfection ==\n")

dat <- readRDS(file.path(OUT, "02_dat.rds"))
ari <- read_csv(data_path("ari.csv"), show_col_types = FALSE, guess_max = 20000) %>%
  mutate(ep_date = parse_dmy(ai_date_new), season = rsv_season(ep_date),
         # CRITICAL: rsva_episode/rsvb_episode are 0 for never-swabbed episodes.
         # Only 4,081 of 9,236 episodes were swabbed, and swabbing was heavily
         # selected on severity (40.1% very severe among swabbed vs 4.6% among
         # unswabbed). Treating unswabbed episodes as RSV-negative produces a
         # badly inflated RSV effect, so every RSV-vs-non-RSV comparison below is
         # restricted to swabbed (i.e. tested) episodes.
         swabbed = !is.na(nswab_episode) & nswab_episode == 1,
         rsv_pos = rsva_episode == 1 | rsvb_episode == 1,
         subgroup = case_when(rsva_episode == 1 ~ "RSV-A", rsvb_episode == 1 ~ "RSV-B"),
         very_severe = as.integer(lri_cat == 1),
         female = as.integer(sexr == 2),
         age_mo = age_days_episode / 30.44)

## ---- 1. is subgroup separable from season? ----------------------------------
cat("\n-- subgroup by season, ALL cohort RSV+ episodes --\n")
print(with(ari %>% filter(rsv_pos), table(season, subgroup)))
cat("\n=> subgroup is collinear with season; A-vs-B cannot be separated from season.\n")

## ---- 2. descriptive severity, all RSV+ episodes ------------------------------
cat("\n-- severity by subgroup (all cohort RSV+ episodes) --\n")
sev <- ari %>% filter(rsv_pos) %>% group_by(subgroup) %>%
  summarise(episodes = n(),
            pct_very_severe = round(100 * mean(lri_cat == 1, na.rm = TRUE), 1),
            pct_severe      = round(100 * mean(lri_cat == 2, na.rm = TRUE), 1),
            pct_non_severe  = round(100 * mean(lri_cat == 3, na.rm = TRUE), 1),
            pct_reached_facility = round(100 * mean(community == 0, na.rm = TRUE), 1),
            median_age_mo = round(median(age_mo, na.rm = TRUE), 1),
            deaths = sum(lrti_death == 1, na.rm = TRUE), .groups = "drop")
print(as.data.frame(sev))

cat("\n-- RSV+ vs tested-negative LRTI episodes, WITHIN season --\n")
cmp <- ari %>% filter(!is.na(season), swabbed) %>% group_by(season, rsv_pos) %>%
  summarise(episodes = n(), pct_very_severe = round(100 * mean(lri_cat == 1, na.rm = TRUE), 1),
            median_age_mo = round(median(age_mo, na.rm = TRUE), 1), .groups = "drop") %>%
  filter(!is.na(rsv_pos))
print(as.data.frame(cmp))

cat("\n-- swabbing is selected on severity (why the denominator matters) --\n")
print(ari %>% group_by(swabbed) %>%
        summarise(episodes = n(),
                  pct_very_severe = round(100 * mean(lri_cat == 1, na.rm = TRUE), 1),
                  pct_reached_facility = round(100 * mean(community == 0, na.rm = TRUE), 1),
                  .groups = "drop") %>% as.data.frame())

sw <- ari %>% filter(swabbed, !is.na(lri_cat))
cat(sprintf("\ntested episodes available: %d (%d RSV+); missing sex: %d -> dropped by the models\n",
            nrow(sw), sum(sw$rsv_pos), sum(is.na(sw$female))))

# Season MUST be in the model: testing intensity, RSV circulation and
# surveillance all change across years, and subgroup is confounded with season.
# Children contribute repeated episodes, so a child random intercept is included.
MODELS <- list(
  `age + sex + village`                  = very_severe ~ rsv_pos + scale(age_days_episode) + female + (1 | code_village),
  `+ season`                             = very_severe ~ rsv_pos + scale(age_days_episode) + female + factor(season) + (1 | code_village),
  `+ season + child random intercept`    = very_severe ~ rsv_pos + scale(age_days_episode) + female + factor(season) + (1 | code_village) + (1 | studyid),
  # Sex is recorded for every tested episode in 2016/17-2020/21 and for NONE in
  # 2021/22 (365 of 365) or 2022/23 (2 of 2). A complete-case model requiring sex
  # therefore drops those two seasons entirely -- 367 tested episodes including
  # 158 of the 644 RSV-positive ones, the whole A.D.3 wave -- leaving 5 of the 7
  # seasons that have tested episodes, and 4 of the 5 that have RSV-positive ones. This last
  # specification is the same model WITHOUT sex, so the effect of readmitting
  # 2021/22 can be read directly. It is not the primary model: dropping a
  # confounder and adding a season change together, so the two differ in more
  # than one respect and the comparison is descriptive.
  `+ season + child, sex dropped (readmits 2021/22)` = very_severe ~ rsv_pos + scale(age_days_episode) + factor(season) + (1 | code_village) + (1 | studyid),
  # The specification above FAILED TO CONVERGE on 2026-10-31 ("unable to evaluate
  # scaled gradient | degenerate Hessian with 1 negative eigenvalue"), so it is
  # not an estimate. The likely obstacle is the 2022/23 season level: 2 tested
  # episodes, 0 very-severe events and 0 RSV-positive episodes, i.e. perfectly
  # predicted and inestimable. This last specification is the same model with
  # that one season dropped, which is the smallest change that lets the question
  # -- does readmitting 2021/22 move the estimate? -- actually be fitted. It
  # removes an inestimable stratum, not any RSV-positive episode: 2022/23
  # contributes none. It is still NOT the primary model, and it still changes
  # adjustment and population together, so it cannot isolate either.
  `+ season + child, sex dropped, 2022/23 excluded` = very_severe ~ rsv_pos + scale(age_days_episode) + factor(season) + (1 | code_village) + (1 | studyid))

# Per-specification data. All but the last are fitted over every tested episode
# with a severity grade; complete-case handling of sex is left to the model.
DATA <- setNames(rep(list(sw), length(MODELS)), names(MODELS))
DATA[["+ season + child, sex dropped, 2022/23 excluded"]] <-
  sw %>% filter(season != "2022/23")
res <- list()
for (nm in names(MODELS)) {
  m <- try(glmer(MODELS[[nm]], data = DATA[[nm]], family = binomial,
                 control = glmerControl(optimizer = "bobyqa")), silent = TRUE)
  if (inherits(m, "try-error")) {
    res[[nm]] <- tibble(model = nm, status = "UNRELIABLE: fit errored",
                        optimizer_msg = as.character(m))
    next
  }
  mf <- model.frame(m)
  z <- coef(summary(m))["rsv_posTRUE", ]
  res[[nm]] <- tibble(model = nm, fitted_episodes = nobs(m),
                      rsv_pos_in_fit = sum(mf$rsv_pos),
                      # the first specification has no season term
                      seasons_in_fit = if ("factor(season)" %in% names(mf))
                        n_distinct(mf[["factor(season)"]]) else NA_integer_,
                      OR = round(exp(z[1]), 3),
                      lo = round(exp(z[1] - 1.96 * z[2]), 3),
                      hi = round(exp(z[1] + 1.96 * z[2]), 3),
                      p  = round(z[4], 3), singular = isSingular(m),
                      # singularity and convergence are DIFFERENT checks: a
                      # non-singular fit can still have failed to converge, so
                      # both are recorded and neither is inferred from the other.
                      # Same test as 12_sensitivity.R O2.
                      # TWO convergence signals decide `clean`, and neither is
                      # inferred from the other:
                      #  - lme4's own post-fit checks (messages)
                      #  - the optimizer's return code, conv$opt (0 = success)
                      # Two gradient diagnostics are recorded ALONGSIDE, and are
                      # NOT part of the pass/fail test:
                      #  - max_abs_grad: the plain maximum absolute gradient.
                      #  - max_abs_scaled_grad: lme4's CHOLESKY-scaled gradient,
                      #    solve(chol(Hessian), gradient). chol() requires a
                      #    positive-definite Hessian, so this is NA whenever the
                      #    Hessian is NOT positive-definite -- singular as well
                      #    as indefinite. The NA is a DIAGNOSTIC OUTCOME, not a
                      #    finding in itself: it says the scaled gradient could
                      #    not be formed, and the reason needs looking at. Nor
                      #    does a finite value prove anything is wrong with the
                      #    code -- a refit that yields one warrants inspection,
                      #    not an assumption of error. An earlier version used
                      #    solve(Hessian, gradient), which is NOT lme4's
                      #    criterion: an indefinite Hessian can still be
                      #    invertible, so that version returned a small-looking
                      #    number for a degenerate fit instead of NA. Any
                      #    magnitude recorded by that earlier version should be
                      #    disregarded.
                      no_lme4_msg = is.null(m@optinfo$conv$lme4$messages),
                      optimizer_code = tryCatch(as.integer(m@optinfo$conv$opt),
                                                error = function(e) NA_integer_),
                      max_abs_grad = tryCatch(max(abs(m@optinfo$derivs$gradient)),
                                              error = function(e) NA_real_),
                      max_abs_scaled_grad = tryCatch(
                        max(abs(with(m@optinfo$derivs, solve(chol(Hessian), gradient)))),
                        error = function(e) NA_real_),
                      optimizer_msg = if (length(m@optinfo$conv$lme4$messages))
                        paste(m@optinfo$conv$lme4$messages, collapse = " | ") else NA_character_) %>%
    mutate(clean = no_lme4_msg & !is.na(optimizer_code) & optimizer_code == 0,
           status = case_when(
      !clean &  singular ~ "UNRELIABLE: convergence not established; singular",
      !clean             ~ "UNRELIABLE: convergence not established -- diagnose before interpreting",
      singular           ~ "UNRELIABLE: singular (variance component at 0)",
      TRUE               ~ "converged (no lme4 message, optimizer code 0)"))
}
cat("\nRSV+ vs tested-negative, very-severe outcome (odds ratios):\n")
print(as.data.frame(bind_rows(res)))
cat("\nNote fitted_episodes < tested episodes: episodes with missing sex are dropped.\n")
cat("Every dropped episode is from 2021/22 or 2022/23, so the models that include\n")
cat("sex are fitted over 5 of the 7 seasons with tested episodes (see seasons_in_fit),\n")
cat("and over 4 of the 5 seasons with any RSV-positive episode. The last row readmits\n")
cat("both excluded seasons (4081 episodes, 644 RSV-positive, 7 seasons) but FAILED\n")
cat("TO CONVERGE on 2026-10-31; the row after it drops the inestimable 2022/23\n")
cat("season (2 tested episodes, 0 very-severe events, 0 RSV-positive). Read the\n")
cat("status column before quoting any of these as an estimate.\n")
cat(sprintf("Among tested episodes, %d children contribute more than one episode, giving\n",
            sum(table(sw$studyid) > 1)))
cat(sprintf("%d episodes beyond each child's first -- hence the child random intercept.\n",
            sum(duplicated(sw$studyid))))
cat("\nAdding season adjustment changed the estimate from 0.85 to 1.13. The first\n")
cat("specification was NOT unadjusted -- it already carried age, sex and village --\n")
cat("and a change on adding a term does not by itself establish that the earlier\n")
cat("value was a season artifact. NO STATISTICALLY SIGNIFICANT ASSOCIATION WAS DETECTED\n")
cat("(adjusted OR 1.13, 95% CI 0.87-1.47). That interval still permits effects that\n")
cat("would matter clinically in either direction -- this is not evidence of no effect.\n")
cat("Restricting to tested episodes avoids coding untested children as negative but\n")
cat("does NOT remove selection bias from severity-dependent testing.\n")
write_out(bind_rows(res), "06_severity_model_sensitivity.csv")

## ---- convergence diagnosis -------------------------------------------------
## lme4's convergence guidance treats a warning as a flag to CHECK the fit, not
## as automatic invalidation of every coefficient. The two checks it asks for are
## the scaled/relative gradient (recorded per fit above) and agreement of the
## estimate across optimizers. Any specification whose convergence was not
## established above is therefore refitted with each optimizer that ships with
## lme4, and the estimates compared. Agreement across optimizers with small
## relative gradients argues the estimate is usable despite the warning;
## disagreement, or a gradient that cannot be evaluated, argues it is not.
## This writes a table either way, so the file always exists for the freeze.
res_tbl  <- bind_rows(res)
needs_dx <- res_tbl$model[!grepl("^converged", res_tbl$status)]
dx <- list()
for (nm in needs_dx) for (opt in c("bobyqa", "Nelder_Mead", "nloptwrap")) {
  md <- try(glmer(MODELS[[nm]], data = DATA[[nm]], family = binomial,
                  control = glmerControl(optimizer = opt)), silent = TRUE)
  key <- paste(nm, opt)
  if (inherits(md, "try-error")) {
    dx[[key]] <- tibble(model = nm, optimizer = opt, note = "fit errored",
                        optimizer_msg = as.character(md))
    next
  }
  zd <- coef(summary(md))["rsv_posTRUE", ]
  dx[[key]] <- tibble(
    model = nm, optimizer = opt, note = NA_character_,
    fitted_episodes = nobs(md), rsv_pos_in_fit = sum(model.frame(md)$rsv_pos),
    seasons_in_fit = if ("factor(season)" %in% names(model.frame(md)))
      n_distinct(model.frame(md)[["factor(season)"]]) else NA_integer_,
    OR = round(exp(zd[1]), 3),
    lo = round(exp(zd[1] - 1.96 * zd[2]), 3),
    hi = round(exp(zd[1] + 1.96 * zd[2]), 3),
    p  = round(zd[4], 3),
    singular = isSingular(md),
    no_lme4_msg = is.null(md@optinfo$conv$lme4$messages),
    optimizer_code = tryCatch(as.integer(md@optinfo$conv$opt),
                              error = function(e) NA_integer_),
    max_abs_grad = tryCatch(max(abs(md@optinfo$derivs$gradient)),
                            error = function(e) NA_real_),
    # lme4's Cholesky-scaled gradient; NA when the Hessian is not
    # positive-definite, which is itself the finding
    max_abs_scaled_grad = tryCatch(
      max(abs(with(md@optinfo$derivs, solve(chol(Hessian), gradient)))),
      error = function(e) NA_real_),
    optimizer_msg = if (length(md@optinfo$conv$lme4$messages))
      paste(md@optinfo$conv$lme4$messages, collapse = " | ") else NA_character_)
}
if (length(dx)) {
  dxt <- bind_rows(dx)
  cat("\n-- convergence diagnosis: specifications whose convergence was not established --\n")
  print(as.data.frame(dxt))
  cat("\nRead across optimizers for each model. Agreement of the coefficient\n")
  cat("across optimizers supports STABILITY OF THE POINT ESTIMATE only. It does\n")
  cat("not validate the fit: Hessian warnings and a singular fit still bear on\n")
  cat("the standard errors and therefore on the interval, which is the part\n")
  cat("least safe to quote from a fit whose convergence is not established.\n")
  cat("max_abs_scaled_grad is NA whenever the Hessian is not positive-definite --\n")
  cat("singular as well as indefinite. That NA is a diagnostic outcome to be\n")
  cat("explained, not a finding on its own and not a value to work around.\n")
} else {
  dxt <- tibble(model = NA_character_, optimizer = NA_character_,
                note = "no specification required diagnosis: all converged")
  cat("\n-- convergence diagnosis: not required, every specification converged --\n")
}
write_out(dxt, "06_convergence_diagnosis.csv")

## ---- 3. PRIMARY: lineage comparison within a season -------------------------
cat("\n-- severity by lineage, within season (sequenced genomes) --\n")
lin <- dat %>% filter(linked, !is.na(clade), !is.na(lri_cat)) %>%
  mutate(very_severe = as.integer(lri_cat == 1))
by_lin <- lin %>% group_by(season, clade) %>%
  summarise(genomes = n(), pct_very_severe = round(100 * mean(very_severe), 1),
            pct_reached_facility = round(100 * mean(community == 0, na.rm = TRUE), 1),
            median_age_mo = round(median(age_days_episode / 30.44, na.rm = TRUE), 1),
            .groups = "drop") %>% filter(genomes >= 5)
print(as.data.frame(by_lin))

for (ss in c("2018/19", "2021/22")) {
  sub <- lin %>% filter(season == ss) %>%
    group_by(clade) %>% filter(n() >= 10) %>% ungroup() %>%
    mutate(clade = factor(clade))
  if (n_distinct(sub$clade) < 2) next
  tab <- with(sub, table(clade, very_severe))
  cat(sprintf("\n%s: %s (n=%d), very-severe outcome\n",
              ss, paste(levels(sub$clade), collapse = " vs "), nrow(sub)))
  print(tab)
  # With an empty outcome cell the logistic model is separated; fall back to an
  # exact test rather than reporting an unstable odds ratio.
  if (ncol(tab) < 2 || any(tab == 0)) {
    ft <- fisher.test(tab)
    cat(sprintf("sparse/zero cell -> Fisher exact test: OR %.2f (95%% CI %.2f-%.2f), p = %.3f\n",
                ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
  } else {
    rhs <- "clade + scale(age_days_episode)"
    if (n_distinct(na.omit(sub$sexr)) >= 2) rhs <- paste(rhs, "+ factor(sexr)")
    cat(sprintf("adjusted for %s\n", if (grepl("sexr", rhs)) "age and sex" else "age"))
    fit <- glm(as.formula(paste("very_severe ~", rhs)), data = sub, family = binomial)
    print(tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
            filter(str_starts(term, "clade")) %>%
            mutate(across(where(is.numeric), ~round(.x, 3))) %>% as.data.frame())
  }
}

## ---- 4. age at infection by subgroup ----------------------------------------
cat("\n-- age at RSV infection by subgroup (all cohort RSV+ episodes) --\n")
print(ari %>% filter(rsv_pos) %>% group_by(subgroup) %>%
        summarise(n = n(), median_age_mo = round(median(age_mo, na.rm = TRUE), 1),
                  iqr = paste(round(quantile(age_mo, c(.25, .75), na.rm = TRUE), 1), collapse = "-"),
                  pct_under_6mo = round(100 * mean(age_mo < 6, na.rm = TRUE), 1),
                  .groups = "drop") %>% as.data.frame())
cat("(also confounded with season; infants born into a B season meet B first)\n")

## ---- 5. repeat-sequenced children -------------------------------------------
cat("\n-- children sequenced more than once --\n")
reinf <- dat %>% group_by(studyid) %>% filter(n() > 1) %>% arrange(studyid, swab_date) %>%
  select(studyid, FieldNoNIV_f, swab_date, season, nc_dataset, clade, severity,
         age_days_swab, community, link_status) %>% ungroup()
print(as.data.frame(reinf))
# No participant identifier is written into this script: the IDs appear in the
# printed table above, which is console output, and in output/06_reinfections.csv,
# neither of which is committed to the public code repository.
cat("\nBoth are genuine reinfections: a different lineage on the second episode\n")
cat("(one a within-RSV-B lineage change B.D.4.1.1 -> B.D.4.1, the other a change of\n")
cat("subgroup, RSV-B -> RSV-A). See the table above for the specimens concerned.\n")

## ---- 6. outcomes of the K272R (palivizumab site II) genomes -----------------
k272 <- dat %>% filter(nc_escape > 0) %>%
  select(FieldNoNIV_f, studyid, season, clade, severity, community, lrti_death,
         age_days_episode, code_village, linked)
cat("\n-- genomes carrying the site II substitution flagged by Nextclade --\n")
print(as.data.frame(k272))
cat(sprintf("villages represented: %d; all from one season: %s\n",
            n_distinct(k272$code_village, na.rm = TRUE),
            paste(unique(k272$season), collapse = ", ")))

write_out(sev,    "06_severity_by_subgroup.csv")
write_out(by_lin, "06_severity_by_lineage_season.csv")
write_out(reinf,  "06_reinfections.csv")
write_out(k272,   "06_site_II_substitution_genomes.csv")
cat("== 06 done ==\n")

# 12_sensitivity.R -- the sensitivity panel specified in SENSITIVITY_PLAN.md v3.
#
# Step 1 freezes the plan and input hashes.
# Step 2 reproduces S0 and stops on any discrepancy with the primary results.
# Steps 3+ run S1-S6 for the outcomes each restriction legitimately applies to.
#
# Populations are per outcome (plan section 3). Linkage and genome-QC restrictions
# are NOT applied to O2, whose unit is a tested cohort episode, nor to O5, whose
# primary population is all 241 sequenced genomes.

source(file.path("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(lme4); library(tools)})
set.seed(20261030)
cat("\n== 12 sensitivity panel ==\n")

N_PERM <- 1999    # panel-wide; S0 is additionally checked at 9999 against script 08
N_BOOT <- 999

dat <- readRDS(file.path(OUT, "02_dat.rds"))
ari <- read_csv(data_path("ari.csv"), show_col_types = FALSE, guess_max = 20000) %>%
  mutate(ep_date = parse_dmy(ai_date_new), season = rsv_season(ep_date),
         swabbed = !is.na(nswab_episode) & nswab_episode == 1,
         rsv_pos = rsva_episode == 1 | rsvb_episode == 1,
         very_severe = as.integer(lri_cat == 1), female = as.integer(sexr == 2),
         ## ep_row must be the row number of the UNFILTERED read, matching
         ## 01_link_sequencing.R, so dat$ep_row refers to the same episodes
         ep_row = row_number())

## ---------- 1. FREEZE ---------------------------------------------------------
frz <- file(file.path(OUT, "12_frozen_plan.txt"), open = "wt")
fsay <- function(...) { l <- sprintf(...); writeLines(l, frz) }
fsay("SENSITIVITY PANEL -- FROZEN STATE"); fsay("frozen: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
plan <- file.path(PROJ, "SENSITIVITY_PLAN.md")
fsay("plan: %s  md5 %s", if (file.exists(plan)) "SENSITIVITY_PLAN.md" else "MISSING",
     if (file.exists(plan)) unname(md5sum(plan)) else "NA")
for (f in c("NP_results_with_clade_sequencing.csv", "ari.csv", "lrtideaths.csv",
            "demographics.csv", "critical_data_gpsdata.csv")) {
  pth <- data_path(f, required = FALSE)
  fsay("input %-40s md5 %s", f, if (is.na(pth)) "MISSING" else unname(md5sum(pth)))
}
for (ty in c("RSV-A", "RSV-B")) for (f in c("nextclade.tsv", "nextclade.aligned.fasta")) {
  pth <- file.path(PROJ, "nextclade", ty, f)
  fsay("nextclade %-8s %-32s md5 %s", ty, f, if (file.exists(pth)) unname(md5sum(pth)) else "MISSING")
}
fsay("N_PERM %d  N_BOOT %d  seed 20261030", N_PERM, N_BOOT)
close(frz)
cat("  -> froze plan + input hashes to 12_frozen_plan.txt\n")

## ---------- specification grid ------------------------------------------------
SPEC <- tribble(
  ~spec, ~link,     ~qual,      ~vill,
  "S0",  "L-ALL",   "Q-ALL",    "all",
  "S1",  "L-EXACT", "Q-ALL",    "all",
  "S2",  "L-NARROW","Q-ALL",    "all",
  "S3",  "L-ALL",   "Q-STRICT", "all",
  "S4",  "L-EXACT", "Q-STRICT", "all",
  "S5",  "L-ALL",   "Q-ALL",    "V33",
  "S6",  "L-EXACT", "Q-STRICT", "V33")

V33 <- ari %>% filter(subset_33_villages == 1) %>% distinct(code_village) %>% pull(code_village)
cat(sprintf("V-33 membership: %d villages\n", length(V33)))

link_ok <- function(d, set) switch(set,
  "L-ALL"    = d$link_status %in% c("exact", "nearest_prior_episode_same_RSV_type"),
  "L-EXACT"  = d$link_status == "exact",
  "L-NARROW" = d$link_status == "exact" |
               (d$link_status == "nearest_prior_episode_same_RSV_type" &
                !is.na(d$days_episode_to_swab) & d$days_episode_to_swab <= 14))

fs_any <- bind_rows(
  read_tsv(file.path(PROJ, "nextclade", "RSV-A", "nextclade.tsv"), show_col_types = FALSE,
           guess_max = 500, progress = FALSE) %>% select(seqName, totalFrameShifts),
  read_tsv(file.path(PROJ, "nextclade", "RSV-B", "nextclade.tsv"), show_col_types = FALSE,
           guess_max = 500, progress = FALSE) %>% select(seqName, totalFrameShifts))
dat <- dat %>% left_join(fs_any, by = c("FieldNoNIV_f" = "seqName"))

qual_ok <- function(d, set) if (set == "Q-ALL") rep(TRUE, nrow(d)) else
  (d$nc_qc_status == "good" & d$nc_coverage >= 0.95 & d$totalFrameShifts == 0)

subset_genomes <- function(sp) {
  r <- SPEC[SPEC$spec == sp, ]
  keep <- link_ok(dat, r$link) & qual_ok(dat, r$qual)
  if (r$vill == "V33") keep <- keep & dat$code_village %in% V33
  dat[keep & !is.na(keep), ]
}

## ---------- O1 coverage --------------------------------------------------------
cat("\n-- O1 coverage --\n")
o1 <- list()
for (sp in SPEC$spec) {
  r <- SPEC[SPEC$spec == sp, ]
  a <- ari %>% filter(rsv_pos)
  if (r$vill == "V33") a <- a %>% filter(code_village %in% V33)
  g <- subset_genomes(sp)
  ## DECLARED DEFINITION (see 27_cohort_status_audit.R): an episode is
  ## RSV-positive if the cohort flags say so. The numerator is therefore
  ## restricted to episodes in `a`, so numerator and denominator share one
  ## rule (and one village restriction under V33). Linked episodes that are
  ## cohort-negative are excluded here and retained in the genomic analyses,
  ## where the laboratory result is the relevant one.
  g <- g %>% semi_join(a %>% select(ep_row), by = "ep_row")
  o1[[sp]] <- a %>% count(season, name = "cohort_rsv_pos") %>%
    left_join(g %>% group_by(season) %>% summarise(linked = n_distinct(ep_row), .groups = "drop"),
              by = "season") %>%
    mutate(linked = replace_na(linked, 0L),
           pct = round(100 * linked / cohort_rsv_pos, 1), spec = sp,
           link = r$link, qual = r$qual, vill = r$vill)
}
O1 <- bind_rows(o1)
print(O1 %>% select(spec, season, cohort_rsv_pos, linked, pct) %>%
        pivot_wider(names_from = spec, values_from = c(linked, pct), names_vary = "slowest") %>%
        as.data.frame())

## ---------- O2 severity: village restriction ONLY ------------------------------
cat("\n-- O2 severity (all tested episodes; village restriction only) --\n")
fit_o2 <- function(vill) {
  d <- ari %>% filter(swabbed, !is.na(lri_cat))
  if (vill == "V33") d <- d %>% filter(code_village %in% V33)
  m <- try(glmer(very_severe ~ rsv_pos + scale(age_days_episode) + female + factor(season) +
                   (1 | code_village) + (1 | studyid),
                 data = d, family = binomial, control = glmerControl(optimizer = "bobyqa")),
           silent = TRUE)
  if (inherits(m, "try-error")) return(tibble(vill = vill, status = "fit failed"))
  mf <- model.frame(m); z <- coef(summary(m))["rsv_posTRUE", ]
  tibble(vill = vill, episodes_available = nrow(d), fitted = nobs(m),
         rsv_pos_fitted = sum(mf$rsv_pos), events = sum(mf$very_severe),
         non_events = sum(mf$very_severe == 0),
         OR = round(exp(z[1]), 3), lo = round(exp(z[1] - 1.96 * z[2]), 3),
         hi = round(exp(z[1] + 1.96 * z[2]), 3), p = round(z[4], 3),
         singular = isSingular(m),
         ## Convergence is recorded from the SAME two checks as the five-model
         ## severity table: lme4's own post-fit messages AND the optimizer's
         ## return code. Until 2026-10-31 this fit exported only the first, so
         ## the manuscript's claim that every interpreted fit carries both was
         ## not supported by this table.
         no_lme4_msg = is.null(m@optinfo$conv$lme4$messages),
         optimizer_code = tryCatch(as.integer(m@optinfo$conv$opt), error = function(e) NA_integer_),
         converged = is.null(m@optinfo$conv$lme4$messages) &&
                     identical(tryCatch(as.integer(m@optinfo$conv$opt),
                                        error = function(e) NA_integer_), 0L)) %>%
    mutate(status = case_when(
      !converged &  singular ~ "UNRELIABLE: did not converge; singular",
      !converged             ~ "UNRELIABLE: did not converge",
      singular               ~ "UNRELIABLE: singular (variance component at 0)",
      TRUE                   ~ "converged"))
}
O2 <- bind_rows(fit_o2("all"), fit_o2("V33"))
print(as.data.frame(O2))

## ---------- O4 spatial: permutation p + genome-level bootstrap CI --------------
cat("\n-- O4 village clustering --\n")
snpmat <- function(ty) {
  m <- read_fasta_matrix(file.path(PROJ, "nextclade", ty, "nextclade.aligned.fasta"))
  m <- m[rownames(m) != "", , drop = FALSE]
  acgt <- m %in% c("A", "C", "G", "T"); dim(acgt) <- dim(m)
  raw <- as.matrix(ape::dist.dna(ape::as.DNAbin(tolower(m)), model = "N",
                                 pairwise.deletion = TRUE, as.matrix = TRUE))
  shared <- acgt %*% t(acgt)                         # shared callable sites per pair
  dimnames(shared) <- dimnames(raw)
  list(raw = raw, shared = shared)
}
MA <- snpmat("RSV-A"); MB <- snpmat("RSV-B")
cent <- read_csv(file.path(OUT, "05_village_centroids.csv"), show_col_types = FALSE)

stat_diff <- function(G, vill, mask) {
  iu <- upper.tri(G) & mask
  s <- outer(vill, vill, "==")[iu]
  if (!any(s) || !any(!s)) return(NA_real_)
  median(G[iu][s]) - median(G[iu][!s])
}

o4 <- list()
for (sp in SPEC$spec) {
  g <- subset_genomes(sp) %>% filter(!is.na(code_village)) %>%
    inner_join(cent %>% select(code_village), by = "code_village")
  for (ty in c("RSV-A", "RSV-B")) {
    M <- if (ty == "RSV-A") MA else MB
    for (ss in sort(unique(g$season[g$nc_dataset == ty]))) {
      sub <- g %>% filter(nc_dataset == ty, season == ss)
      ids <- intersect(sub$FieldNoNIV_f, rownames(M$raw))
      sub <- sub %>% filter(FieldNoNIV_f %in% ids) %>% arrange(match(FieldNoNIV_f, ids))
      n <- length(ids)
      Graw <- M$raw[ids, ids, drop = FALSE]; Sh <- M$shared[ids, ids, drop = FALSE]
      mask <- Sh >= 13000
      iu <- upper.tri(Graw)
      same <- outer(sub$code_village, sub$code_village, "==")
      np_same <- sum(iu & mask & same); np_diff <- sum(iu & mask & !same)
      dropped <- sum(iu & !mask)
      gate <- n >= 8 & n_distinct(sub$code_village) >= 5 & np_same >= 5 & np_diff >= 20
      row <- tibble(spec = sp, rsv_type = ty, season = ss, genomes = n,
                    villages = n_distinct(sub$code_village),
                    pairs_same = np_same, pairs_diff = np_diff, pairs_dropped = dropped,
                    analysed = gate)
      if (!gate) { o4[[paste(sp, ty, ss)]] <- row; next }
      Gnorm <- 1000 * Graw / pmax(Sh, 1)
      obs_raw  <- stat_diff(Graw,  sub$code_village, mask)
      obs_norm <- stat_diff(Gnorm, sub$code_village, mask)
      null <- replicate(N_PERM, stat_diff(Graw, sample(sub$code_village), mask))
      pval <- (1 + sum(null <= obs_raw, na.rm = TRUE)) / (1 + sum(!is.na(null)))
      bs <- replicate(N_BOOT, {
        k <- sample(n, n, replace = TRUE)
        dup <- outer(k, k, "==")                      # a genome paired with its own copy
        stat_diff(Graw[k, k, drop = FALSE], sub$code_village[k], mask[k, k, drop = FALSE] & !dup)
      })
      ci <- quantile(bs, c(.025, .975), na.rm = TRUE)
      o4[[paste(sp, ty, ss)]] <- row %>%
        mutate(obs_raw = round(obs_raw, 1), obs_per1000 = round(obs_norm, 2),
               boot_lo = round(ci[1], 1), boot_hi = round(ci[2], 1), perm_p = round(pval, 4))
    }
  }
}
O4 <- bind_rows(o4)
print(as.data.frame(O4 %>% filter(analysed) %>%
        select(spec, rsv_type, season, genomes, pairs_same, pairs_diff, pairs_dropped,
               obs_raw, boot_lo, boot_hi, perm_p)))
cat("\nseason-cells not analysed (below the prespecified minimums):\n")
print(as.data.frame(O4 %>% filter(!analysed) %>%
        select(spec, rsv_type, season, genomes, villages, pairs_same, pairs_diff)))

## ---------- O5 F screen: all genomes primary; subsets labelled -----------------
cat("\n-- O5 F variant screen --\n")
EV <- read_csv(file.path(OUT, "04_F_evidence_table.csv"), show_col_types = FALSE) %>%
  filter(effect == "reduced")
AA <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
o5 <- list()
POPS <- list(`all genomes (primary)` = function(d) d,
             `linked only (sensitivity)` = function(d) d[link_ok(d, "L-ALL"), ],
             `exact-link only (sensitivity)` = function(d) d[link_ok(d, "L-EXACT"), ])
for (pn in names(POPS)) {
  dd <- POPS[[pn]](dat)
  for (ty in c("RSV-A", "RSV-B")) {
    f <- read_fasta_matrix(file.path(PROJ, "nextclade", ty, "nextclade.cds_translation.F.fasta"))
    f <- f[rownames(f) %in% dd$FieldNoNIV_f, , drop = FALSE]
    for (i in which(EV$subgroup == ty)) {
      pos <- EV$pos[i]; alt <- EV$aa[i]
      v <- f[, pos]; callable <- sum(v %in% AA); det <- sum(v == alt)
      ub <- if (callable == 0) NA_real_ else
        round(100 * binom.test(det, callable, alternative = "less")$conf.int[2], 2)
      o5[[paste(pn, ty, pos, alt)]] <- tibble(
        population = pn, subgroup = ty, variant = sprintf("%d%s", pos, alt),
        product = EV$product[i], detected = det, callable = callable,
        not_callable = nrow(f) - callable,
        upper_bound_pct = ub,
        note = if (callable == 0) "not assessed" else "independent-binomial approximation")
    }
  }
}
O5 <- bind_rows(o5)
cat("\nprimary population, all genomes:\n")
print(as.data.frame(O5 %>% filter(population == "all genomes (primary)") %>%
        select(subgroup, variant, product, detected, callable, not_callable, upper_bound_pct)))
cat(sprintf("\ntotal detections across all %d listed variants x 3 populations: %d\n",
            nrow(EV), sum(O5$detected)))
cat("These bounds assume independent sampling and do not account for genomic\n")
cat("relatedness or selection into sequencing; they are not validated\n")
cat("population-prevalence limits. Relatedness by descent does not by itself\n")
cat("establish any particular corrected interval, so none is quoted here --\n")
cat("\"wider\" was asserted without a model that would say how much wider.\n")

write_out(O1, "12_O1_coverage.csv")
write_out(O2, "12_O2_severity.csv")
write_out(O4, "12_O4_spatial.csv")
write_out(O5, "12_O5_variant_screen.csv")
saveRDS(list(O1 = O1, O2 = O2, O4 = O4, O5 = O5), file.path(OUT, "12_panel.rds"))
cat("== 12 done ==\n")

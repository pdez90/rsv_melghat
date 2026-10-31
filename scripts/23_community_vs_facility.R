# 23_community_vs_facility.R -- what does community-based sequencing add over
# facility-based sequencing?
#
# VARIABLE VALIDATION FIRST (the reason this script starts with a check, not an
# analysis):
#
#   community, reached_facility and hosp_lri are PERFECTLY COLLINEAR in these
#   data -- community == 1 exactly when reached_facility == 0 and hosp_lri == 0.
#   They are one variable under three names. Only one is used below; the others
#   must never be entered into the same model.
#
#   referred_hosp_episode is NOT the same variable and does not behave as its
#   name suggests: 97 episodes are referred but did not reach a facility, while
#   only 1 referred episode reached one, and 57 reached without referral. That
#   is either a severe access-to-care signal or a coding problem. It is reported
#   and NOT used until resolved.
#
# The comparison is made with equal-sized subsampling, because the community
# group is larger (145 vs 58) and both lineage richness and mean pairwise
# distance are sensitive to sample size.
#
# NOTE ON THE VALIDATION ABOVE: reached_facility was CONSTRUCTED from community,
# so their agreement is guaranteed by construction and is a consistency check,
# not evidence about what the coding means clinically. Agreement with hosp_lri
# is informative but still does not establish the clinical definition.

source(file.path("scripts", "00_setup.R"))
set.seed(20261030)
cat("\n== 23 community vs facility ==\n")

B <- 999   # subsampling replicates

dat <- readRDS(file.path(OUT, "02_dat.rds"))

## ---- 1. variable validation, reported every run ----------------------------
cat("\n-- variable validation --\n")
v <- dat %>% filter(!is.na(community))
cat(sprintf("  community==1 <-> reached_facility==0 : %s\n",
            all(v$community == 1) == all(v$reached_facility == 0) &&
            all((v$community == 1) == (v$reached_facility == 0))))
cat(sprintf("  community==1 <-> hosp_lri==0         : %s\n",
            all((v$community == 1) == (v$hosp_lri == 0))))
cat("  -> collinear; using `community` only.\n")
print(with(dat, table(referred = referred_hosp_episode, reached = reached_facility,
                      useNA = "ifany")))
cat("  -> referred_hosp_episode NOT used (see header).\n")

d <- dat %>%
  filter(link_status %in% c("exact", "nearest_prior_episode_same_RSV_type"),
         !is.na(community), !is.na(clade)) %>%
  mutate(setting = if_else(community == 1, "community-managed", "facility-attending"))

cat(sprintf("\n  %d linked sequenced episodes with a setting: %d community, %d facility\n",
            nrow(d), sum(d$setting == "community-managed"),
            sum(d$setting == "facility-attending")))

## ---- 2. lineage representation, by season ----------------------------------
cat("\n-- lineages captured, by season and setting --\n")
tab <- d %>% count(nc_dataset, season, setting, clade) %>%
  group_by(nc_dataset, season, setting) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup()
print(as.data.frame(tab), right = FALSE)

## rarefied lineage richness: subsample both settings to the smaller n
cat("\n-- lineage richness, rarefied to the smaller group (median over", B, "draws) --\n")
rich <- d %>% group_by(nc_dataset, season) %>% group_modify(function(x, key) {
  nc <- sum(x$setting == "community-managed"); nf <- sum(x$setting == "facility-attending")
  k <- min(nc, nf)
  if (k < 5) return(tibble(n_comm = nc, n_fac = nf, k = k,
                           rich_comm = NA_real_, rich_fac = NA_real_,
                           note = "insufficient: smaller group < 5"))
  rc <- replicate(B, n_distinct(sample(x$clade[x$setting == "community-managed"], k)))
  rf <- replicate(B, n_distinct(sample(x$clade[x$setting == "facility-attending"], k)))
  tibble(n_comm = nc, n_fac = nf, k = k,
         rich_comm = median(rc), rich_fac = median(rf),
         note = sprintf("community %.0f vs facility %.0f lineages at n=%d", median(rc), median(rf), k))
}) %>% ungroup()
print(as.data.frame(rich), right = FALSE)

## ---- 3. genetic diversity, equal-sized subsampling -------------------------
snpmat <- function(ty) {
  m <- read_fasta_matrix(file.path(PROJ, "nextclade", ty, "nextclade.aligned.fasta"))
  m <- m[rownames(m) != "", , drop = FALSE]
  acgt <- m %in% c("A", "C", "G", "T"); dim(acgt) <- dim(m)
  raw <- as.matrix(ape::dist.dna(ape::as.DNAbin(tolower(m)), model = "N",
                                 pairwise.deletion = TRUE, as.matrix = TRUE))
  sh <- acgt %*% t(acgt); dimnames(sh) <- dimnames(raw)
  list(raw = raw, shared = sh)
}
MA <- snpmat("RSV-A"); MB <- snpmat("RSV-B")

cat("\n-- mean pairwise distance within setting, rarefied to equal n --\n")
cat("   ordering_consistency is DESCRIPTIVE; inference is p_label_perm\n")
div <- d %>% group_by(nc_dataset, season) %>% group_modify(function(x, key) {
  M <- if (key$nc_dataset == "RSV-A") MA else MB
  ids <- intersect(x$FieldNoNIV_f, rownames(M$raw))
  x <- x %>% filter(FieldNoNIV_f %in% ids)
  gc <- x$FieldNoNIV_f[x$setting == "community-managed"]
  gf <- x$FieldNoNIV_f[x$setting == "facility-attending"]
  k <- min(length(gc), length(gf))
  if (k < 5) return(tibble(n_comm = length(gc), n_fac = length(gf), k = k,
                           d_comm = NA_real_, d_fac = NA_real_, p_diff = NA_real_,
                           note = "insufficient: smaller group < 5"))
  mp <- function(s) {
    sub <- M$raw[s, s, drop = FALSE]; sh <- M$shared[s, s, drop = FALSE]
    iu <- upper.tri(sub) & (sh >= 13000)
    if (!any(iu)) return(NA_real_)
    mean(1000 * sub[iu] / pmax(sh[iu], 1))
  }
  ## Statistic: DETERMINISTIC difference in mean pairwise distance, computed on
  ## all usable pairs within each setting. Unlike lineage richness, the
  ## expectation of mean pairwise distance does not depend on sample size, so no
  ## rarefaction is needed for the test and none is used -- rarefaction is kept
  ## below as a DESCRIPTIVE summary only.
  stat <- function(labels) {
    a <- x$FieldNoNIV_f[labels == "community-managed"]
    b <- x$FieldNoNIV_f[labels != "community-managed"]
    if (length(a) < 5 || length(b) < 5) return(NA_real_)
    mp(a) - mp(b)
  }
  obs <- stat(x$setting)
  ## NULL: permute the SETTING LABEL within this season stratum.
  ## The earlier version used mean(dc >= df) over subsamples drawn separately
  ## within the observed groups. That is subsampling consistency conditional on
  ## these samples, not a null of "no setting difference", and is not a p-value.
  ##
  ## Exchangeability caveat: this tests setting labels as exchangeable WITHIN a
  ## season. It does not adjust for collection time, village or lineage, any of
  ## which could differ by setting and drive the statistic.
  NPERM <- 9999
  nulld <- replicate(NPERM, stat(sample(x$setting)))
  okp <- sum(!is.na(nulld))
  pperm <- if (okp < 200) NA_real_ else
    (1 + sum(abs(nulld) >= abs(obs), na.rm = TRUE)) / (1 + okp)
  mc_se <- if (is.na(pperm)) NA_real_ else sqrt(pperm * (1 - pperm) / okp)
  ## descriptive rarefaction only
  dc <- replicate(B, mp(sample(gc, k))); df <- replicate(B, mp(sample(gf, k)))
  consistency <- mean(dc >= df, na.rm = TRUE)
  tibble(n_comm = length(gc), n_fac = length(gf), k = k,
         ## TWO DIFFERENT QUANTITIES, kept in separately named columns. Until
         ## 2026-10-31 the rarefied medians were exported under d_comm/d_fac
         ## beside obs_diff, which is the difference of ALL-PAIR means, so the
         ## columns did not subtract to the difference and a table built from
         ## them mixed a descriptive with an inferential statistic.
         d_comm_allpair = round(mp(x$FieldNoNIV_f[x$setting == "community-managed"]), 3),
         d_fac_allpair  = round(mp(x$FieldNoNIV_f[x$setting != "community-managed"]), 3),
         d_comm_rarefied_median = round(median(dc, na.rm = TRUE), 3),
         d_fac_rarefied_median  = round(median(df, na.rm = TRUE), 3),
         obs_diff = round(obs, 3),
         ordering_consistency = round(consistency, 4),
         p_label_perm = round(pperm, 4), p_mc_se = round(mc_se, 4), perms_usable = okp,
         note = "obs_diff is deterministic (all pairs, no rarefaction); ordering_consistency is DESCRIPTIVE, not a p-value; p_label_perm is a two-sided setting-label permutation within season, not adjusted for time, village or lineage")
}) %>% ungroup()
print(as.data.frame(div), right = FALSE)

## ---- 4. what would facility-only sequencing have missed? -------------------
cat("\n-- lineages seen ONLY in community-managed episodes --\n")
only <- d %>% group_by(nc_dataset, season, clade) %>%
  summarise(n_comm = sum(setting == "community-managed"),
            n_fac = sum(setting == "facility-attending"), .groups = "drop") %>%
  filter(n_fac == 0, n_comm > 0)
if (nrow(only)) print(as.data.frame(only), right = FALSE) else cat("  none\n")
cat("  NOTE: with 145 community vs 58 facility genomes, a lineage absent from the\n")
cat("  facility group may reflect sample size alone. Compare against the rarefied\n")
cat("  richness above rather than reading this table on its own.\n")

write_csv(tab,  file.path(OUT, "23_lineage_by_setting.csv"))
write_csv(rich, file.path(OUT, "23_rarefied_richness.csv"))
write_csv(div,  file.path(OUT, "23_diversity_by_setting.csv"))
cat("\n  -> wrote 23_lineage_by_setting.csv, 23_rarefied_richness.csv, 23_diversity_by_setting.csv\n")
cat("== 23 done ==\n")

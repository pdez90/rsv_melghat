# 14_village_influence.R -- leave-one-village-out influence analysis for the
# village-clustering result, under spec S0.
#
# Purpose: the within-village comparisons in each stratum are concentrated in a
# few villages. This quantifies that concentration directly: for every village
# that contributes same-village pairs, drop all its genomes and recompute the
# statistic and its permutation p-value on the remaining data.
#
# Reports BOTH raw SNP differences and normalised (per 1,000 callable sites)
# differences, since script 12 carried inference on raw distances only.

source(file.path("scripts", "00_setup.R"))
set.seed(20261030)
cat("\n== 14 village influence (leave-one-village-out, S0) ==\n")

N_PERM <- 1999
THRESH <- 13000

snpmat <- function(ty) {
  m <- read_fasta_matrix(file.path(PROJ, "nextclade", ty, "nextclade.aligned.fasta"))
  m <- m[rownames(m) != "", , drop = FALSE]
  acgt <- m %in% c("A", "C", "G", "T"); dim(acgt) <- dim(m)
  raw <- as.matrix(ape::dist.dna(ape::as.DNAbin(tolower(m)), model = "N",
                                 pairwise.deletion = TRUE, as.matrix = TRUE))
  shared <- acgt %*% t(acgt); dimnames(shared) <- dimnames(raw)
  list(raw = raw, shared = shared)
}
MA <- snpmat("RSV-A"); MB <- snpmat("RSV-B")

dat  <- readRDS(file.path(OUT, "02_dat.rds"))
cent <- read_csv(file.path(OUT, "05_village_centroids.csv"), show_col_types = FALSE)

stat_diff <- function(G, vill, mask) {
  iu <- upper.tri(G) & mask
  s <- outer(vill, vill, "==")[iu]
  if (!any(s) || !any(!s)) return(NA_real_)
  median(G[iu][s]) - median(G[iu][!s])
}

# S0 = all links, all QC, all villages
g0 <- dat %>%
  filter(link_status %in% c("exact", "nearest_prior_episode_same_RSV_type"),
         !is.na(code_village)) %>%
  inner_join(cent %>% select(code_village), by = "code_village")

STRATA <- tribble(~ty, ~ss,
                  "RSV-A", "2016/17", "RSV-A", "2019/20",
                  "RSV-A", "2021/22", "RSV-B", "2018/19")

# one evaluation of a stratum on a given genome subset
evaluate <- function(sub, M) {
  ids  <- intersect(sub$FieldNoNIV_f, rownames(M$raw))
  sub  <- sub %>% filter(FieldNoNIV_f %in% ids) %>% arrange(match(FieldNoNIV_f, ids))
  n    <- length(ids)
  if (n < 8) return(NULL)
  Graw <- M$raw[ids, ids, drop = FALSE]; Sh <- M$shared[ids, ids, drop = FALSE]
  mask <- Sh >= THRESH
  iu   <- upper.tri(Graw); same <- outer(sub$code_village, sub$code_village, "==")
  np_same <- sum(iu & mask & same); np_diff <- sum(iu & mask & !same)
  if (n_distinct(sub$code_village) < 5 || np_same < 5 || np_diff < 20) return(NULL)
  Gnorm <- 1000 * Graw / pmax(Sh, 1)
  o_raw  <- stat_diff(Graw,  sub$code_village, mask)
  o_norm <- stat_diff(Gnorm, sub$code_village, mask)
  nr <- replicate(N_PERM, stat_diff(Graw,  sample(sub$code_village), mask))
  nn <- replicate(N_PERM, stat_diff(Gnorm, sample(sub$code_village), mask))
  tibble(genomes = n, villages = n_distinct(sub$code_village),
         pairs_same = np_same, pairs_diff = np_diff,
         obs_raw  = round(o_raw, 1),
         p_raw    = round((1 + sum(nr <= o_raw,  na.rm = TRUE)) / (1 + sum(!is.na(nr))), 4),
         obs_norm = round(o_norm, 3),
         p_norm   = round((1 + sum(nn <= o_norm, na.rm = TRUE)) / (1 + sum(!is.na(nn))), 4))
}

res <- list()
for (i in seq_len(nrow(STRATA))) {
  ty <- STRATA$ty[i]; ss <- STRATA$ss[i]
  M  <- if (ty == "RSV-A") MA else MB
  sub_full <- g0 %>% filter(nc_dataset == ty, season == ss)
  full <- evaluate(sub_full, M)
  if (is.null(full)) next
  res[[paste(ty, ss, "FULL")]] <- full %>%
    mutate(rsv_type = ty, season = ss, dropped_village = "(none: full stratum)",
           genomes_dropped = 0L, same_pairs_contributed = NA_integer_, .before = 1)

  # villages that actually contribute within-village pairs
  tb <- table(sub_full$code_village)
  contrib <- names(tb)[tb >= 2]
  for (v in contrib) {
    out <- evaluate(sub_full %>% filter(code_village != v), M)
    row <- if (is.null(out))
      tibble(genomes = NA_integer_, villages = NA_integer_, pairs_same = NA_integer_,
             pairs_diff = NA_integer_, obs_raw = NA_real_, p_raw = NA_real_,
             obs_norm = NA_real_, p_norm = NA_real_)
    else out
    res[[paste(ty, ss, v)]] <- row %>%
      mutate(rsv_type = ty, season = ss, dropped_village = v,
             genomes_dropped = as.integer(tb[[v]]),
             same_pairs_contributed = as.integer(choose(tb[[v]], 2)), .before = 1)
  }
}

INF <- bind_rows(res) %>%
  group_by(rsv_type, season) %>%
  mutate(full_obs = obs_raw[dropped_village == "(none: full stratum)"][1],
         shift    = round(obs_raw - full_obs, 1),
         verdict  = case_when(
           dropped_village == "(none: full stratum)" ~ "-",
           is.na(obs_raw)                            ~ "stratum no longer estimable",
           p_raw >= 0.05                             ~ "signal lost (p >= 0.05)",
           TRUE                                      ~ "signal retained")) %>%
  ungroup() %>%
  select(rsv_type, season, dropped_village, genomes_dropped, same_pairs_contributed,
         genomes, pairs_same, pairs_diff, obs_raw, shift, p_raw, obs_norm, p_norm, verdict)

print(as.data.frame(INF), right = FALSE)
write_csv(INF, file.path(OUT, "14_village_influence.csv"))
cat("\n  -> wrote 14_village_influence.csv\n== 14 done ==\n")

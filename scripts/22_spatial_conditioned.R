# 22_spatial_conditioned.R -- does village-associated genetic similarity survive
# conditioning on WHEN and WHICH virus was sampled?
#
# Leave-one-village-out (script 14) showed no single village drives the signal.
# It does not address confounding: genomes collected days apart, or belonging to
# the same lineage, are more similar for reasons that have nothing to do with
# village. If sampling within a village is also clustered in time or lineage --
# which it will be, because a village outbreak is a short burst of one lineage --
# the unconditioned permutation test cannot separate the two.
#
# Three approaches, all on the same statistic (median same-village SNP distance
# minus median different-village SNP distance, negative = villages more similar):
#
#   U   unconditioned        permute village labels freely        (baseline)
#   C   conditioned          permute village labels ONLY WITHIN
#                            lineage x time-bin blocks
#   R   restricted pairs     compute the statistic using only pairs matched on
#                            lineage, and on lineage + collection window
#
# Also closes two deferred sensitivity items: the pair-inclusion threshold is run
# at 12,000 / 13,000 / 14,000 shared callable sites, and the bootstrap is run on
# the NORMALISED statistic (substitutions per 1,000 callable sites) as well as
# the raw one.
#
# WHAT THIS IS NOT. Changing the permutation scheme conditions the NULL. It does
# not convert the statistic, or its bootstrap interval, into a lineage-adjusted
# effect estimate. The unconditioned and same-lineage rows are medians over
# DIFFERENT SETS OF PAIRS, so the drop from one to the other is a descriptive
# contrast, not a decomposition of the effect into lineage and village parts.
#
# Permutation freedom is measured empirically below, by counting how often each
# pair actually changes same/different-village status across draws. Counting
# pairs that happen to lie inside one block does NOT measure this: a pair split
# across blocks can still be broken, because each genome independently draws a
# label from its own block.
#
# Power is limited by design here. Conditioning removes most of the contrasts,
# so "not estimable" is a legitimate and expected result and is reported as such
# rather than as a null finding.

source(file.path("scripts", "00_setup.R"))
set.seed(20261030)
cat("\n== 22 spatial, conditioned on time and lineage ==\n")

N_PERM  <- 1999
N_BOOT  <- 999
THRESH  <- c(12000, 13000, 14000)
TIMEBIN <- 14      # days per block for the conditioned permutation
WINDOWS <- c(14, 28)

dat <- readRDS(file.path(OUT, "02_dat.rds"))
cent <- read_csv(file.path(OUT, "05_village_centroids.csv"), show_col_types = FALSE)

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

## statistic over a boolean pair-inclusion matrix
stat_diff <- function(G, vill, use) {
  iu <- upper.tri(G) & use
  if (!any(iu)) return(NA_real_)
  s <- outer(vill, vill, "==")[iu]
  if (sum(s) < 5 || sum(!s) < 20) return(NA_real_)
  median(G[iu][s]) - median(G[iu][!s])
}

## permute labels within blocks; singleton blocks keep their label
perm_within <- function(vill, block) {
  out <- vill
  for (b in unique(block)) {
    i <- which(block == b)
    if (length(i) > 1) out[i] <- sample(vill[i])
  }
  out
}

STRATA <- tribble(~ty, ~ss, "RSV-B", "2018/19", "RSV-A", "2021/22")
res <- list()

for (k in seq_len(nrow(STRATA))) {
  ty <- STRATA$ty[k]; ss <- STRATA$ss[k]
  M  <- if (ty == "RSV-A") MA else MB
  g  <- dat %>%
    filter(link_status %in% c("exact", "nearest_prior_episode_same_RSV_type"),
           !is.na(code_village), nc_dataset == ty, season == ss) %>%
    inner_join(cent %>% select(code_village), by = "code_village")
  ids <- intersect(g$FieldNoNIV_f, rownames(M$raw))
  g   <- g %>% filter(FieldNoNIV_f %in% ids) %>% arrange(match(FieldNoNIV_f, ids))
  n   <- length(ids)
  Graw <- M$raw[ids, ids, drop = FALSE]; Sh <- M$shared[ids, ids, drop = FALSE]
  Gnorm <- 1000 * Graw / pmax(Sh, 1)
  vill <- g$code_village; lin <- g$clade; dte <- as.numeric(g$swab_date)
  blk  <- paste(lin, floor((dte - min(dte)) / TIMEBIN))
  same_lin <- outer(lin, lin, "==")
  dgap <- abs(outer(dte, dte, "-"))

  cat(sprintf("\n[%s %s] %d genomes, %d villages, %d lineages, %d-day span\n",
              ty, ss, n, n_distinct(vill), n_distinct(lin),
              as.integer(diff(range(dte)))))
  ## How much freedom does the conditioned permutation actually have?
  ##
  ## NOT by counting pairs that sit inside a block. A same-village pair split
  ## across two blocks is NOT fixed: each genome draws a label from its own
  ## block's pool, so either can acquire a different village and break the pair.
  ## The only honest diagnostic is empirical -- run the permutations and record,
  ## for every pair, how often its same/different-village status actually
  ## changes.
  nv_blk <- tapply(vill, blk, n_distinct)
  iu0 <- upper.tri(matrix(0, n, n))
  obs_same <- outer(vill, vill, "==")[iu0]
  switch_count <- integer(sum(iu0))
  NDIAG <- 499
  for (r in seq_len(NDIAG)) {
    pv <- perm_within(vill, blk)
    switch_count <- switch_count + (outer(pv, pv, "==")[iu0] != obs_same)
  }
  sw <- switch_count / NDIAG
  cat(sprintf("  conditioning blocks (lineage x %dd): %d total, %d with >1 genome, %d with >1 village\n",
              TIMEBIN, n_distinct(blk), sum(table(blk) > 1), sum(nv_blk > 1)))
  cat(sprintf("  permutation freedom, measured over %d draws:\n", NDIAG))
  cat(sprintf("    same-village pairs   : %d, of which %d (%.0f%%) ever change status; median switch rate %.3f\n",
              sum(obs_same), sum(sw[obs_same] > 0),
              100 * mean(sw[obs_same] > 0), median(sw[obs_same])))
  cat(sprintf("    different-village    : %d, of which %d (%.0f%%) ever change status; median switch rate %.3f\n",
              sum(!obs_same), sum(sw[!obs_same] > 0),
              100 * mean(sw[!obs_same] > 0), median(sw[!obs_same])))

  for (th in THRESH) {
    mask <- Sh >= th
    iu <- upper.tri(Graw)
    np_s <- sum(iu & mask & outer(vill, vill, "==")); np_d <- sum(iu & mask & !outer(vill, vill, "=="))

    row <- tibble(rsv_type = ty, season = ss, threshold = th, genomes = n,
                  pairs_same = np_s, pairs_diff = np_d)

    if (np_s < 5 || np_d < 20) {
      cat(sprintf("  thr %5d: same %3d / diff %5d -> NOT ESTIMABLE (needs >=5 same, >=20 diff)\n",
                  th, np_s, np_d))
      res[[paste(ty, ss, th)]] <- row %>% mutate(analysed = FALSE)
      next
    }

    o_raw  <- stat_diff(Graw,  vill, mask)
    o_norm <- stat_diff(Gnorm, vill, mask)

    ## U: unconditioned permutation
    nullU <- replicate(N_PERM, stat_diff(Graw, sample(vill), mask))
    pU <- (1 + sum(nullU <= o_raw, na.rm = TRUE)) / (1 + sum(!is.na(nullU)))

    ## C: conditioned permutation (within lineage x time block)
    nullC <- replicate(N_PERM, stat_diff(Graw, perm_within(vill, blk), mask))
    okC <- sum(!is.na(nullC))
    pC <- if (okC < 100) NA_real_ else (1 + sum(nullC <= o_raw, na.rm = TRUE)) / (1 + okC)

    ## C on the normalised scale
    nullCn <- replicate(N_PERM, stat_diff(Gnorm, perm_within(vill, blk), mask))
    okCn <- sum(!is.na(nullCn))
    pCn <- if (okCn < 100) NA_real_ else (1 + sum(nullCn <= o_norm, na.rm = TRUE)) / (1 + okCn)

    ## bootstrap CI, raw AND normalised, genome-level, self-pairs excluded
    bs <- function(G) replicate(N_BOOT, {
      i <- sample(n, n, replace = TRUE)
      dup <- outer(i, i, "==")
      stat_diff(G[i, i, drop = FALSE], vill[i], mask[i, i, drop = FALSE] & !dup)
    })
    bR <- bs(Graw); bN <- bs(Gnorm)
    ciR <- quantile(bR, c(.025, .975), na.rm = TRUE)
    ciN <- quantile(bN, c(.025, .975), na.rm = TRUE)

    cat(sprintf("  thr %5d: same %3d / diff %5d | raw %7.1f [%.0f, %.0f] pU=%.4f pC=%s\n",
                th, np_s, np_d, o_raw, ciR[1], ciR[2], pU,
                ifelse(is.na(pC), "n/e", sprintf("%.4f", pC))))
    cat(sprintf("             norm %6.2f per 1000 [%.2f, %.2f] pC=%s  (%d/%d conditioned perms usable)\n",
                o_norm, ciN[1], ciN[2],
                ifelse(is.na(pCn), "n/e", sprintf("%.4f", pCn)), okC, N_PERM))

    res[[paste(ty, ss, th)]] <- row %>%
      mutate(analysed = TRUE,
             obs_raw = round(o_raw, 1), raw_lo = round(ciR[1], 1), raw_hi = round(ciR[2], 1),
             obs_norm = round(o_norm, 3), norm_lo = round(ciN[1], 3), norm_hi = round(ciN[2], 3),
             p_uncond = round(pU, 4), p_cond_raw = round(pC, 4), p_cond_norm = round(pCn, 4),
             cond_perms_usable = okC, n_blocks = n_distinct(blk),
             blocks_gt1 = sum(table(blk) > 1))
  }

  ## R: pair-restriction, at the 13,000 threshold only
  mask <- Sh >= 13000
  for (lab in c("same lineage", paste0("same lineage + <=", WINDOWS, "d"))) {
    use <- mask & same_lin
    if (grepl("<=", lab)) use <- use & (dgap <= as.numeric(sub(".*<=(\\d+)d", "\\1", lab)))
    iu <- upper.tri(Graw)
    ns <- sum(iu & use & outer(vill, vill, "==")); nd <- sum(iu & use & !outer(vill, vill, "=="))
    if (ns < 5 || nd < 20) {
      cat(sprintf("  restricted to %-26s: same %3d / diff %5d -> NOT ESTIMABLE\n", lab, ns, nd))
      res[[paste(ty, ss, lab)]] <- tibble(rsv_type = ty, season = ss, threshold = 13000,
        restriction = lab, genomes = n, pairs_same = ns, pairs_diff = nd, analysed = FALSE)
      next
    }
    o <- stat_diff(Graw, vill, use)
    nl <- replicate(N_PERM, stat_diff(Graw, perm_within(vill, blk), use))
    ok <- sum(!is.na(nl))
    p <- if (ok < 100) NA_real_ else (1 + sum(nl <= o, na.rm = TRUE)) / (1 + ok)
    ## genome-level bootstrap ON THE RESTRICTED PAIR SET, so the interval
    ## describes the statistic actually being reported rather than the
    ## unrestricted one
    bsr <- replicate(N_BOOT, {
      i <- sample(n, n, replace = TRUE)
      dup <- outer(i, i, "==")
      stat_diff(Graw[i, i, drop = FALSE], vill[i], use[i, i, drop = FALSE] & !dup)
    })
    nb <- sum(!is.na(bsr))
    ci <- if (nb < 100) c(NA_real_, NA_real_) else quantile(bsr, c(.025, .975), na.rm = TRUE)
    cat(sprintf("  restricted to %-26s: same %3d / diff %5d | raw %7.1f [%s, %s] pC=%s (%d/%d boots usable)\n",
                lab, ns, nd, o,
                ifelse(is.na(ci[1]), "n/e", sprintf("%.1f", ci[1])),
                ifelse(is.na(ci[2]), "n/e", sprintf("%.1f", ci[2])),
                ifelse(is.na(p), "n/e", sprintf("%.4f", p)), nb, N_BOOT))
    res[[paste(ty, ss, lab)]] <- tibble(rsv_type = ty, season = ss, threshold = 13000,
      restriction = lab, genomes = n, pairs_same = ns, pairs_diff = nd, analysed = TRUE,
      obs_raw = round(o, 1), raw_lo = round(ci[1], 1), raw_hi = round(ci[2], 1),
      p_cond_raw = round(p, 4), cond_perms_usable = ok, boots_usable = nb)
  }
}

OUTt <- bind_rows(res)
write_csv(OUTt, file.path(OUT, "22_spatial_conditioned.csv"))
cat("\n  -> wrote 22_spatial_conditioned.csv\n== 22 done ==\n")

# 19_mask_comparison.R -- does excluding uncertain regions change the Melghat
# conclusions? Compares the full-alignment trees with the masked ones
# (low-occupancy columns + the whole 72-nt G duplication block removed).
#
# The comparison is deliberately about the quantities that would be reported --
# Melghat clade membership and the support on Melghat-relevant branches -- not
# about whole-tree topological distance, which is dominated by 2,000 contextual
# tips nobody is making claims about.

source(file.path("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(ape)})
cat("\n== 19 mask comparison ==\n")
dat <- readRDS(file.path(OUT, "02_dat.rds"))

descendants <- function(phy) {
  n_tip <- Ntip(phy); n_node <- Nnode(phy)
  kids <- vector("list", n_tip + n_node)
  for (i in seq_len(nrow(phy$edge))) kids[[phy$edge[i,1]]] <- c(kids[[phy$edge[i,1]]], phy$edge[i,2])
  out <- vector("list", n_tip + n_node)
  for (i in seq_len(n_tip)) out[[i]] <- i
  repeat { done <- TRUE
    for (nd in (n_tip+1):(n_tip+n_node)) {
      if (!is.null(out[[nd]])) next
      ch <- kids[[nd]]
      if (all(!sapply(out[ch], is.null))) out[[nd]] <- unlist(out[ch]) else done <- FALSE }
    if (done) break }
  out
}

profile <- function(phy, mel_names) {
  n_tip <- Ntip(phy); is_mel <- phy$tip.label %in% mel_names
  sup <- suppressWarnings(as.numeric(phy$node.label))
  des <- descendants(phy); nid <- (n_tip+1):(n_tip+Nnode(phy))
  parent <- integer(n_tip + Nnode(phy)); parent[phy$edge[,2]] <- phy$edge[,1]
  pure <- sapply(nid, function(nd) all(is_mel[des[[nd]]]))
  maximal <- nid[pure][sapply(nid[pure], function(nd) {
    p <- parent[nd]; p == 0 || !pure[match(p, nid)] || is.na(match(p, nid)) })]
  cl <- lapply(maximal, function(nd) sort(phy$tip.label[des[[nd]]]))
  inside <- unique(unlist(cl))
  place <- sapply(maximal, function(nd) sup[match(parent[nd], nid)])
  within <- sup[match(maximal, nid)]
  list(clades = cl, clades_all = list(phy$tip.label[is_mel]),
       n_clades = length(cl), n_inside = length(inside),
       n_single = sum(is_mel) - length(inside),
       place = place[!is.na(place)], within = within[!is.na(within)],
       all = sup[!is.na(sup)])
}

rows <- list()
for (ty in c("RSV-A", "RSV-B")) {
  s <- substr(ty, 5, 5)
  f1 <- file.path(PROJ, "contextual", sprintf("tree_RSV-%s.nwk", s))
  f2 <- file.path(PROJ, "contextual", sprintf("tree_masked_RSV-%s.nwk", s))
  if (!file.exists(f1) || !file.exists(f2)) { cat(sprintf("  [%s] missing tree\n", ty)); next }
  mel <- dat$FieldNoNIV_f[dat$nc_dataset == ty]
  A <- profile(read.tree(f1), mel); B <- profile(read.tree(f2), mel)

  ## how many maximal Melghat clades are identical sets in both trees?
  keyA <- sapply(A$clades, paste, collapse = "|")
  keyB <- sapply(B$clades, paste, collapse = "|")
  shared <- length(intersect(keyA, keyB))

  ## Pair-level comparison over EVERY Melghat genome present in both trees.
  ## A genome that sits singly among contextual tips is its OWN group -- the
  ## previous version dropped such genomes before computing agreement, which
  ## excluded exactly the cases where membership changed and guaranteed a high
  ## number. Singletons get distinct negative group ids so they can never be
  ## counted as co-clustered with anything.
  grp <- function(P, genomes) {
    g <- setNames(rep(NA_integer_, length(genomes)), genomes)
    for (i in seq_along(P$clades)) g[intersect(P$clades[[i]], genomes)] <- i
    singles <- which(is.na(g))
    if (length(singles)) g[singles] <- -seq_along(singles)
    g
  }
  gg <- intersect(mel, intersect(unlist(A$clades_all), unlist(B$clades_all)))
  gA <- grp(A, gg); gB <- grp(B, gg)
  pa <- outer(gA, gA, "=="); pb <- outer(gB, gB, "==")
  iu <- upper.tri(pa)
  a <- pa[iu]; b <- pb[iu]
  retained   <- sum(a & b)        # together in both
  separated  <- sum(a & !b)       # together in full, apart after masking
  newly      <- sum(!a & b)       # apart in full, together after masking
  apart_both <- sum(!a & !b)
  agree      <- (retained + apart_both) / length(a)
  ## membership changes at genome level
  moved <- sum(mapply(function(x, y) {
    # did this genome's co-set change at all?
    !identical(sort(names(gA)[gA == gA[x]]), sort(names(gB)[gB == gB[x]]))
  }, seq_along(gg), seq_along(gg)))

  cat(sprintf("\n[%s]\n", ty))
  cat(sprintf("  maximal Melghat clades : full %d -> masked %d (%d identical sets)\n",
              A$n_clades, B$n_clades, shared))
  cat(sprintf("  genomes inside a clade : full %d -> masked %d (singletons %d -> %d)\n",
              A$n_inside, B$n_inside, A$n_single, B$n_single))
  cat(sprintf("  Melghat genomes compared: %d (singletons counted as their own group)\n", length(gg)))
  cat(sprintf("  pairs together in BOTH (retained)      : %d\n", retained))
  cat(sprintf("  pairs SEPARATED by masking             : %d\n", separated))
  cat(sprintf("  pairs NEWLY JOINED by masking          : %d\n", newly))
  cat(sprintf("  genomes whose co-set changed at all    : %d of %d\n", moved, length(gg)))
  cat(sprintf("  overall pair agreement                 : %.2f%%\n", 100 * agree))
  f <- function(x) sprintf("median %.3f, >=0.95 %.0f%%", median(x), 100*mean(x >= 0.95))
  cat(sprintf("  support, all nodes     : full %s | masked %s\n", f(A$all), f(B$all)))
  cat(sprintf("  support, placement     : full %s | masked %s\n", f(A$place), f(B$place)))
  # NOTE: this is the CLADE-DEFINING node of each maximal Melghat clade -- i.e.
  # support for the grouping. It is NOT support for branching order inside the
  # clade; script 16 reports that separately and it is much lower.
  cat(sprintf("  support, clade-DEFINING: full %s | masked %s\n", f(A$within), f(B$within)))

  rows[[ty]] <- tibble(rsv_type = ty,
    clades_full = A$n_clades, clades_masked = B$n_clades, clades_identical = shared,
    inside_full = A$n_inside, inside_masked = B$n_inside,
    singletons_full = A$n_single, singletons_masked = B$n_single,
    n_compared = length(gg), pairs_retained = retained,
    pairs_separated = separated, pairs_newly_joined = newly,
    genomes_coset_changed = moved, pair_agreement = round(agree, 4),
    all_median_full = round(median(A$all),3),   all_median_masked = round(median(B$all),3),
    place_median_full = round(median(A$place),3), place_median_masked = round(median(B$place),3),
    clade_defining_median_full = round(median(A$within),3),
    clade_defining_median_masked = round(median(B$within),3))
}
if (length(rows)) {
  write_csv(bind_rows(rows), file.path(OUT, "19_mask_comparison.csv"))
  cat("\n  -> wrote 19_mask_comparison.csv\n")
}
cat("== 19 done ==\n")

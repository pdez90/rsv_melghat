# 16_tree_annotation.R -- read the exploratory FastTree trees and ask the
# questions that actually bear on Melghat placement.
#
# Why this exists: an overall support distribution says nothing about whether
# the branches relevant to THESE genomes are resolved. A tree can be 60%
# unresolved across 2,000 contextual tips and still place every Melghat genome
# on a well-supported branch, or the reverse. So this script separates:
#
#   (a) nodes whose descendants are entirely Melghat  -- within-Melghat structure
#   (b) the parent of each maximal Melghat clade      -- its PLACEMENT in context
#   (c) every node in the tree                        -- the baseline to compare to
#
# Node values are SH-like LOCAL support (Price, Dehal & Arkin 2010,
# PLoS ONE 5:e9490), not bootstrap. "Bad splits" in the FastTree log is an
# optimisation diagnostic -- the number of splits where NNI could not be
# resolved to a local optimum -- NOT a count of biologically wrong branches.

source(file.path("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(ape)})
cat("\n== 16 tree annotation ==\n")

dat <- readRDS(file.path(OUT, "02_dat.rds"))

## descendant tip sets, post-order, without extra dependencies
descendants <- function(phy) {
  n_tip <- Ntip(phy); n_node <- Nnode(phy)
  kids <- vector("list", n_tip + n_node)
  for (i in seq_len(nrow(phy$edge))) {
    p <- phy$edge[i, 1]; c <- phy$edge[i, 2]
    kids[[p]] <- c(kids[[p]], c)
  }
  out <- vector("list", n_tip + n_node)
  ord <- rev(postorder <- reorder(phy, "postorder")$edge[, 2])
  for (i in seq_len(n_tip)) out[[i]] <- i
  # process internal nodes in reverse topological order
  stack <- (n_tip + 1):(n_tip + n_node)
  repeat {
    done <- TRUE
    for (nd in stack) {
      if (!is.null(out[[nd]])) next
      ch <- kids[[nd]]
      if (all(!sapply(out[ch], is.null))) out[[nd]] <- unlist(out[ch]) else done <- FALSE
    }
    if (done) break
  }
  out
}

analyse <- function(ty) {
  f <- file.path(PROJ, "contextual", sprintf("tree_RSV-%s.nwk", substr(ty, 5, 5)))
  if (!file.exists(f)) { cat(sprintf("  [%s] tree missing\n", ty)); return(NULL) }
  phy <- read.tree(f)
  n_tip <- Ntip(phy)
  mel_names <- dat$FieldNoNIV_f[dat$nc_dataset == ty]
  is_mel <- phy$tip.label %in% mel_names
  sup <- suppressWarnings(as.numeric(phy$node.label))

  cat(sprintf("\n[%s] %d tips, %d Melghat, %d internal nodes (%d with support)\n",
              ty, n_tip, sum(is_mel), Nnode(phy), sum(!is.na(sup))))

  des <- descendants(phy)
  node_ids <- (n_tip + 1):(n_tip + Nnode(phy))
  n_desc   <- sapply(node_ids, function(nd) length(des[[nd]]))
  n_mel    <- sapply(node_ids, function(nd) sum(is_mel[des[[nd]]]))
  pure_mel <- n_mel == n_desc

  ## parent lookup
  parent <- integer(n_tip + Nnode(phy))
  parent[phy$edge[, 2]] <- phy$edge[, 1]

  ## maximal pure-Melghat clades: pure, and parent is not pure
  idx_pure <- which(pure_mel)
  maximal <- idx_pure[sapply(idx_pure, function(i) {
    nd <- node_ids[i]; p <- parent[nd]
    p == 0 || !pure_mel[match(p, node_ids)] || is.na(match(p, node_ids))
  })]

  ## Melghat tips not inside any pure clade -> singleton placements
  in_clade <- unique(unlist(lapply(maximal, function(i) des[[node_ids[i]]])))
  singles  <- setdiff(which(is_mel), in_clade)

  ## (b) placement branches = parent node of each maximal clade, and parent of each singleton tip
  place_nodes <- c(sapply(maximal, function(i) parent[node_ids[i]]),
                   parent[singles])
  place_nodes <- place_nodes[place_nodes > n_tip]
  place_sup   <- sup[match(place_nodes, node_ids)]

  ## (a) within-Melghat structure.
  ## Two different things, kept separate because conflating them overstates what
  ## the tree shows:
  ##   grouping  -- the node that DEFINES each maximal Melghat-only clade.
  ##                High values say these genomes group together.
  ##   internal  -- every node strictly INSIDE those clades. These say whether
  ##                the branching order among them is resolved. The first does
  ##                not imply the second.
  group_sup <- sup[maximal][n_desc[maximal] >= 2]
  inside_nodes <- unlist(lapply(maximal, function(i) {
    nd <- node_ids[i]
    d  <- des[[nd]]
    kept <- nid_in <- node_ids[sapply(node_ids, function(x)
      x != nd && all(des[[x]] %in% d))]
    kept
  }))
  inside_nodes <- unique(inside_nodes)
  within_sup <- sup[match(inside_nodes, node_ids)]

  q <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return("none")
    sprintf("n=%d median %.3f, IQR %.3f-%.3f, >=0.95 %.0f%%",
            length(x), median(x), quantile(x, .25), quantile(x, .75),
            100 * mean(x >= 0.95))
  }
  cat("  (c) ALL nodes                      ", q(sup), "\n")
  cat("  (b) Melghat PLACEMENT branches     ", q(place_sup), "\n")
  cat("  (a1) clade-DEFINING nodes (grouping)", q(group_sup), "\n")
  cat("  (a2) nodes INSIDE clades (branching order)", q(within_sup), "\n")

  cat(sprintf("\n  %d maximal Melghat-only clades covering %d of %d Melghat genomes;",
              length(maximal), length(in_clade), sum(is_mel)))
  cat(sprintf(" %d sit singly among contextual genomes\n", length(singles)))

  ## annotate each maximal clade
  meta <- dat %>% filter(nc_dataset == ty) %>%
    select(FieldNoNIV_f, season, clade, code_village, link_status)
  rows <- lapply(seq_along(maximal), function(j) {
    i <- maximal[j]; nd <- node_ids[i]
    tips <- phy$tip.label[des[[nd]]]
    m <- meta %>% filter(FieldNoNIV_f %in% tips)
    lab <- function(x) { t <- sort(table(x), decreasing = TRUE)
      paste(sprintf("%s(%d)", names(t), as.integer(t)), collapse = " ") }
    tibble(rsv_type = ty, clade_id = j, n_genomes = length(tips),
           support = round(sup[i], 3),
           placement_support = round(sup[match(parent[nd], node_ids)], 3),
           seasons = lab(m$season), lineages = lab(m$clade),
           villages = lab(m$code_village),
           n_villages = n_distinct(m$code_village))
  })
  out <- bind_rows(rows) %>% arrange(desc(n_genomes))
  print(as.data.frame(out %>% filter(n_genomes >= 2) %>%
          select(clade_id, n_genomes, support, placement_support,
                 n_villages, seasons, lineages)), right = FALSE)

  ## singletons, with their placement support
  sing <- tibble(rsv_type = ty, clade_id = NA_integer_, n_genomes = 1L,
                 support = NA_real_,
                 placement_support = round(sup[match(parent[singles], node_ids)], 3)) %>%
    bind_cols(meta[match(phy$tip.label[singles], meta$FieldNoNIV_f),
                   c("season", "clade", "code_village")] %>%
                rename(seasons = season, lineages = clade, villages = code_village)) %>%
    mutate(n_villages = 1L)

  list(summary = tibble(rsv_type = ty,
                        all_median = median(sup, na.rm = TRUE),
                        all_ge95 = mean(sup >= 0.95, na.rm = TRUE),
                        place_median = median(place_sup, na.rm = TRUE),
                        place_ge95 = mean(place_sup >= 0.95, na.rm = TRUE),
                        group_median = median(group_sup, na.rm = TRUE),
                        group_ge95 = mean(group_sup >= 0.95, na.rm = TRUE),
                        inside_n = sum(!is.na(within_sup)),
                        inside_median = median(within_sup, na.rm = TRUE),
                        inside_ge95 = mean(within_sup >= 0.95, na.rm = TRUE),
                        n_maximal_clades = length(maximal),
                        n_in_clades = length(in_clade),
                        n_singletons = length(singles)),
       clades = bind_rows(out, sing))
}

res <- lapply(c("RSV-A", "RSV-B"), analyse)
res <- res[!sapply(res, is.null)]
if (length(res)) {
  S <- bind_rows(lapply(res, `[[`, "summary"))
  C <- bind_rows(lapply(res, `[[`, "clades"))
  write_csv(S, file.path(OUT, "16_support_by_relevance.csv"))
  write_csv(C, file.path(OUT, "16_melghat_clades.csv"))
  cat("\n  -> wrote 16_support_by_relevance.csv, 16_melghat_clades.csv\n")
}
cat("== 16 done ==\n")

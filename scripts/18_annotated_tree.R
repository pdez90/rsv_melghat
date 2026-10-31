# 18_annotated_tree.R -- readable annotated view of where Melghat genomes sit.
#
# A 2,312-tip tree is not a figure. This prunes to every Melghat genome plus the
# contextual tips each Melghat clade is actually sister to, then annotates:
# season (colour), lineage, village, and candidate F/G states. SH-like support is
# printed ONLY on branches that bear on Melghat placement or internal structure --
# printing it everywhere is what makes these plots unreadable and invites the
# overall-support summary that does not answer the question.

source(file.path("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(ape)})
cat("\n== 18 annotated tree ==\n")

dat <- readRDS(file.path(OUT, "02_dat.rds"))
CTXN <- 3   # contextual neighbours to retain per Melghat clade

SEASON_COL <- c("2016/17" = "#2a78d6", "2017/18" = "#eb6834", "2018/19" = "#1baf7a",
                "2019/20" = "#eda100", "2021/22" = "#6250d6")

descendants <- function(phy) {
  n_tip <- Ntip(phy); n_node <- Nnode(phy)
  kids <- vector("list", n_tip + n_node)
  for (i in seq_len(nrow(phy$edge))) kids[[phy$edge[i,1]]] <- c(kids[[phy$edge[i,1]]], phy$edge[i,2])
  out <- vector("list", n_tip + n_node)
  for (i in seq_len(n_tip)) out[[i]] <- i
  repeat {
    done <- TRUE
    for (nd in (n_tip+1):(n_tip+n_node)) {
      if (!is.null(out[[nd]])) next
      ch <- kids[[nd]]
      if (all(!sapply(out[ch], is.null))) out[[nd]] <- unlist(out[ch]) else done <- FALSE
    }
    if (done) break
  }
  out
}

build <- function(ty) {
  f <- file.path(PROJ, "contextual", sprintf("tree_RSV-%s.nwk", substr(ty, 5, 5)))
  if (!file.exists(f)) { cat(sprintf("  [%s] tree missing\n", ty)); return(invisible(NULL)) }
  phy <- read.tree(f); n_tip <- Ntip(phy)
  meta <- dat %>% filter(nc_dataset == ty)
  mel  <- intersect(meta$FieldNoNIV_f, phy$tip.label)
  is_mel <- phy$tip.label %in% mel
  sup <- suppressWarnings(as.numeric(phy$node.label))
  des <- descendants(phy); nid <- (n_tip+1):(n_tip+Nnode(phy))
  parent <- integer(n_tip + Nnode(phy)); parent[phy$edge[,2]] <- phy$edge[,1]
  pure <- sapply(nid, function(nd) { d <- des[[nd]]; all(is_mel[d]) })
  maximal <- nid[pure][sapply(nid[pure], function(nd) {
    p <- parent[nd]; p == 0 || !pure[match(p, nid)] || is.na(match(p, nid)) })]

  ## contextual neighbours: non-Melghat tips under each maximal clade's parent
  ctx_keep <- unique(unlist(lapply(maximal, function(nd) {
    p <- parent[nd]; if (p == 0) return(NULL)
    d <- des[[p]]; nm <- d[!is_mel[d]]
    if (!length(nm)) return(NULL)
    phy$tip.label[head(nm, CTXN)]
  })))
  ## singleton Melghat tips keep their own neighbours too
  singles <- setdiff(which(is_mel), unique(unlist(lapply(maximal, function(nd) des[[nd]]))))
  ctx_keep <- unique(c(ctx_keep, unlist(lapply(singles, function(t) {
    p <- parent[t]; d <- des[[p]]; nm <- d[!is_mel[d]]
    phy$tip.label[head(nm, CTXN)] }))))

  keep <- unique(c(mel, ctx_keep))
  sub <- keep.tip(phy, keep)
  subsup <- suppressWarnings(as.numeric(sub$node.label))
  cat(sprintf("  [%s] pruned %d -> %d tips (%d Melghat + %d contextual)\n",
              ty, n_tip, Ntip(sub), length(mel), Ntip(sub) - length(mel)))

  ## tip annotation
  m <- meta[match(sub$tip.label, meta$FieldNoNIV_f), ]
  is_m <- !is.na(m$FieldNoNIV_f)
  # contextual tips: accession_Country_Date_Subgroup
  ctx_lab <- sapply(strsplit(sub$tip.label, "_"), function(p)
    if (length(p) >= 3) sprintf("%s %s", p[2], substr(p[3], 1, 4)) else p[1])
  lab <- ifelse(is_m,
                sprintf("%s | %s | %s | %s", sub$tip.label, m$season, m$clade,
                        sub("^[0-9]+-", "", m$code_village)),
                paste0("  ", ctx_lab))
  col <- ifelse(is_m, SEASON_COL[as.character(m$season)], "#9a9a95")
  col[is.na(col)] <- "#9a9a95"

  ## which internal nodes matter: those with >=1 Melghat descendant and >=2 tips
  dsub <- descendants(sub); nsub <- (Ntip(sub)+1):(Ntip(sub)+Nnode(sub))
  relevant <- sapply(nsub, function(nd) {
    d <- dsub[[nd]]; any(is_m[d]) && length(d) >= 2 })
  shown <- rep(NA_character_, Nnode(sub))
  shown[relevant] <- ifelse(is.na(subsup[relevant]), "",
                            sprintf("%.2f", subsup[relevant]))

  h <- max(7, Ntip(sub) * 0.105)
  png(file.path(OUT, sprintf("18_annotated_tree_%s.png", ty)),
      width = 1500, height = h * 100, res = 130)
  par(mar = c(3, 0.5, 3, 0.5), bg = "#fcfcfb")
  sub$tip.label <- lab
  plot(sub, show.tip.label = TRUE, tip.color = col, cex = 0.42,
       label.offset = 0.00012, edge.width = 0.7, edge.color = "#5c5c58")
  nodelabels(shown, frame = "none", cex = 0.34, col = "#b03a3a", adj = c(1.1, -0.4))
  title(main = sprintf("%s: Melghat genomes in global context (exploratory)", ty),
        cex.main = 0.95, adj = 0, font.main = 2)
  mtext(paste("Tip colour = season. Grey = contextual (country, year).",
              "Red numbers = SH-like local support on Melghat-relevant branches",
              "-- not bootstrap."),
        side = 1, line = 1, cex = 0.55, adj = 0, col = "#5c5c58")
  legend("bottomright", legend = names(SEASON_COL), text.col = SEASON_COL,
         bty = "n", cex = 0.6)
  dev.off()
  cat(sprintf("  -> 18_annotated_tree_%s.png (%.0f x %.0f px)\n", ty, 1500, h*100))
}

invisible(lapply(c("RSV-A", "RSV-B"), build))
cat("== 18 done ==\n")

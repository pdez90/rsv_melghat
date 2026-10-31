# 05_phylogeography.R -- does RSV genetic relatedness track village geography?
#
# Question: within a season, are genomes from the same village more closely
# related than genomes from different villages?
#
# NOTE ON INFERENCE: the medians computed here are DESCRIPTIVE. Pairwise distances
# share genomes, so pair-level tests (e.g. Wilcoxon) are invalid and have been
# removed. The valid test is the genome-label permutation in 08_robustness.R.
# NOTE ON DISTANCE: dist.dna(model="N", pairwise.deletion=TRUE) is an unnormalised
# mismatch count, so a low-coverage genome can look artificially close because
# fewer sites were compared. Shared callable length per pair is reported below. If yes, transmission has a local
# structure; if no, the season looks like one well-mixed epidemic across Melghat.
# Also estimates a minimum number of introductions per season by counting
# genetic clusters.

source(file.path("scripts", "00_setup.R"))
cat("\n== 05 phylogeography ==\n")
set.seed(20261030)

dat <- readRDS(file.path(OUT, "02_dat.rds"))

## ---- village centroids from household GPS -----------------------------------
# The village survey carries no coordinates, so village position is the median
# household GPS point. Raw GPS has data-entry errors (latitude up to 166624), so
# points are first clipped to the Melghat bounding box.
gps <- read_csv(data_path("critical_data_gpsdata.csv"),
                show_col_types = FALSE, guess_max = 5000, progress = FALSE)
BOX <- list(lat = c(20.9, 21.9), lon = c(76.3, 77.6))
gps_ok <- gps %>%
  filter(!is.na(latitude), !is.na(longitude),
         between(latitude, BOX$lat[1], BOX$lat[2]),
         between(longitude, BOX$lon[1], BOX$lon[2]))
cat(sprintf("GPS households: %d total, %d inside the Melghat box (%d dropped as out-of-range)\n",
            nrow(gps), nrow(gps_ok), nrow(gps) - nrow(gps_ok)))

cent <- gps_ok %>% group_by(code_village) %>%
  summarise(lat = median(latitude), lon = median(longitude), households = n(), .groups = "drop")
cat(sprintf("village centroids: %d villages\n", nrow(cent)))

## ---- genomes with a village and an episode ----------------------------------
d <- dat %>% filter(linked, !is.na(code_village), !is.na(clade)) %>%
  left_join(cent, by = c("code_village" = "code_village"))
cat(sprintf("linked genomes: %d; with village centroid: %d (%d villages)\n",
            sum(dat$linked), sum(!is.na(d$lat)), n_distinct(d$code_village[!is.na(d$lat)])))
d <- d %>% filter(!is.na(lat))

## ---- pairwise genetic distance ----------------------------------------------
snp_dist <- function(path) {
  m <- read_fasta_matrix(path)
  m <- m[rownames(m) != "", , drop = FALSE]
  as.matrix(ape::dist.dna(ape::as.DNAbin(tolower(m)), model = "N",
                          pairwise.deletion = TRUE, as.matrix = TRUE))
}
DA <- snp_dist(file.path(PROJ, "nextclade", "RSV-A", "nextclade.aligned.fasta"))
DB <- snp_dist(file.path(PROJ, "nextclade", "RSV-B", "nextclade.aligned.fasta"))
cat(sprintf("SNP distance matrices: RSV-A %d x %d, RSV-B %d x %d\n",
            nrow(DA), ncol(DA), nrow(DB), ncol(DB)))

haversine <- function(lat1, lon1, lat2, lon2) {
  R <- 6371; p <- pi/180
  a <- sin((lat2-lat1)*p/2)^2 + cos(lat1*p)*cos(lat2*p)*sin((lon2-lon1)*p/2)^2
  2*R*asin(pmin(1, sqrt(a)))
}

## ---- per season: same-village vs different-village genetic distance ---------
mantel_perm <- function(G, Km, nperm = 9999) {
  iu <- upper.tri(G)
  obs <- suppressWarnings(cor(G[iu], Km[iu], method = "spearman"))
  if (is.na(obs)) return(c(rho = NA, p = NA))
  n <- nrow(G)
  null <- replicate(nperm, { s <- sample(n); suppressWarnings(cor(G[s, s][iu], Km[iu], method = "spearman")) })
  c(rho = obs, p = (1 + sum(abs(null) >= abs(obs), na.rm = TRUE)) / (nperm + 1))
}

results <- list(); clusters <- list()
for (ty in c("RSV-A", "RSV-B")) {
  D <- if (ty == "RSV-A") DA else DB
  for (ss in sort(unique(d$season[d$nc_dataset == ty]))) {
    sub <- d %>% filter(nc_dataset == ty, season == ss)
    ids <- intersect(sub$FieldNoNIV_f, rownames(D))
    if (length(ids) < 8) next
    sub <- sub %>% filter(FieldNoNIV_f %in% ids) %>% arrange(match(FieldNoNIV_f, ids))
    G <- D[ids, ids]
    Km <- outer(seq_along(ids), seq_along(ids),
                Vectorize(function(i, j) haversine(sub$lat[i], sub$lon[i], sub$lat[j], sub$lon[j])))
    iu <- upper.tri(G)
    same <- outer(sub$code_village, sub$code_village, "==")[iu]
    # Same-lineage control: villages tend to see one lineage, so an apparent
    # same-village effect could simply be a same-lineage effect. Restricting to
    # pairs of the same lineage separates the two.
    samelin <- outer(sub$clade, sub$clade, "==")[iu]
    mt <- mantel_perm(G, Km)
    results[[paste(ty, ss)]] <- tibble(
      rsv_type = ty, season = ss, genomes = length(ids),
      villages = n_distinct(sub$code_village),
      median_snp_overall      = median(G[iu]),
      median_snp_same_village = if (any(same)) median(G[iu][same]) else NA_real_,
      median_snp_diff_village = if (any(!same)) median(G[iu][!same]) else NA_real_,
      n_same_village_pairs    = sum(same),

      # same-lineage pairs only
      median_snp_same_village_samelin = if (any(same & samelin)) median(G[iu][same & samelin]) else NA_real_,
      median_snp_diff_village_samelin = if (any(!same & samelin)) median(G[iu][!same & samelin]) else NA_real_,

      mantel_rho = unname(mt["rho"]), mantel_p = unname(mt["p"]))
    # Connected components at a SNP threshold. These are THRESHOLD-DEFINED
    # GENETIC CLUSTERS, not a minimum number of introductions: one introduction
    # that then diversifies splits into several components, and separate
    # introductions of similar genomes merge into one.
    for (thr in c(3, 5, 10)) {
      adj <- G <= thr
      comp <- rep(NA_integer_, length(ids)); k <- 0L
      for (i in seq_along(ids)) if (is.na(comp[i])) {
        k <- k + 1L; queue <- i
        while (length(queue)) {
          v <- queue[1]; queue <- queue[-1]
          if (is.na(comp[v])) { comp[v] <- k; queue <- c(queue, which(adj[v, ] & is.na(comp))) }
        }
      }
      clusters[[paste(ty, ss, thr)]] <- tibble(rsv_type = ty, season = ss, threshold_snps = thr,
                                               genomes = length(ids), clusters = k,
                                               largest_cluster = max(table(comp)))
    }
  }
}
res <- bind_rows(results); clu <- bind_rows(clusters)

cat("\n-- genetic distance within vs between villages, by season --\n")
print(as.data.frame(res), digits = 3)
cat("\n-- threshold-defined genetic clusters per season --\n")
cat("NOT a count of introductions; see the note in the script header.\n")
print(as.data.frame(clu))

## ---- overall pooled test ----------------------------------------------------
cat("\n-- pooled across seasons (same-season pairs only) --\n")
pooled <- res %>% summarise(
  seasons_tested = n(),
  seasons_same_village_closer = sum(median_snp_same_village < median_snp_diff_village, na.rm = TRUE),
  seasons_mantel_p_lt_05 = sum(mantel_p < 0.05, na.rm = TRUE))
cat("Wilcoxon tests on pairwise distances have been REMOVED: overlapping pairs are\n")
cat("not independent, so their p-values were invalid. Inference is the genome-label\n")
cat("permutation test in 08_robustness.R. Medians below are descriptive only.\n")
print(as.data.frame(pooled))

write_out(res,  "05_phylogeography_by_season.csv")
write_out(clu,  "05_genetic_clusters.csv")
write_out(cent, "05_village_centroids.csv")
cat("== 05 done ==\n")

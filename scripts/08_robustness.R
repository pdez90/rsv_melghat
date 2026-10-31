# 08_robustness.R -- checks prompted by external review of the first results.
#
# 1. Village clustering re-tested with a permutation that respects the
#    non-independence of pairwise comparisons (Wilcoxon on pairs does not).
# 2. Whether same-village pairs are merely close in TIME.
# 3. Whether the Dhodara K272R samples are a genetic cluster or just share one
#    substitution.
# 4. Whether the antigenic-site positions were actually covered by sequence.
# 5. Reinfection intervals and sequence quality.

source(file.path("scripts", "00_setup.R"))
cat("\n== 08 robustness checks ==\n")
set.seed(20261030)

dat <- readRDS(file.path(OUT, "02_dat.rds"))

snp_dist <- function(path) {
  m <- read_fasta_matrix(path); m <- m[rownames(m) != "", , drop = FALSE]
  as.matrix(ape::dist.dna(ape::as.DNAbin(tolower(m)), model = "N",
                          pairwise.deletion = TRUE, as.matrix = TRUE))
}
DA <- snp_dist(file.path(PROJ, "nextclade", "RSV-A", "nextclade.aligned.fasta"))
DB <- snp_dist(file.path(PROJ, "nextclade", "RSV-B", "nextclade.aligned.fasta"))

cent <- read_csv(file.path(OUT, "05_village_centroids.csv"), show_col_types = FALSE)
d <- dat %>% filter(linked, !is.na(code_village), !is.na(clade)) %>%
  left_join(cent, by = "code_village") %>% filter(!is.na(lat))

## ---- 1 & 2. village clustering, permutation test ----------------------------
# The unit of independence is the GENOME, not the pair. Permuting village labels
# across genomes within a season preserves the pair-dependence structure and
# gives a valid null for "same-village pairs are closer".
cat("\n-- village clustering: permutation test (village labels shuffled, 9999x) --\n")
out <- list()
for (ty in c("RSV-A", "RSV-B")) {
  D <- if (ty == "RSV-A") DA else DB
  for (ss in sort(unique(d$season[d$nc_dataset == ty]))) {
    sub <- d %>% filter(nc_dataset == ty, season == ss)
    ids <- intersect(sub$FieldNoNIV_f, rownames(D))
    if (length(ids) < 8) next
    sub <- sub %>% filter(FieldNoNIV_f %in% ids) %>% arrange(match(FieldNoNIV_f, ids))
    G <- D[ids, ids]; iu <- upper.tri(G)
    vill <- sub$code_village
    stat <- function(v) {
      s <- outer(v, v, "==")[iu]
      if (!any(s) || !any(!s)) return(NA_real_)
      median(G[iu][s]) - median(G[iu][!s])
    }
    obs <- stat(vill)
    if (is.na(obs)) next
    null <- replicate(9999, stat(sample(vill)))
    pval <- (1 + sum(null <= obs, na.rm = TRUE)) / (1 + sum(!is.na(null)))
    # time: are same-village pairs also closer in date?
    dt <- abs(outer(as.numeric(sub$swab_date), as.numeric(sub$swab_date), "-"))[iu]
    s <- outer(vill, vill, "==")[iu]
    out[[paste(ty, ss)]] <- tibble(
      rsv_type = ty, season = ss, genomes = length(ids),
      obs_diff_snps = round(obs, 1),
      null_mean = round(mean(null, na.rm = TRUE), 1),
      perm_p = round(pval, 4),
      median_days_apart_same_village = median(dt[s]),
      median_days_apart_diff_village = median(dt[!s]))
  }
}
perm <- bind_rows(out)
print(as.data.frame(perm))
cat("\nobs_diff_snps = median(same-village pairs) - median(different-village pairs);\n")
cat("negative means same-village genomes are more closely related. perm_p is one-sided.\n")

## ---- 3. the Dhodara K272R samples -------------------------------------------
cat("\n-- are the Dhodara K272R genomes a genetic cluster? --\n")
k <- dat %>% filter(nc_escape > 0) %>%
  select(FieldNoNIV_f, studyid, code_village, swab_date, clade, nc_qc_status)
ids <- intersect(k$FieldNoNIV_f, rownames(DB))
Gk <- DB[ids, ids]
print(Gk)
dh <- k$FieldNoNIV_f[k$code_village == "028-Dhodara"]
cat(sprintf("\nDhodara samples (n=%d): %s\n", length(dh), paste(dh, collapse = ", ")))
cat(sprintf("pairwise SNPs within Dhodara: median %.0f, range %.0f-%.0f\n",
            median(Gk[dh, dh][upper.tri(Gk[dh, dh])]), min(Gk[dh, dh][upper.tri(Gk[dh, dh])]),
            max(Gk[dh, dh][upper.tri(Gk[dh, dh])])))
oth <- setdiff(ids, dh)
cat(sprintf("Dhodara vs the other K272R samples: median %.0f SNPs\n", median(Gk[dh, oth])))
# background: all 2018/19 B.D.4.1 pairs
bg <- dat %>% filter(season == "2018/19", clade == "B.D.4.1", linked) %>% pull(FieldNoNIV_f)
bg <- intersect(bg, rownames(DB))
Gb <- DB[bg, bg]
cat(sprintf("background, all 2018/19 B.D.4.1 pairs (n=%d): median %.0f SNPs\n",
            length(bg), median(Gb[upper.tri(Gb)])))
cat(sprintf("collection dates in Dhodara: %s\n",
            paste(format(sort(k$swab_date[k$code_village == "028-Dhodara"]), "%d %b %Y"), collapse = ", ")))

## ---- 4. were the antigenic positions actually sequenced? --------------------
cat("\n-- coverage of the F antigenic sites (X = unknown amino acid) --\n")
SITE_POS <- c(62:69, 196:212, 254:277, 376:389, 46:54, 301:311, 422:471, 55:61, 146:194)
for (ty in c("RSV-A", "RSV-B")) {
  f <- read_fasta_matrix(file.path(PROJ, "nextclade", sub("RSV-", "RSV-", ty),
                                   "nextclade.cds_translation.F.fasta"))
  pos <- SITE_POS[SITE_POS <= ncol(f)]
  unknown <- rowSums(f[, pos, drop = FALSE] == "X")
  cat(sprintf("%s: %d genomes; antigenic-site residues unknown (X): %d genomes affected, max %d of %d residues\n",
              ty, nrow(f), sum(unknown > 0), max(unknown), length(pos)))
  if (any(unknown > 0)) print(sort(unknown[unknown > 0], decreasing = TRUE))
}

## ---- 5. reinfections ---------------------------------------------------------
cat("\n-- reinfection intervals and quality --\n")
re <- dat %>% group_by(studyid) %>% filter(n() > 1) %>% arrange(studyid, swab_date) %>%
  mutate(interval_days = as.integer(swab_date - lag(swab_date))) %>%
  select(studyid, FieldNoNIV_f, swab_date, interval_days, nc_dataset, clade,
         nc_qc_status, nc_coverage, CTvalue_f) %>% ungroup()
print(as.data.frame(re))

write_out(perm, "08_village_clustering_permutation.csv")
write_out(re,   "08_reinfection_detail.csv")
cat("== 08 done ==\n")

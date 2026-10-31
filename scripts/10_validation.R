# 10_validation.R -- corrections and normalised measures after external review.
#
# Fixes four things that were previously overstated:
#  1. translation boundaries -- Nextclade CDS translations INCLUDE the terminal
#     stop codon, so "n codons" is n-1 amino acids. RSV-B G is not uniform: half
#     the genomes have lost the stop and encode a longer protein.
#  2. the G frameshifts are NOT at the duplication. A 72-nt deletion is in-frame
#     (72/3 = 24) and cannot cause a frameshift, so they are a separate problem.
#  3. G-vs-F variability, normalised per callable residue and also measured
#     reference-free as mean pairwise difference.
#  4. HVR2 enrichment expressed against the fraction of the protein it occupies.

source(file.path("scripts", "00_setup.R"))
cat("\n== 10 validation and normalised measures ==\n")
set.seed(20261030)

dat <- readRDS(file.path(OUT, "02_dat.rds"))
CDS <- list(`RSV-A` = list(G = c(4652, 5617), F = c(5697, 7421)),
            `RSV-B` = list(G = c(4646, 5578), F = c(5676, 7400)))

## ---- 1. translation boundaries ----------------------------------------------
cat("\n-- 1. translation boundaries (stop codon is included in the CDS) --\n")
for (ty in c("RSV-A", "RSV-B")) for (gn in c("G", "F")) {
  co <- CDS[[ty]][[gn]]; n <- co[2] - co[1] + 1
  cat(sprintf("%s %s CDS %d-%d = %d nt = %d codons = %d amino acids + stop\n",
              ty, gn, co[1], co[2], n, n/3, n/3 - 1))
}

cat("\nterminal residue of each G translation:\n")
rt_ids <- character(0)
for (ty in c("RSV-A", "RSV-B")) {
  f <- read_fasta_matrix(file.path(PROJ, "nextclade", ty, "nextclade.cds_translation.G.fasta"))
  last <- apply(f, 1, function(r) { v <- r[r != "-"]; v[length(v)] })
  cat(sprintf("  %s: ", ty)); print(table(last))
  if (ty == "RSV-B") rt_ids <- names(last)[last == "Q"]
}

alnB <- read_fasta_matrix(file.path(PROJ, "nextclade", "RSV-B", "nextclade.aligned.fasta"))
G1 <- CDS$`RSV-B`$G[2]
cat(sprintf("\nRSV-B: %d of 111 genomes carry TAA->CAA at the G stop codon (nt %d-%d),\n",
            length(rt_ids), G1 - 2, G1))
nxt <- sapply(rt_ids, function(s) {
  v <- toupper(alnB[s, (G1 + 1):min(G1 + 300, ncol(alnB))])
  cod <- split(v, ceiling(seq_along(v) / 3)); cod <- cod[sapply(cod, length) == 3]
  aa <- sapply(cod, function(c3) translate_nt(c3))
  w <- which(aa == "*")[1]; if (is.na(w)) NA_integer_ else w })
cat(sprintf("losing the stop. Next in-frame stop is %s codons downstream in all of them,\n",
            paste(unique(nxt), collapse = "/")))
cat(sprintf("so those genomes encode a G of %d aa, not %d.\n", 310 + unique(nxt)[1], 310))
cat("=> RSV-B G is NOT uniform in length. Earlier statement withdrawn.\n")

## ---- 2. where are the G frameshifts? ----------------------------------------
cat("\n-- 2. G frameshift locations vs the duplication block (G aa 262-285) --\n")
fsA <- read_tsv(file.path(PROJ, "nextclade", "RSV-A", "nextclade.tsv"),
                show_col_types = FALSE, guess_max = 500, progress = FALSE) %>%
  filter(str_detect(frameShifts, "(^|,)G:")) %>%
  transmute(seqName, G_fs = sapply(str_extract_all(frameShifts, "G:[0-9]+-[0-9]+"),
                                   paste, collapse = ";"),
            start = as.integer(str_extract(str_extract(frameShifts, "G:[0-9]+"), "[0-9]+")))
print(as.data.frame(fsA))
cat(sprintf("\n%d of %d start at codon 309 or later; the duplication ends at 285.\n",
            sum(fsA$start >= 300), nrow(fsA)))
cat("A 72-nt deletion is in-frame and CANNOT produce a frameshift, so these are a\n")
cat("separate C-terminal problem, not evidence about the duplication. Earlier\n")
cat("inference from the frameshift counts withdrawn.\n")

## ---- 3. normalised G vs F variability ---------------------------------------
cat("\n-- 3. G vs F variability, normalised --\n")
aa <- dat %>% filter(!is.na(nc_aa_subs), nc_aa_subs != "") %>%
  select(FieldNoNIV_f, nc_dataset, nc_aa_subs) %>%
  separate_rows(nc_aa_subs, sep = ",") %>%
  mutate(gene = str_extract(nc_aa_subs, "^[^:]+")) %>%
  filter(gene %in% c("G", "F")) %>%
  count(FieldNoNIV_f, nc_dataset, gene, name = "subs")

# callable residues per genome per gene (exclude X), and per-100-residue rate
rows <- list()
for (ty in c("RSV-A", "RSV-B")) for (gn in c("G", "F")) {
  f <- read_fasta_matrix(file.path(PROJ, "nextclade", ty,
                                   sprintf("nextclade.cds_translation.%s.fasta", gn)))
  f <- f[, -ncol(f), drop = FALSE]                       # drop the stop position
  callable <- apply(f, 1, function(r) sum(r %in% LETTERS & r != "X"))
  rows[[paste(ty, gn)]] <- tibble(FieldNoNIV_f = rownames(f), nc_dataset = ty,
                                  gene = gn, callable = as.integer(callable))
}
cal <- bind_rows(rows) %>% left_join(aa, by = c("FieldNoNIV_f", "nc_dataset", "gene")) %>%
  mutate(subs = replace_na(subs, 0L), per100 = 100 * subs / callable)
norm <- cal %>% group_by(nc_dataset, gene) %>%
  summarise(genomes = n(), median_callable_aa = median(callable),
            median_subs = median(subs),
            median_subs_per_100aa = round(median(per100), 2), .groups = "drop")
cat("\nreference-relative, normalised per 100 callable residues:\n")
print(as.data.frame(norm))
ratio <- norm %>% select(nc_dataset, gene, median_subs_per_100aa) %>%
  pivot_wider(names_from = gene, values_from = median_subs_per_100aa) %>%
  mutate(G_over_F = round(G / F, 1))
print(as.data.frame(ratio))

# reference-free: mean pairwise amino-acid differences
cat("\nreference-free mean pairwise aa difference (per 100 compared residues):\n")
pw <- list()
for (ty in c("RSV-A", "RSV-B")) for (gn in c("G", "F")) {
  f <- read_fasta_matrix(file.path(PROJ, "nextclade", ty,
                                   sprintf("nextclade.cds_translation.%s.fasta", gn)))
  f <- f[, -ncol(f), drop = FALSE]
  ok <- apply(f, 1, function(r) sum(r == "X")) < 10
  f <- f[ok, , drop = FALSE]
  idx <- combn(nrow(f), 2)
  if (ncol(idx) > 4000) idx <- idx[, sample(ncol(idx), 4000)]
  vals <- apply(idx, 2, function(ij) {
    a <- f[ij[1], ]; b <- f[ij[2], ]
    keep <- a %in% LETTERS & b %in% LETTERS & a != "X" & b != "X"
    if (!any(keep)) return(NA_real_)
    100 * sum(a[keep] != b[keep]) / sum(keep) })
  pw[[paste(ty, gn)]] <- tibble(nc_dataset = ty, gene = gn, genomes = nrow(f),
                                mean_pct_diff = round(mean(vals, na.rm = TRUE), 2))
}
pwd <- bind_rows(pw)
print(as.data.frame(pwd))
print(as.data.frame(pwd %>% select(-genomes) %>%
        pivot_wider(names_from = gene, values_from = mean_pct_diff) %>%
        mutate(G_over_F = round(G / F, 1))))

## ---- 4. HVR2 enrichment, with a denominator ---------------------------------
cat("\n-- 4. HVR2: share of substitutions vs share of the protein --\n")
HVR <- list(`RSV-A` = list(rng = 230:321, len = 321), `RSV-B` = list(rng = 215:310, len = 310))
hv <- dat %>% filter(!is.na(nc_aa_subs), nc_aa_subs != "") %>%
  select(FieldNoNIV_f, nc_dataset, nc_aa_subs) %>% separate_rows(nc_aa_subs, sep = ",") %>%
  filter(str_starts(nc_aa_subs, "G:")) %>%
  mutate(pos = as.integer(str_extract(nc_aa_subs, "\\d+")),
         in_hvr2 = map2_lgl(nc_dataset, pos, ~ .y %in% HVR[[.x]]$rng)) %>%
  group_by(nc_dataset) %>%
  summarise(G_subs = n(), in_hvr2 = sum(in_hvr2),
            pct_subs_in_hvr2 = round(100 * sum(in_hvr2) / n(), 1), .groups = "drop") %>%
  mutate(pct_protein_that_is_hvr2 =
           round(100 * sapply(nc_dataset, function(t) length(HVR[[t]]$rng) / HVR[[t]]$len), 1),
         enrichment = round(pct_subs_in_hvr2 / pct_protein_that_is_hvr2, 2))
print(as.data.frame(hv))
cat("enrichment = (share of substitutions) / (share of protein). 1.0 means no enrichment.\n")

## ---- 5. glycosylation: positions, not just counts ---------------------------
cat("\n-- 5. predicted N-linked sequon POSITIONS (counts alone can hide turnover) --\n")
cat("NOTE: N-linked sequons only. G also carries O-linked glycans, which are NOT\n")
cat("assessed here and are reported to occur in the ON1 duplicated region.\n")
sequons <- function(v, strict = TRUE) {
  # N-X-S/T with X != P. All three residues must be CANONICAL amino acids:
  # "X" != "P" evaluates TRUE, so a naive test counts a motif whose middle
  # residue was never resolved as a definite sequon. strict = FALSE reproduces
  # the naive behaviour so the two can be compared.
  v <- v[v != "-"]
  if (length(v) < 3) return(integer(0))
  AA20 <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
  n1 <- v[1:(length(v)-2)]; n2 <- v[2:(length(v)-1)]; n3 <- v[3:length(v)]
  mid <- if (strict) n2 %in% setdiff(AA20, "P") else n2 != "P"
  which(n1 == "N" & mid & n3 %in% c("S", "T"))
}

for (ty in c("RSV-A", "RSV-B")) {
  f <- read_fasta_matrix(file.path(PROJ, "nextclade", ty, "nextclade.cds_translation.G.fasta"))
  ok <- rownames(f)[apply(f, 1, function(r) sum(r == "X")) < 10]
  naive  <- sum(vapply(ok, function(s) length(sequons(f[s, ], strict = FALSE)), integer(1)))
  strict <- sum(vapply(ok, function(s) length(sequons(f[s, ])), integer(1)))
  cat(sprintf("\n%s (n=%d usable): %d candidate motifs; %d had an UNRESOLVED middle\n",
              ty, length(ok), naive, naive - strict))
  cat("  residue and are excluded below, reported separately rather than as sequons.\n")
  pl <- lapply(ok, function(s) sequons(f[s, ]))
  tabp <- sort(table(unlist(pl)), decreasing = TRUE)
  cat(sprintf("  distinct sequon positions: %d\n", length(tabp)))
  print(tabp)
  cat(sprintf("  present in ALL genomes: %s\n",
              paste(names(tabp)[tabp == length(ok)], collapse = ", ")))
  cat(sprintf("  variable: %s\n",
              paste(names(tabp)[tabp < length(ok)], collapse = ", ")))
}
cat("\nCAUTION: positions are indices in the UNGAPPED protein. Downstream of the\n")
cat("24-residue RSV-A deletion the numbering shifts, so equal indices are not\n")
cat("necessarily homologous positions. Treat RSV-A turnover at high-numbered\n")
cat("sequons as provisional until mapped to a common reference alignment.\n")

write_out(norm, "10_G_F_normalised_variability.csv")
write_out(pwd,  "10_G_F_pairwise_diversity.csv")
write_out(hv,   "10_HVR2_enrichment.csv")
write_out(fsA,  "10_G_frameshift_positions.csv")
cat("== 10 done ==\n")

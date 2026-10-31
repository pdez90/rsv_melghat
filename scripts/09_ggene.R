# 09_ggene.R -- variation in the RSV attachment (G) glycoprotein.
#
# G is the most variable RSV protein and the classical basis of genotyping. It is
# NOT a target of the licensed monoclonals or the maternal vaccine (those target
# F), so this analysis describes viral diversity and possible immune adaptation
# rather than product susceptibility.
#
# Everything is done SEPARATELY for RSV-A and RSV-B, because the two subgroups
# differ so much in G that pooling them would swamp any within-subgroup signal.
#
# IMPORTANT: amino-acid variation is reported against each subgroup's own
# Nextclade reference (RSV-A hRSV/A/England/397/2017, G = 322 aa;
# RSV-B hRSV/B/Australia/VIC-RCH056/2019, G = 311 aa). A substitution means
# "differs from that reference", not "escapes antibodies". Glycosylation sites
# are sequence-motif PREDICTIONS (N-X-S/T, X != P), not evidence that sugar is
# attached.

source(file.path("scripts", "00_setup.R"))
cat("\n== 09 G gene ==\n")
set.seed(20261030)

dat <- readRDS(file.path(OUT, "02_dat.rds"))

GENE <- list(
  `RSV-A` = list(dir = "RSV-A", cds_start = 4652, ref_len = 322, hvr2 = 230:322),
  `RSV-B` = list(dir = "RSV-B", cds_start = 4646, ref_len = 311, hvr2 = 215:311))

## ---- 1. G protein length: the duplication polymorphism ----------------------
cat("\n-- G protein length --\n")
lens <- list()
for (ty in names(GENE)) {
  f <- read_fasta_matrix(file.path(PROJ, "nextclade", GENE[[ty]]$dir,
                                   "nextclade.cds_translation.G.fasta"))
  L <- apply(f, 1, function(r) sum(r != "-"))
  lens[[ty]] <- tibble(FieldNoNIV_f = names(L), G_len = as.integer(L), rsv_sub = ty)
  cat(sprintf("%s (reference %d translated positions, = %d aa + stop): ",
              ty, GENE[[ty]]$ref_len, GENE[[ty]]$ref_len - 1))
  print(table(L))
}
glen <- bind_rows(lens)
cat("\nNOTE lengths above are TRANSLATED POSITIONS and include the terminal stop.\n")
cat("RSV-B looks uniform here ONLY because the translation window is fixed at the\n")
cat("reference CDS: 53 of 111 RSV-B genomes have lost the stop codon and encode a\n")
cat("longer G. See 10_validation.R section 1.\n")

# RSV-A splits into two length classes 24 aa apart -- the size of the ON1
# 72-nucleotide duplication in the C-terminal hypervariable region.
glen <- glen %>% mutate(
  G_dup = case_when(rsv_sub == "RSV-A" & G_len >= 310 ~ "72-nt duplication present",
                    rsv_sub == "RSV-A" & G_len <  310 ~ "72-nt duplication absent",
                    TRUE ~ "n/a (RSV-B uniform)"))
d <- dat %>% left_join(glen, by = "FieldNoNIV_f")

cat("\n-- RSV-A duplication status by lineage and season --\n")
print(with(d %>% filter(rsv_sub == "RSV-A"), table(clade, G_dup)))
print(with(d %>% filter(rsv_sub == "RSV-A"), table(season, G_dup)))

## ---- 2. is the 'absence' real, or an assembly artifact? ---------------------
# Tandem repeats are a classic failure mode for consensus assembly: a collapsed
# duplication looks exactly like a clean deletion. Three checks.
cat("\n-- is the duplication loss real? --\n")
alnA <- read_fasta_matrix(file.path(PROJ, "nextclade", "RSV-A", "nextclade.aligned.fasta"))
nt <- GENE$`RSV-A`$cds_start + (262 - 1) * 3 + 0:(24 * 3 - 1)
ab <- d$FieldNoNIV_f[d$G_dup == "72-nt duplication absent" & !is.na(d$G_dup)]
blk <- t(sapply(ab, function(s) { v <- alnA[s, nt]
  c(gap = sum(v == "-"), N = sum(v == "N"), base = sum(v %in% c("A","C","G","T"))) }))
cat(sprintf("(a) block at nt %d-%d in the %d 'absent' genomes: median gap %d, median N %d\n",
            min(nt), max(nt), length(ab), median(blk[, "gap"]), median(blk[, "N"])))
flankN <- sapply(ab, function(s) sum(alnA[s, c((min(nt)-30):(min(nt)-1),
                                               (max(nt)+1):(max(nt)+30))] == "N"))
cat(sprintf("    flanking 30 nt either side: median %d N -> not a coverage dropout\n", median(flankN)))

DA <- as.matrix(ape::dist.dna(ape::as.DNAbin(tolower(alnA)), model = "N",
                              pairwise.deletion = TRUE, as.matrix = TRUE))
ad3 <- intersect(d$FieldNoNIV_f[d$clade == "A.D.3"], rownames(DA))
g1 <- intersect(ad3, ab); g2 <- setdiff(ad3, ab)
if (length(g1) > 2 && length(g2) > 2) {
  w1 <- median(DA[g1,g1][upper.tri(DA[g1,g1])]); w2 <- median(DA[g2,g2][upper.tri(DA[g2,g2])])
  b  <- median(DA[g1,g2]); obs <- b - mean(c(w1, w2))
  all3 <- c(g1, g2); lab <- c(rep("a", length(g1)), rep("b", length(g2)))
  null <- replicate(4999, { s <- sample(lab); x <- all3[s=="a"]; y <- all3[s=="b"]
    median(DA[x,y]) - mean(c(median(DA[x,x][upper.tri(DA[x,x])]),
                             median(DA[y,y][upper.tri(DA[y,y])]))) })
  cat(sprintf("(b) within lineage A.D.3 (%d absent vs %d present): between %d SNPs vs within %d/%d;\n",
              length(g1), length(g2), b, w1, w2))
  cat(sprintf("    between-minus-within %.1f, permutation p = %.4f -> phylogenetically coherent\n",
              obs, (1 + sum(null >= obs)) / 5000))
}
fsG <- read_tsv(file.path(PROJ, "nextclade", "RSV-A", "nextclade.tsv"),
                show_col_types = FALSE, guess_max = 500, progress = FALSE) %>%
  filter(str_detect(frameShifts, "(^|,)G:")) %>% pull(seqName)
cat(sprintf("(c) %d genomes carry a G frameshift; %d of them are in the 'absent' group.\n",
            length(fsG), sum(fsG %in% ab)))
cat("    This is NOT evidence about the duplication: 72 nt is divisible by 3, so the\n")
cat("    deletion is in-frame and cannot cause a frameshift, and 8 of the 10 frameshifts\n")
cat("    start at codon 309 -- past the duplication, which ends at 285. Separate problem.\n")
cat("    See 10_validation.R section 2.\n")
cat("\n    VERDICT: unresolved. Neither a real deletion nor a collapsed tandem repeat\n")
cat("    need produce Ns, so the consensus alignment cannot distinguish them. Decisive\n")
cat("    evidence would be reads spanning the repeat junctions, or targeted PCR.\n")
cat("    Report as 'an apparent 72-nt length difference in the consensus alignment'.\n")

## ---- 3. amino-acid variation in G, by subgroup ------------------------------
cat("\n-- amino-acid substitutions in G (vs each subgroup's own reference) --\n")
aa <- dat %>% filter(!is.na(nc_aa_subs), nc_aa_subs != "") %>%
  select(FieldNoNIV_f, season, nc_dataset, clade, nc_aa_subs) %>%
  separate_rows(nc_aa_subs, sep = ",") %>%
  mutate(gene = str_extract(nc_aa_subs, "^[^:]+"),
         sub  = str_remove(nc_aa_subs, "^[^:]+:"),
         pos  = as.integer(str_extract(sub, "\\d+"))) %>%
  filter(gene %in% c("G", "F"), !is.na(pos))

gf <- aa %>% count(FieldNoNIV_f, nc_dataset, gene, name = "subs") %>%
  group_by(nc_dataset, gene) %>%
  summarise(genomes_with_subs = n(), median_subs = median(subs),
            max_subs = max(subs), .groups = "drop")
cat("\nsubstitution burden, G vs F:\n"); print(as.data.frame(gf))
cat("\nRaw counts only -- G is about half the length of F, so these are not comparable.\n")
cat("See 10_validation.R section 3 for rates per 100 callable residues and for a\n")
cat("reference-free pairwise measure (G exceeds F by 5- to 27-fold depending on\n")
cat("subgroup and measure).\n")

cat("\n-- most variable G positions, by subgroup --\n")
topg <- aa %>% filter(gene == "G") %>% count(nc_dataset, sub, pos, name = "genomes") %>%
  group_by(nc_dataset) %>% slice_max(genomes, n = 12, with_ties = FALSE) %>%
  arrange(nc_dataset, desc(genomes)) %>% ungroup()
print(as.data.frame(topg))
cat("NOTE RSV-B '*311Q' is loss of the terminal stop codon, not an ordinary\n")
cat("substitution -- those genomes encode a longer G (10_validation.R section 1).\n")

cat("\n-- G substitutions inside HVR2 (share only; for enrichment vs the share of the\n")
cat("   protein HVR2 occupies, see 10_validation.R section 4) --\n")
hvr <- aa %>% filter(gene == "G") %>%
  mutate(in_hvr2 = map2_lgl(nc_dataset, pos, ~ .y %in% GENE[[.x]]$hvr2)) %>%
  group_by(nc_dataset) %>%
  summarise(all_G_subs = n(), n_in_hvr2 = sum(in_hvr2),
            pct_in_hvr2 = round(100 * sum(in_hvr2) / n(), 1), .groups = "drop")
print(as.data.frame(hvr))

## ---- 4. predicted N-glycosylation sites -------------------------------------
# Sequon N-X-S/T with X != P. A PREDICTION from sequence only.
cat("\n-- predicted N-glycosylation sequons in G (N-X-S/T, X != P) --\n")
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

gly <- list()
for (ty in names(GENE)) {
  f <- read_fasta_matrix(file.path(PROJ, "nextclade", GENE[[ty]]$dir,
                                   "nextclade.cds_translation.G.fasta"))
  usable <- rownames(f)[apply(f, 1, function(r) sum(r == "X")) < 10]
  n <- sapply(usable, function(s) length(sequons(f[s, ])))
  gly[[ty]] <- tibble(FieldNoNIV_f = usable, n_sequons = as.integer(n), rsv_sub = ty)
  cat(sprintf("%s (n=%d usable): median %d sequons, range %d-%d\n",
              ty, length(usable), median(n), min(n), max(n)))
}
glyd <- bind_rows(gly) %>% left_join(dat %>% select(FieldNoNIV_f, season, clade), by = "FieldNoNIV_f")
cat("\nsequon count by lineage:\n")
print(glyd %>% group_by(rsv_sub, clade) %>%
        summarise(n = n(), median_sequons = median(n_sequons),
                  range = paste0(min(n_sequons), "-", max(n_sequons)), .groups = "drop") %>%
        as.data.frame())
cat("\nsequon count by RSV-A duplication status:\n")
print(glyd %>% filter(rsv_sub == "RSV-A") %>%
        left_join(glen %>% select(FieldNoNIV_f, G_dup), by = "FieldNoNIV_f") %>%
        group_by(G_dup) %>% summarise(n = n(), median_sequons = median(n_sequons),
                                      .groups = "drop") %>% as.data.frame())

cat("\nCounts alone can hide positional turnover, and they do here: see\n")
cat("10_validation.R section 5 for the sequon POSITIONS, which vary between genomes.\n")
cat("N-linked sequons only; O-linked glycosylation is not assessed anywhere.\n")

write_out(glen %>% left_join(dat %>% select(FieldNoNIV_f, season, clade), by = "FieldNoNIV_f"),
          "09_G_length_duplication.csv")
write_out(aa %>% filter(gene == "G"), "09_G_substitutions.csv")
write_out(topg, "09_G_top_positions.csv")
write_out(glyd, "09_G_glycosylation.csv")
cat("== 09 done ==\n")

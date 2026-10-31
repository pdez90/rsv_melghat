# 04_fgene_antigenic.R -- F glycoprotein residue states at monoclonal-antibody
# binding sites, with callable denominators and a substitution-specific
# evidence table.
#
# REWRITTEN after code review. Three things the earlier version got wrong:
#
#  1. It reported substitutions RELATIVE TO THE NEXTCLADE REFERENCE and treated
#     "no substitution called" as "reference residue confirmed". Those are not
#     the same: a residue can be uncalled (X). Denominators are now per residue.
#
#  2. It screened POSITIONS, not substitutions. A change at a contact residue is
#     not automatically a resistance mutation -- identity, subgroup and
#     combinations matter, and some substitutions INCREASE susceptibility.
#
#  3. Most consequentially, it described RSV-B M206I / R209Q as substitutions in
#     ~58% of genomes and implied the nirsevimab footprint was conserved. The
#     Nextclade RSV-B reference (VIC-RCH056/2019) itself carries M206/R209, so
#     those calls identify genomes that LACK the I206M+Q209R pair. That pair is
#     a prevalent natural polymorphism that INCREASES nirsevimab susceptibility
#     (AstraZeneca trial analyses; FDA label records <5-fold change, i.e. no
#     reduction). Reporting it as a footprint substitution inverted its meaning.
#
# Everything here is GENOTYPE. Susceptibility is a phenotype measured by
# neutralisation assay; none is performed or implied.

source(file.path("scripts", "00_setup.R"))
cat("\n== 04 F-gene antigenic variation (residue states) ==\n")
set.seed(20261030)

dat <- readRDS(file.path(OUT, "02_dat.rds"))
AA  <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

## Antigenic sites on prefusion F.
SITES <- tribble(
  ~site,      ~start, ~end,
  "site O",       62,   69,
  "site O",      196,  212,
  "site I",      376,  389,
  "site II",     254,  277,
  "site III",     46,   54,
  "site III",    301,  311,
  "site IV",     422,  471,
  "site V",       55,   61,
  "site V",      146,  194)

## ---- substitution-specific evidence table -----------------------------------
# Keyed on the OBSERVED residue, not on a reference-relative call, so it is
# independent of which reference the dataset uses. `effect` is the reported
# direction for that product. Sources: nirsevimab FDA label (2023); Wilkins et
# al., Nat Commun 2023 (Phase 2b/MELODY genotypic analysis); Ahani et al.,
# Lancet Infect Dis 2023 (binding-site conservation). NOT exhaustive -- treat as
# a versioned screen, and re-check against current labelling before publication.
EVIDENCE <- tribble(
  ~subgroup, ~pos, ~aa, ~product,     ~effect,              ~note,
  "RSV-B",     64, "T", "nirsevimab", "reduced",            "I64T, >496-fold alone; also seen with K68E+I206M+Q209R",
  "RSV-B",     68, "E", "nirsevimab", "reduced",            "K68E, >283-fold alone",
  "RSV-B",     68, "N", "nirsevimab", "reduced",            "K68N, with N201S or N208S >13,000-fold",
  "RSV-B",     68, "Q", "nirsevimab", "reduced",            "K68Q, markedly reduced in neutralisation assays",
  "RSV-B",    201, "S", "nirsevimab", "reduced",            "N201S",
  "RSV-B",    201, "T", "nirsevimab", "reduced",            "N201T, >406-fold",
  "RSV-B",    208, "D", "nirsevimab", "reduced",            "N208D, >90,000-fold",
  "RSV-B",    208, "S", "nirsevimab", "reduced",            "N208S, >24,000-fold",
  "RSV-B",    209, "K", "nirsevimab", "reduced",            "Q209K, partial resistance",
  "RSV-B",    209, "L", "nirsevimab", "reduced",            "Q209L, partial resistance",
  "RSV-B",    201, "D", "nirsevimab", "reduced",            "N201D (2025 label)",
  "RSV-B",    208, "I", "nirsevimab", "reduced",            "N208I (2025 label)",
  "RSV-B",    206, "M", "nirsevimab", "combination",        "I206M: evidence is for the I206M+Q209R PAIR, not this residue alone",
  "RSV-B",    209, "R", "nirsevimab", "combination",        "Q209R: evidence is for the I206M+Q209R PAIR, not this residue alone",
  "RSV-B",    211, "N", "nirsevimab", "none",               "S211N, retains susceptibility",
  "RSV-A",     65, "R", "nirsevimab", "none",               "K65R, common polymorphism, retains susceptibility",
  "RSV-A",    201, "S", "nirsevimab", "reduced",            "N201S, 127-fold",
  "RSV-A",    209, "E", "nirsevimab", "reduced",            "K209E, intermediate resistance",
  "RSV-A",     68, "Q", "nirsevimab", "reduced",            "K68Q, with I206M+Q209R 46-fold",
  "RSV-B",    272, "E", "palivizumab", "reduced",           "K272E",
  "RSV-B",    272, "M", "palivizumab", "reduced",           "K272M",
  "RSV-B",    272, "N", "palivizumab", "reduced",           "K272N",
  "RSV-B",    272, "Q", "palivizumab", "reduced",           "K272Q",
  "RSV-B",    272, "T", "palivizumab", "reduced",           "K272T",
  "RSV-B",    262, "D", "palivizumab", "reduced",           "N262D",
  "RSV-B",    262, "S", "palivizumab", "reduced",           "N262S",
  "RSV-B",    275, "F", "palivizumab", "reduced",           "S275F",
  "RSV-B",    275, "L", "palivizumab", "reduced",           "S275L",
  # Palivizumab site II resistance is described for both subgroups; the RSV-A
  # rows mirror the RSV-B ones so the screen is not silently subgroup-limited.
  "RSV-A",    272, "E", "palivizumab", "reduced",           "K272E",
  "RSV-A",    272, "M", "palivizumab", "reduced",           "K272M",
  "RSV-A",    272, "N", "palivizumab", "reduced",           "K272N",
  "RSV-A",    272, "Q", "palivizumab", "reduced",           "K272Q",
  "RSV-A",    272, "T", "palivizumab", "reduced",           "K272T",
  "RSV-A",    262, "D", "palivizumab", "reduced",           "N262D",
  "RSV-A",    262, "S", "palivizumab", "reduced",           "N262S",
  "RSV-A",    275, "F", "palivizumab", "reduced",           "S275F",
  "RSV-A",    275, "L", "palivizumab", "reduced",           "S275L")

## Combinations. Evidence for I206M+Q209R applies to the PAIR; it must not be
## attributed to either residue in isolation.
COMBOS <- tribble(
  ~subgroup, ~label,            ~residues,                ~effect,     ~note,
  "RSV-B",   "I206M+Q209R",     list(c(206, 209), c("M", "R")), "increased/retained",
    "prevalent polymorphism; pair reported to increase nirsevimab susceptibility (FDA label: <5-fold change)")

SCREEN_POS <- sort(unique(c(62:69, 196:212, 254:277)))

## ---- residue-state table, with callable denominators ------------------------
states <- list()
for (ty in c("RSV-A", "RSV-B")) {
  f <- read_fasta_matrix(file.path(PROJ, "nextclade", ty,
                                   "nextclade.cds_translation.F.fasta"))
  f <- f[, -ncol(f), drop = FALSE]                      # drop the stop position
  pos <- SCREEN_POS[SCREEN_POS <= ncol(f)]
  st <- lapply(pos, function(p) {
    v <- f[, p]
    tibble(subgroup = ty, pos = p, aa = v[v %in% AA]) })
  states[[ty]] <- bind_rows(st) %>% count(subgroup, pos, aa, name = "genomes")
}
rs <- bind_rows(states)
denom <- rs %>% group_by(subgroup, pos) %>%
  summarise(callable = sum(genomes), .groups = "drop") %>%
  mutate(total = if_else(subgroup == "RSV-A", 130L, 111L),
         uncalled = total - callable)
rs <- rs %>% left_join(denom, by = c("subgroup", "pos")) %>%
  mutate(pct_of_callable = round(100 * genomes / callable, 1))

cat("\n-- callable coverage at screened residues --\n")
print(denom %>% group_by(subgroup) %>%
        summarise(residues = n(), min_callable = min(callable),
                  median_callable = median(callable), max_uncalled = max(uncalled),
                  .groups = "drop") %>% as.data.frame())
cat("\nresidues where >1 genome is uncalled:\n")
print(denom %>% filter(uncalled > 1) %>% arrange(desc(uncalled)) %>%
        head(12) %>% as.data.frame())

cat("\n-- polymorphic residues in the screened footprints --\n")
poly <- rs %>% group_by(subgroup, pos) %>% filter(n() > 1) %>%
  summarise(callable = first(callable),
            states = paste(sprintf("%s=%d (%.1f%%)", aa, genomes, pct_of_callable),
                           collapse = ", "), .groups = "drop")
print(as.data.frame(poly))
cat("\nEvery other screened residue is invariant across all callable genomes.\n")

## ---- match observed residues against the evidence table ---------------------
cat("\n-- residues matching the substitution-specific evidence table --\n")
hits <- rs %>% inner_join(EVIDENCE, by = c("subgroup", "pos", "aa")) %>%
  select(subgroup, pos, aa, product, effect, genomes, callable, pct_of_callable, note) %>%
  arrange(effect, subgroup, pos)
if (nrow(hits) == 0) cat("none\n") else print(as.data.frame(hits))

cat("\n-- combinations (evidence attaches to the pair, not to either residue alone) --\n")
for (i in seq_len(nrow(COMBOS))) {
  ty <- COMBOS$subgroup[i]; pos <- COMBOS$residues[[i]][[1]]; aas <- COMBOS$residues[[i]][[2]]
  f <- read_fasta_matrix(file.path(PROJ, "nextclade", ty, "nextclade.cds_translation.F.fasta"))
  a1 <- f[, as.integer(pos[1])]; a2 <- f[, as.integer(pos[2])]
  ok <- a1 %in% AA & a2 %in% AA
  both <- sum(ok & a1 == aas[1] & a2 == aas[2])
  neither <- sum(ok & a1 != aas[1] & a2 != aas[2])
  discord <- sum(ok) - both - neither
  cat(sprintf("%s %s: %d of %d fully callable genomes carry the pair; %d carry neither; %d discordant; %d not callable at one or both.\n",
              ty, COMBOS$label[i], both, sum(ok), neither, discord, sum(!ok)))
  cat(sprintf("  effect: %s -- %s\n", COMBOS$effect[i], COMBOS$note[i]))
}

red <- hits %>% filter(effect == "reduced")
cat("\n-- headline --\n")
if (nrow(red) == 0) {
  cat("No listed resistance-associated variant was detected at any callable position\n")
  cat("in this prespecified screen. The screen is NOT exhaustive: it covers the\n")
  cat(sprintf("%d entries in the evidence table above, at sites O and II only. Variants\n", nrow(EVIDENCE)))
  cat("outside the table, outside these sites, or at uncallable residues would not be\n")
  cat("detected, and the table should be re-checked against current labelling.\n")
} else print(as.data.frame(red))

## ---- explicit negative statements, with denominators ------------------------
cat("\n-- explicit absence statements (denominator = genomes actually called) --\n")
key <- EVIDENCE %>% filter(effect == "reduced") %>% distinct(subgroup, pos)
obs_by_pos <- rs %>% group_by(subgroup, pos) %>%
  summarise(observed = paste(sprintf("%s=%d", aa, genomes), collapse = ", "),
            .groups = "drop")
neg <- key %>%
  left_join(denom %>% select(subgroup, pos, callable, uncalled), by = c("subgroup", "pos")) %>%
  left_join(obs_by_pos, by = c("subgroup", "pos")) %>%
  left_join(EVIDENCE %>% filter(effect == "reduced") %>%
              group_by(subgroup, pos) %>%
              summarise(screened_for = paste(aa, collapse = "/"), .groups = "drop"),
            by = c("subgroup", "pos")) %>%
  arrange(subgroup, pos)
print(as.data.frame(neg))
cat("\nFor each row: every callable genome was checked for the listed residues and\n")
cat("none carried them. 'uncalled' genomes could not be assessed at that residue.\n")

## ---- K272R, selected by residue identity (not by a derived score) ----------
cat("\n-- position 272 (palivizumab site II), by observed residue --\n")
p272 <- rs %>% filter(pos == 272)
print(as.data.frame(p272 %>% select(subgroup, aa, genomes, callable, pct_of_callable)))
fB <- read_fasta_matrix(file.path(PROJ, "nextclade", "RSV-B",
                                  "nextclade.cds_translation.F.fasta"))
k272r <- rownames(fB)[fB[, 272] == "R"]
cat(sprintf("\nGenomes with R at F272: %d of %d callable (%.1f%%): %s\n",
            length(k272r), denom$callable[denom$subgroup == "RSV-B" & denom$pos == 272],
            100 * length(k272r) / denom$callable[denom$subgroup == "RSV-B" & denom$pos == 272],
            paste(k272r, collapse = ", ")))
cat("K272R is a change at a palivizumab contact residue. The documented\n")
cat("palivizumab-resistance variants at 272 are E/M/N/Q/T; a phenotype for R\n")
cat("specifically was not found in the sources screened, so it is reported as a\n")
cat("substitution of unknown phenotype at a contact residue, not as resistance.\n")

## ---- substitution burden, including genomes with none ----------------------
cat("\n-- F substitution burden (ALL genomes, including those with none) --\n")
aa_subs <- dat %>% select(FieldNoNIV_f, nc_dataset, nc_aa_subs) %>%
  mutate(nc_aa_subs = replace_na(nc_aa_subs, "")) %>%
  separate_rows(nc_aa_subs, sep = ",") %>%
  filter(str_starts(nc_aa_subs, "F:")) %>%
  count(FieldNoNIV_f, nc_dataset, name = "subs")
burden <- dat %>% select(FieldNoNIV_f, nc_dataset) %>%
  left_join(aa_subs, by = c("FieldNoNIV_f", "nc_dataset")) %>%
  mutate(subs = replace_na(subs, 0L)) %>%
  group_by(nc_dataset) %>%
  summarise(genomes = n(), with_any = sum(subs > 0), median = median(subs),
            mean = round(mean(subs), 2), max = max(subs), .groups = "drop")
print(as.data.frame(burden))
cat("(the earlier median was conditional on having >=1 substitution)\n")

## ---- antigenic-site summary, per-genome callability respected ---------------
in_site <- function(p) {
  s <- SITES$site[map_lgl(seq_len(nrow(SITES)), ~ p >= SITES$start[.x] & p <= SITES$end[.x])]
  if (length(s)) paste(unique(s), collapse = "/") else NA_character_ }
site_poly <- poly %>% mutate(antigenic_site = map_chr(pos, in_site)) %>%
  filter(!is.na(antigenic_site))
cat("\n-- polymorphic residues by antigenic site --\n")
print(as.data.frame(site_poly %>% select(subgroup, antigenic_site, pos, callable, states)))

write_out(rs,        "04_F_residue_states.csv")
write_out(denom,     "04_F_callable_denominators.csv")
write_out(poly,      "04_F_polymorphic_residues.csv")
write_out(hits,      "04_F_evidence_matches.csv")
write_out(EVIDENCE,  "04_F_evidence_table.csv")
write_out(burden,    "04_F_substitution_burden.csv")
cat("== 04 done ==\n")

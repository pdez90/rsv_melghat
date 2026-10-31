# 02_clade_assign.R -- Nextclade-equivalent clade assignment, implemented in R.
#
# The Nextclade CLI could not be executed in this environment, so this script
# reproduces its logic using Nextclade's own RSV-A and RSV-B reference datasets
# (nextstrain/nextclade_data): reference genome, genome annotation, and the
# clade-labelled reference tree.
#
#   1. queries are aligned to the dataset reference with MAFFT --keeplength,
#      so all coordinates stay in reference space;
#   2. every tip of the reference tree is reconstructed from the tree's root
#      sequence by applying the nucleotide mutations along its branch path;
#   3. each query is assigned the clade (and G-clade) of its nearest reference
#      tip by Hamming distance over all variable sites it actually covers.
#
# This is Nextclade's nearest-reference-node placement. It does not reproduce
# Nextclade's QC scoring or private-mutation handling.

source(file.path("scripts", "00_setup.R"))
cat("\n== 02 clade assignment ==\n")

linked <- readRDS(file.path(OUT, "01_linked.rds"))
seqdir <- file.path(OUT, "seq"); dir.create(seqdir, showWarnings = FALSE, recursive = TRUE)

DATASETS <- list(
  `RSV-A` = list(dir = file.path(REF, "rsv_a"), tag = "a"),
  `RSV-B` = list(dir = file.path(REF, "rsv_b"), tag = "b"))

## ---- write query FASTA per subgroup (sequences live in the CSV) -------------
seqs <- read_csv(data_path("NP_results_with_clade_sequencing.csv"),
                 show_col_types = FALSE, guess_max = 1000) %>%
  mutate(rsv_type = case_when(rsva == 1 ~ "RSV-A", rsvb == 1 ~ "RSV-B", TRUE ~ NA_character_),
         seq = toupper(str_replace_all(toupper(FASTA_f), "[^ACGTN\\-]", "N")))
cat(sprintf("genomes: %d (median length %d bp)\n", nrow(seqs), median(nchar(seqs$seq))))

write_fasta <- function(ids, sq, path) {
  writeLines(as.vector(rbind(paste0(">", ids), sq)), path)
}

## ---- reconstruct reference tips from the Nextclade tree ---------------------
# Each node carries branch_attrs$mutations$nuc entries like "A84G" (reference
# base A at position 84 -> G). Walking the tree root->tip and accumulating these
# gives every tip's genome as a sparse diff from the root/reference sequence.
collect_tips <- function(node, acc = list()) {
  muts <- node$branch_attrs$mutations$nuc
  if (!is.null(muts) && length(muts)) {
    m   <- unlist(muts)
    pos <- as.integer(str_sub(m, 2, nchar(m) - 1))
    alt <- str_sub(m, nchar(m), nchar(m))
    keep <- !is.na(pos) & alt %in% c("A", "C", "G", "T", "-")
    if (any(keep)) acc[as.character(pos[keep])] <- alt[keep]
  }
  if (is.null(node$children) || !length(node$children)) {
    cl <- node$node_attrs$clade_membership$value %||% NA_character_
    gc <- node$node_attrs$G_clade$value %||% NA_character_
    return(list(list(name = node$name, clade = cl, G_clade = gc, muts = acc)))
  }
  unlist(lapply(node$children, collect_tips, acc = acc), recursive = FALSE)
}

assign_subgroup <- function(type) {
  ds <- DATASETS[[type]]
  cat(sprintf("\n-- %s --\n", type))

  q <- seqs %>% filter(rsv_type == type, nchar(seq) > 1000)
  cat(sprintf("queries: %d\n", nrow(q)))
  qfa   <- file.path(seqdir, sprintf("queries_%s.fasta", ds$tag))
  alnfa <- file.path(seqdir, sprintf("aligned_%s.fasta", ds$tag))
  write_fasta(q$FieldNoNIV_f, q$seq, qfa)

  # MAFFT --keeplength keeps reference coordinates (insertions vs ref dropped)
  ref_path <- file.path(ds$dir, "reference.fasta")
  cmd <- sprintf("mafft --quiet --thread -1 --keeplength --add %s %s > %s",
                 shQuote(qfa), shQuote(ref_path), shQuote(alnfa))
  t0 <- Sys.time(); system(cmd)
  cat(sprintf("aligned in %.0fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  aln <- read_fasta_matrix(alnfa)
  refseq <- aln[1, ]; qaln <- aln[-1, , drop = FALSE]
  cat(sprintf("alignment: %d seqs x %d sites (reference coordinates)\n",
              nrow(aln), ncol(aln)))

  tj   <- fromJSON(file.path(ds$dir, "tree.json"), simplifyVector = FALSE)
  tips <- collect_tips(tj$tree)
  root <- toupper(str_split(tj$root_sequence$nuc, "")[[1]])
  cat(sprintf("reference tree: %d tips; root sequence %d bp\n", length(tips), length(root)))

  # variable sites = any site mutated in the tree, or differing in any query
  tip_pos <- unique(as.integer(unlist(lapply(tips, function(t) names(t$muts)))))
  qdiff   <- which(apply(qaln, 2, function(col) any(col %in% c("A","C","G","T") & col != refseq)))
  V <- sort(unique(c(tip_pos, qdiff)))
  V <- V[V >= 1 & V <= length(root)]
  cat(sprintf("variable sites compared: %d\n", length(V)))

  # tip base matrix over V
  tipMat <- matrix(root[V], nrow = length(tips), ncol = length(V), byrow = TRUE)
  Vidx <- setNames(seq_along(V), as.character(V))
  for (i in seq_along(tips)) {
    m <- tips[[i]]$muts
    if (length(m)) {
      j <- Vidx[names(m)]; ok <- !is.na(j)
      if (any(ok)) tipMat[i, j[ok]] <- unlist(m)[ok]
    }
  }
  tip_clade <- vapply(tips, function(t) t$clade %||% NA_character_, character(1))
  tip_gclade<- vapply(tips, function(t) t$G_clade %||% NA_character_, character(1))
  tip_name  <- vapply(tips, function(t) t$name, character(1))

  # nearest reference tip per query, over sites the query actually covers
  res <- vector("list", nrow(qaln))
  for (i in seq_len(nrow(qaln))) {
    qv  <- qaln[i, V]
    cov <- qv %in% c("A", "C", "G", "T")
    d   <- rowSums(tipMat[, cov, drop = FALSE] != matrix(qv[cov], nrow = nrow(tipMat),
                                                         ncol = sum(cov), byrow = TRUE))
    k <- which.min(d)
    res[[i]] <- tibble(
      FieldNoNIV_f   = rownames(qaln)[i],
      nextclade_clade   = tip_clade[k],
      nextclade_G_clade = tip_gclade[k],
      nearest_reference = tip_name[k],
      dist_to_nearest   = d[k],
      runner_up_clade   = tip_clade[order(d)[2]],
      runner_up_dist    = sort(d)[2],
      sites_covered     = sum(cov),
      coverage_pct      = round(100 * mean(qaln[i, ] %in% c("A","C","G","T")), 1))
  }
  out <- bind_rows(res) %>% mutate(rsv_type = type)
  saveRDS(list(aln = aln, V = V), file.path(OUT, sprintf("02_alignment_%s.rds", ds$tag)))
  out
}

calls <- bind_rows(lapply(names(DATASETS), assign_subgroup))
saveRDS(calls, file.path(OUT, "02_raw_calls.rds"))

## ---- compare with the clade column supplied with the sequencing file --------
# Lineage names are hierarchical and dotted ("A.D" -> "A.D.2" -> "A.D.2.1"), and
# the supplied column writes the bare parent with a trailing dot ("A.D."). Strip
# trailing dots before comparing, and classify any disagreement as either a
# parent/child difference along one lineage path or a genuine conflict.
norm_clade <- function(x) sub("\\.+$", "", str_trim(x))
is_ancestor <- function(a, b) !is.na(a) & !is.na(b) & a != b &
  (str_starts(b, paste0(a, ".")) | str_starts(a, paste0(b, ".")))

cmp <- linked %>% select(FieldNoNIV_f, studyid, swab_date, season, rsv_type,
                         CLADE_f, CTvalue_f, linked, severity) %>%
  left_join(calls %>% select(-rsv_type), by = "FieldNoNIV_f") %>%
  mutate(supplied = norm_clade(CLADE_f),
         called   = norm_clade(nextclade_clade),
         concordance = case_when(
           is.na(called) | is.na(supplied) ~ NA_character_,
           supplied == called              ~ "identical",
           is_ancestor(supplied, called)   ~ "parent/child on same lineage",
           TRUE                            ~ "conflicting"),
         # a finer call than the supplied column, e.g. supplied "A.D" -> "A.D.3"
         refines_supplied = !is.na(called) & !is.na(supplied) &
           str_starts(called, paste0(supplied, ".")))

cat("\n-- supplied CLADE_f vs R/Nextclade call --\n")
print(table(supplied = cmp$supplied, called = cmp$called, useNA = "ifany"))
cat("\nconcordance:\n"); print(table(cmp$concordance, useNA = "ifany"))
cat(sprintf("\nidentical calls: %d/%d (%.1f%%)\n",
            sum(cmp$concordance == "identical", na.rm = TRUE), sum(!is.na(cmp$concordance)),
            100 * mean(cmp$concordance == "identical", na.rm = TRUE)))
cat(sprintf("calls refining the supplied lineage: %d\n", sum(cmp$refines_supplied, na.rm = TRUE)))
cat(sprintf("genomes still only resolved to the parent lineage 'A.D': %d\n",
            sum(cmp$called == "A.D", na.rm = TRUE)))
if (any(cmp$concordance != "identical", na.rm = TRUE)) {
  cat("\ndisagreements:\n")
  print(as.data.frame(cmp %>% filter(concordance != "identical") %>%
    select(FieldNoNIV_f, season, supplied, called, nearest_reference,
           dist_to_nearest, runner_up_clade, runner_up_dist, coverage_pct)))
}
cat(sprintf("median distance to nearest reference tip: %d SNPs (IQR %d-%d)\n",
            median(cmp$dist_to_nearest, na.rm = TRUE),
            quantile(cmp$dist_to_nearest, .25, na.rm = TRUE),
            quantile(cmp$dist_to_nearest, .75, na.rm = TRUE)))
cat(sprintf("median genome coverage: %.1f%%\n", median(cmp$coverage_pct, na.rm = TRUE)))

write_out(cmp, "02_clade_calls.csv")
saveRDS(cmp, file.path(OUT, "02_clades.rds"))
cat("== 02 done ==\n")

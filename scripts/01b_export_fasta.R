# 01b_export_fasta.R -- write the sequenced genomes out as FASTA for Nextclade.
#
# IMPORTANT: the consensus genomes in NP_results_with_clade_sequencing.csv contain
# IUPAC ambiguity codes (R, Y, W, S, K, M, B) at mixed/uncertain sites -- 3,772
# such bases across 102 of the 241 genomes. These MUST be converted to N, not
# deleted: deleting them shifts the reading frame and makes Nextclade report
# spurious frameshifts and premature stop codons.

source(file.path("scripts", "00_setup.R"))
cat("\n== 01b export FASTA for Nextclade ==\n")

seqs <- read_csv(data_path("NP_results_with_clade_sequencing.csv"),
                 show_col_types = FALSE, guess_max = 1000) %>%
  mutate(rsv_type = case_when(rsva == 1 ~ "RSV-A", rsvb == 1 ~ "RSV-B", TRUE ~ NA_character_),
         raw = toupper(FASTA_f),
         # ambiguity codes -> N (length preserved); keep A/C/G/T/N/- as they are
         seq = str_replace_all(raw, "[^ACGTN\\-]", "N"))

cat(sprintf("genomes: %d\n", nrow(seqs)))
cat(sprintf("ambiguity codes converted to N: %d bases in %d genomes\n",
            sum(str_count(seqs$raw, "[^ACGTN\\-]")),
            sum(str_count(seqs$raw, "[^ACGTN\\-]") > 0)))
stopifnot(all(nchar(seqs$seq) == nchar(seqs$raw)))   # length must be preserved
cat(sprintf("length preserved: min %d, median %d, max %d\n",
            min(nchar(seqs$seq)), median(nchar(seqs$seq)), max(nchar(seqs$seq))))

write_fasta <- function(ids, sq, path) {
  writeLines(as.vector(rbind(paste0(">", ids), sq)), path)
  cat(sprintf("  -> %s (%d sequences)\n", basename(path), length(ids)))
}

for (ty in c("RSV-A", "RSV-B")) {
  q <- seqs %>% filter(rsv_type == ty)
  write_fasta(q$FieldNoNIV_f, q$seq,
              file.path(DATA, sprintf("nextclade_input_%s.fasta", ty)))
}
# combined file, headers carrying field number | supplied clade | swab date
write_fasta(paste(seqs$FieldNoNIV_f, seqs$CLADE_f, seqs$swab_date_new, sep = "|"),
            seqs$seq, file.path(DATA, "melghat_rsv_genomes.fasta"))

cat("== 01b done ==\n")

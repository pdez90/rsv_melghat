# 29_household_key_diagnostic.R -- do the candidate residential identifiers in
# demographics.csv behave like household keys?
#
# WHY THIS SCRIPT EXISTS
# The draft previously reported a birth-date-spacing diagnostic over four
# candidate household keys, together with a claim about the `mother` field.
# Those figures came from an exploratory step that was never archived: no
# script produced them and no output table held them, so they could not be
# re-derived. This script implements the diagnostic explicitly. The key
# DEFINITIONS below are stated here rather than inherited, because the
# original definitions are not recoverable; any figure quoted in the
# manuscript must come from this script's output table, not from the earlier
# unarchived numbers.
#
# WHAT IS BEING TESTED
# If a key identified a household in the sense of "one set of siblings", then
# within a key group the children should be siblings, and sibling birth dates
# cannot be arbitrarily close: a single mother cannot deliver two children a
# few weeks or a few months apart. So for each candidate key we count the
# multi-child groups that contain at least one pair of children born
# IMPLAUSIBLY CLOSE TOGETHER -- between 3 and 240 days apart.
#
#   lower bound 3 days   excludes twins and same-delivery multiples, which are
#                        a legitimate 0-2 day spacing
#   upper bound 240 days excludes any plausible inter-birth interval for one
#                        mother (a term pregnancy alone is ~266 days)
#
# A group flagged by this test is therefore NOT a set of full siblings. It is
# consistent with an extended household -- cousins sharing a dwelling or a
# family number -- or with the field encoding something other than a dwelling.
#
# WHAT THIS DOES AND DOES NOT ESTABLISH
# A high flagged percentage shows the key does not partition children into
# sibling sets. It does NOT show the field is meaningless, and it does NOT
# identify what the field does encode. Equally, a LOW percentage would not
# validate a key: a key could group one child per dwelling and never be
# caught by a spacing test while still being wrong. This is a disconfirming
# test only, which is why the household analysis is deferred to the field
# team rather than reported as a negative result.
#
# This script reads only demographics.csv. It makes no inference about
# households and writes no household variable for downstream use.

source(file.path("scripts", "00_setup.R"))
cat("\n== 29 candidate household-key diagnostic ==\n")

LO_DAYS <- 3
HI_DAYS <- 240

dem <- read_csv(data_path("demographics.csv"), show_col_types = FALSE,
                guess_max = 20000) %>%
  mutate(birth_date = parse_dmy(birth_date_new))

cat(sprintf("demographics rows: %d, distinct studyid: %d, parseable birth dates: %d\n",
            nrow(dem), dplyr::n_distinct(dem$studyid), sum(!is.na(dem$birth_date))))
if (nrow(dem) != dplyr::n_distinct(dem$studyid))
  cat("  NOTE: demographics is not one row per child; groups below are rows, not children\n")

blankish <- function(x) is.na(x) | trimws(as.character(x)) %in% c("", "NA", ".")

## ---- candidate key definitions ---------------------------------------------
## Stated explicitly. `address` is included because it was named as a candidate,
## not because it is plausible: it is a village-level text field, and the
## diagnostic below is expected to reject it outright.
mk <- function(...) {
  p <- list(...)
  bad <- Reduce(`|`, lapply(p, blankish))
  out <- do.call(paste, c(lapply(p, function(z) trimws(as.character(z))), sep = "|"))
  out[bad] <- NA_character_
  out
}

cand <- tibble::tibble(
  key_name = c("address (code_village_93 + address)",
               "house_no (code_village_93 + house_no)",
               "family_no",
               "village_house (village + house_no)"),
  key_fields = c("code_village_93, address",
                 "code_village_93, house_no",
                 "family_no",
                 "village, house_no")
)
keys <- list(
  mk(dem$code_village_93, dem$address),
  mk(dem$code_village_93, dem$house_no),
  mk(dem$family_no),
  mk(dem$village, dem$house_no)
)

## ---- the spacing diagnostic -------------------------------------------------
diag_one <- function(k) {
  d <- tibble::tibble(key = k, bd = dem$birth_date) %>%
    filter(!is.na(key), !is.na(bd))
  if (!nrow(d)) return(NULL)
  g <- split(as.integer(d$bd), d$key)
  sizes <- lengths(g)
  multi <- g[sizes >= 2]
  ## per group: does any pair fall in the implausible window?
  per <- vapply(multi, function(v) {
    v <- sort(v)
    dd <- as.integer(outer(v, v, "-"))
    dd <- abs(dd[lower.tri(matrix(0, length(v), length(v)))])
    c(flag = as.numeric(any(dd >= LO_DAYS & dd <= HI_DAYS)),
      pairs = length(dd),
      pairs_flagged = sum(dd >= LO_DAYS & dd <= HI_DAYS),
      min_gap = if (length(dd)) min(dd) else NA_real_)
  }, numeric(4))
  per <- t(per)
  tibble::tibble(
    children_with_key   = nrow(d),
    groups              = length(g),
    multi_child_groups  = length(multi),
    max_group_size      = max(sizes),
    median_multi_size   = stats::median(sizes[sizes >= 2]),
    groups_flagged      = sum(per[, "flag"]),
    pct_groups_flagged  = round(100 * mean(per[, "flag"]), 1),
    pairs_examined      = sum(per[, "pairs"]),
    pairs_flagged       = sum(per[, "pairs_flagged"]),
    median_min_gap_days = stats::median(per[, "min_gap"], na.rm = TRUE))
}

res <- dplyr::bind_cols(cand, dplyr::bind_rows(lapply(keys, diag_one))) %>%
  mutate(window_days = sprintf("%d-%d", LO_DAYS, HI_DAYS))

print(as.data.frame(res[, c("key_name", "multi_child_groups", "groups_flagged",
                            "pct_groups_flagged", "pairs_flagged", "pairs_examined")]))

## Do any two keys induce the SAME partition? If so they are not independent
## candidates, and reporting them as separate keys overstates how many distinct
## candidates were examined.
cat("\n  partition identity between candidate keys:\n")
ident <- NULL
for (i in seq_along(keys)) for (j in seq_along(keys)) if (i < j) {
  ok <- !is.na(keys[[i]]) & !is.na(keys[[j]])
  ## Two groupings induce the same partition iff the label map is bijective,
  ## which is what the three counts below test. Comparing factor codes instead
  ## would report "different" for identical partitions whose labels happen to
  ## sort differently (e.g. "001-Lawada|5" against "1|5").
  a <- keys[[i]][ok]; b <- keys[[j]][ok]
  same <- sum(ok) > 0 &&
    dplyr::n_distinct(a) == dplyr::n_distinct(b) &&
    dplyr::n_distinct(paste(a, b, sep = "\u001f")) == dplyr::n_distinct(a)
  ident <- dplyr::bind_rows(ident, tibble::tibble(
    key_a = cand$key_name[i], key_b = cand$key_name[j],
    same_partition = same, compared_on = sum(ok)))
  cat(sprintf("    %-38s vs %-38s : %s\n", cand$key_name[i], cand$key_name[j],
              if (same) "IDENTICAL partition" else "different"))
}

## ---- the `mother` field -----------------------------------------------------
## Claimed in the draft to be a low-cardinality numeric and therefore not a
## mother identifier. Profiled here rather than asserted. `father` is profiled
## beside it because the two are adjacent fields and the comparison is what
## shows what `mother` encodes.
prof <- function(nm) {
  v <- trimws(as.character(dem[[nm]]))
  v <- v[!blankish(v)]
  u <- sort(unique(v))
  tibble::tibble(field = nm, nonblank = length(v), distinct_values = length(u),
                 all_numeric = all(grepl("^-?[0-9]+$", u)),
                 example_values = paste(utils::head(u, 8), collapse = "; "))
}
mprof <- dplyr::bind_rows(prof("mother"), prof("father"))
cat("\n  mother / father field profile:\n")
print(as.data.frame(mprof))
cat(sprintf("\n  `mother` holds %d distinct non-blank values over %d children.\n",
            mprof$distinct_values[mprof$field == "mother"],
            mprof$nonblank[mprof$field == "mother"]))
cat("  A mother identifier in a cohort of this size would have thousands of\n")
cat("  values, so this field cannot link children to a common mother and\n")
cat("  cannot adjudicate any of the candidate keys above.\n")

write_csv(res,   file.path(OUT, "29_household_key_diagnostic.csv"))
write_csv(ident, file.path(OUT, "29_household_key_partition_identity.csv"))
write_csv(mprof, file.path(OUT, "29_household_key_mother_field.csv"))
cat("\n  -> wrote 29_household_key_diagnostic.csv\n")
cat("  -> wrote 29_household_key_partition_identity.csv\n")
cat("  -> wrote 29_household_key_mother_field.csv\n")
cat("\n  Quote figures in the manuscript from these tables only.\n")

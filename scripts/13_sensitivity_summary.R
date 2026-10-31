# 13_sensitivity_summary.R -- one comparison table and one effect-estimate plot
# from the panel written by 12_sensitivity.R. Reads only 12_* outputs.

source(file.path("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(ggplot2); library(scales)})
cat("\n== 13 sensitivity summary ==\n")

P  <- readRDS(file.path(OUT, "12_panel.rds"))
O1 <- P$O1; O2 <- P$O2; O4 <- P$O4; O5 <- P$O5

SPEC_LAB <- c(
  S0 = "S0  all links, all QC, all villages (primary)",
  S1 = "S1  exact links only",
  S2 = "S2  exact + narrow window",
  S3 = "S3  strict genome QC",
  S4 = "S4  exact links + strict QC",
  S5 = "S5  common 33 villages",
  S6 = "S6  exact + strict QC + 33 villages")

# ---- O2: convergence-aware status ------------------------------------------
# 12 sets status = "ok" whenever the fit does not error; singular/converged are
# carried separately. Re-derive an honest label here so the table cannot be read
# as if every estimate were trustworthy.
O2 <- O2 %>%
  mutate(fit_flag = case_when(
    status != "ok"            ~ status,
    !converged &  singular    ~ "did not converge; singular",
    !converged                ~ "did not converge",
    singular                  ~ "singular (variance component at 0)",
    TRUE                      ~ "converged"))

# ---- comparison table -------------------------------------------------------
cov_by_spec <- O1 %>%
  group_by(spec) %>%
  summarise(cohort_rsv_pos = sum(cohort_rsv_pos), linked = sum(linked),
            pct = round(100 * linked / cohort_rsv_pos, 1), .groups = "drop")

spat_by_spec <- O4 %>%
  group_by(spec) %>%
  summarise(strata_total    = n(),
            strata_analysed = sum(analysed),
            # does the genome-level bootstrap interval exclude zero?
            ci_excl_zero    = sum(analysed & !is.na(boot_lo) &
                                    (boot_lo > 0 | boot_hi < 0)),
            perm_p_lt_05    = sum(analysed & !is.na(perm_p) & perm_p < 0.05),
            .groups = "drop")

sev_by_spec <- tibble(spec = names(SPEC_LAB)) %>%
  mutate(vill = if_else(spec %in% c("S5", "S6"), "V33", "all")) %>%
  left_join(O2 %>% select(vill, OR, lo, hi, p, fit_flag), by = "vill")

esc_by_spec <- O5 %>%
  filter(population == first(population)) %>%          # primary population
  summarise(variants = n_distinct(paste(subgroup, variant)),
            detections = sum(detected)) %>%
  mutate(spec = "all specs")

TAB <- tibble(spec = names(SPEC_LAB), label = unname(SPEC_LAB)) %>%
  left_join(cov_by_spec, by = "spec") %>%
  left_join(spat_by_spec, by = "spec") %>%
  left_join(sev_by_spec %>% select(spec, OR, lo, hi, p, fit_flag), by = "spec") %>%
  mutate(coverage       = sprintf("%d/%d (%.1f%%)", linked, cohort_rsv_pos, pct),
         severity_OR    = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi),
         spatial_strata = sprintf("%d of %d", strata_analysed, strata_total),
         spatial_signal = sprintf("%d bootstrap / %d permutation",
                                  ci_excl_zero, perm_p_lt_05)) %>%
  select(spec, label, coverage, severity_OR, severity_p = p, severity_fit = fit_flag,
         spatial_strata, spatial_signal)

cat("\n-- Comparison table --\n")
print(as.data.frame(TAB), right = FALSE)
write_csv(TAB, file.path(OUT, "13_comparison_table.csv"))

cat("\n-- O5 escape-variant screen, primary population --\n")
cat(sprintf("%d listed variants screened, %d detected\n",
            esc_by_spec$variants, esc_by_spec$detections))

# ---- effect-estimate plot ---------------------------------------------------
# Panel A: severity OR (log scale, null = 1). Panel B: same- minus
# different-village SNP difference (null = 0), one row per spec x stratum.
PAL <- c("RSV-A" = "#2a78d6", "RSV-B" = "#eb6834")
SURF <- "#fcfcfb"
thm <- theme_minimal(base_size = 10) +
  theme(panel.background = element_rect(fill = SURF, colour = NA),
        plot.background  = element_rect(fill = SURF, colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "#e5e5e2", linewidth = .3),
        axis.title.y = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(colour = "#5c5c58", size = 9),
        legend.position = "top", legend.title = element_blank())

pA <- TAB %>%
  left_join(sev_by_spec %>% select(spec, OR, lo, hi), by = "spec") %>%
  mutate(spec = factor(spec, levels = rev(names(SPEC_LAB))),
         unstable = severity_fit != "converged") %>%
  ggplot(aes(x = OR, y = spec)) +
  geom_vline(xintercept = 1, colour = "#9a9a95", linewidth = .4) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = .9, colour = "#2a78d6") +
  geom_point(aes(shape = unstable), size = 2.6, colour = "#2a78d6",
             fill = SURF, stroke = .9) +
  scale_shape_manual(values = c(`FALSE` = 19, `TRUE` = 21), guide = "none") +
  scale_x_log10(breaks = c(.5, .75, 1, 1.5, 2)) +
  labs(title = "A. Very-severe LRI given RSV positivity",
       subtitle = "Odds ratio, 95% Wald CI. Open point = model singular or non-convergent.",
       x = "Odds ratio (log scale)") + thm

pB <- O4 %>%
  filter(analysed) %>%
  mutate(stratum = paste(rsv_type, season),
         spec = factor(spec, levels = rev(names(SPEC_LAB)))) %>%
  ggplot(aes(x = obs_raw, y = spec, colour = rsv_type)) +
  geom_vline(xintercept = 0, colour = "#9a9a95", linewidth = .4) +
  geom_linerange(aes(xmin = boot_lo, xmax = boot_hi), linewidth = .9) +
  geom_point(size = 2.4) +
  facet_wrap(~stratum, ncol = 2, scales = "free_x") +
  scale_colour_manual(values = PAL) +
  labs(title = "B. Within- minus between-village genetic distance",
       subtitle = paste("Median SNP difference, 95% genome-level bootstrap CI",
                        "(999 reps, self-pairs excluded).",
                        "\nNegative = genomes from the same village are more similar."),
       x = "Median SNP difference (same village - different village)") + thm

ggsave(file.path(OUT, "13_effect_estimates_A_severity.png"), pA,
       width = 7.2, height = 3.4, dpi = 200, bg = SURF)
ggsave(file.path(OUT, "13_effect_estimates_B_spatial.png"), pB,
       width = 8.4, height = 6.2, dpi = 200, bg = SURF)
cat("  -> wrote 13_comparison_table.csv, 13_effect_estimates_A_severity.png,",
    "13_effect_estimates_B_spatial.png\n")
cat("== 13 done ==\n")

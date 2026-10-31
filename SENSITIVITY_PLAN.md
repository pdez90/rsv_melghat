# Sensitivity analysis plan — Melghat RSV genomics

**Version 3, 30 October 2026.** Supersedes versions 1 and 2 of the same date. Frozen
before implementation; the frozen copy and input hashes are recorded in
`output/12_frozen_plan.txt`.

**Status of knowledge.** The primary analysis (S0) has already been run and its results
examined; they are in `output/` and summarised in the project notes. The sensitivity
specifications S1–S6 defined here have **not** been run. This is a pre-specification of
the sensitivity panel only, not of the primary analysis. Any rule changed after seeing
S1–S6 results must be recorded here with its reason, and the original rule kept visible.

To be implemented in `scripts/12_sensitivity.R`.

---

## 1. What this panel is for

Three choices could drive the primary findings and have not been tested:

1. **Linkage.** 13 of 203 linked genomes matched a *prior* episode rather than an exact
   one. Wrong attributions would propagate to every episode-level result.
2. **Sequence quality and callability.** Genomes are currently used regardless of QC,
   and callability varies by residue and by gene.
3. **Sampling frame.** 2016/17–2019/20 cover 93 villages; 2021/22 covers a 33-village
   subset, so season comparisons are also village-set comparisons.

---

## 2. Definitions fixed in advance

### 2.1 Linkage sets

Linkage is read from `link_status` in `output/02_dat.rds`, written by
`01_link_sequencing.R`. **It is not re-derived by comparing date columns.** The
sequencing file carries its own `ai_date_new` equal to the swab date, so a test like
`swab_date == ai_date_new` compares a column with itself and would silently retain
near-linked specimens.

| Set | Rule | Expected n |
|---|---|---|
| **L-ALL** (primary) | `link_status %in% c("exact", "nearest_prior_episode_same_RSV_type")` | 203 |
| **L-EXACT** | `link_status == "exact"` | 190 |
| **L-NARROW** | `link_status == "exact"`, or nearest-prior with `days_episode_to_swab <= 14` | to be reported |

`L-NARROW` uses the recorded gap, not recomputed dates.

### 2.2 Residue and sequence state — four categories, tracked separately

A translated position in a genome is exactly one of:

| State | Meaning |
|---|---|
| **canonical** | one of the 20 amino-acid letters |
| **unknown** | `X` — sequence not resolved at this position |
| **deletion** | gap forming a length difference with intact flanking sequence |
| **stop** | `*`, including the terminal stop and any premature stop |

This separation is required. A rule counting a deletion as "uncallable" would exclude
precisely the 70 RSV-A genomes carrying the 24-residue G length difference — removing
the structural state the G analysis exists to study. **A deletion never makes a genome
ineligible for an analysis of that deletion.**

### 2.3 Quality sets

| Set | Rule |
|---|---|
| **Q-ALL** (primary) | all genomes |
| **Q-STRICT** | `qc.overallStatus == "good"` AND `coverage >= 0.95` AND no frameshift anywhere in the genome |

`Q-STRICT` is genome-level and whole-genome because it is applied to whole-genome
distance analyses, where a frameshift anywhere affects the consensus. For single-gene
analyses the gene-level rule in §2.4 applies **instead of**, not in addition to, the
frameshift clause.

### 2.4 Eligibility, per analysis

- **Per-residue screens (O5).** A genome contributes to a residue if that residue is
  *canonical* in it. The denominator is the number of such genomes, reported for every
  residue. Unknown, deletion and stop states are each counted and reported separately;
  none is treated as the reference state.
- **Per-gene protein analyses (F, G).** A genome is eligible for gene *g* if it has no
  frameshift in *g* and fewer than 10 **unknown** residues in the analysed region.
  Deletions and stops do not count toward that threshold.
- **Whole-genome distance (O4).** The quality set named by the specification. S0–S2 and
  S5 use Q-ALL; S3, S4 and S6 use Q-STRICT. Q-STRICT is not imposed on O4 regardless of
  specification — doing so would make S0–S2 and S5 silently identical to S3.

### 2.5 Genetic distance

Two measures, computed and **reported separately, never mixed**:

| Measure | Definition |
|---|---|
| **raw** | count of mismatching sites, pairwise deletion (current primary) |
| **normalised** | mismatches ÷ shared callable sites, as **substitutions per 1,000 callable sites** |

- **Shared callable sites** for a pair = positions where both genomes have an
  unambiguous base (A/C/G/T). Gaps and `N` are excluded from numerator and denominator.
- **Pair inclusion threshold: 13,000 shared callable sites**, about 85% of the ~15.2 kb
  genome. Rationale: the two lowest-coverage genomes (75.8%, 85.8%) fall below it
  against most partners, while genomes at the ~99.5% median retain essentially all
  pairs. This is a declared choice, not a derived optimum; results are reported at
  13,000 and, as a check, at 12,000 and 14,000.
- **Pairs below threshold are dropped, not imputed.** The number dropped is reported per
  season and per specification.
- **In the permutation test**, dropped pairs are excluded from the observed statistic and
  from every permuted replicate, so the null is generated on the same pair set as the
  observed value. Village labels are permuted across genomes; the pair-exclusion mask is
  held fixed.

### 2.6 Common-village set (V-33)

- **Membership:** villages with `subset_33_villages == 1` in `ari.csv`, matched to
  `code_village`.
- **Applied to numerator and denominator alike.** Every rate under V-33 uses only
  episodes, deaths, genomes and covariates from those villages, in every season.
  Restricting only genomes would compare a 33-village numerator with a 93-village
  denominator and is prohibited.
- **Stated in advance:** V-33 equalises the village set across seasons. It does **not**
  remove calendar confounding — surveillance intensity, testing rate and population
  immunity still differ by year — and must not be described as doing so.

---

## 3. Analysis population, defined per outcome

The central correction in version 2. The specifications must not silently change the
scientific population or the question. Linkage and genome-QC restrictions apply **only**
to outcomes whose unit is a sequenced genome.

| Outcome | Unit | Population | Restrictions that apply |
|---|---|---|---|
| **O1 coverage** | episode | cohort RSV-positive episodes (denominator); linked episodes (numerator) | linkage, village. Reported twice — linkage coverage, and coverage by QC-passing genomes — as distinct quantities |
| **O2 severity** | episode | **all tested cohort episodes** (`nswab_episode == 1`), sequenced or not | **village only.** Linkage and genome QC are irrelevant and must not be applied |
| **O3 lineage severity** | genome-episode | linked genomes with a lineage call | linkage, genome QC, village |
| **O4 village clustering** | genome pair | linked genomes with a village centroid | linkage, **the quality set of the specification** (Q-ALL or Q-STRICT), village |
| **O5 F variant screen** | genome | **all 241 sequenced genomes** | **none by default.** Linked and exact subsets reported as labelled sensitivities, never as the primary population |

O2's primary population is 4,081 tested episodes; applying L-ALL would wrongly reduce it
to 203 and change the question from "is RSV-positive illness more severe" to "are
sequenced RSV infections more severe". O5's primary population is 241; applying L-ALL
would silently drop the 38 specimens that carry no cohort episode but are perfectly good
sequence.

---

## 4. Specifications

| # | Linkage | Quality | Villages | Purpose |
|---|---|---|---|---|
| S0 | L-ALL | Q-ALL | all | primary, already run |
| S1 | L-EXACT | Q-ALL | all | does inferred linkage drive anything? |
| S2 | L-NARROW | Q-ALL | all | is the 45-day window too generous? |
| S3 | L-ALL | Q-STRICT | all | does sequence quality drive anything? |
| S4 | L-EXACT | Q-STRICT | all | both restrictions |
| S5 | L-ALL | Q-ALL | V-33 | does the changing village set drive season differences? |
| S6 | L-EXACT | Q-STRICT | V-33 | most restrictive |

Each outcome is computed only under the restrictions §3 says apply to it. For O2 that
means two cells, not seven: all-village and V-33.

---

## 5. Reporting rules

The version-1 rule — "robust if the estimate stays inside S0's confidence interval" —
is **removed**. It is unsound: a wide S0 interval would make large changes look robust,
and losing significance after a sample-size reduction does not by itself demonstrate
specification dependence.

Each outcome × specification cell reports together:

1. **effect estimate** and **confidence interval**;
2. **direction**, where direction is meaningful;
3. **sample loss** — how many genomes, episodes or pairs the restriction removed, and why;
4. **substantive change** — the difference from S0 on the outcome's own scale, as a
   magnitude, not a significance verdict.

Interpretation distinguishes:

- **attenuation** — point estimate moves materially toward the null, precision comparable;
- **loss of precision** — interval widens, point estimate roughly unchanged; a
  consequence of reduced n, not evidence about the finding;
- **sign change** — point estimate crosses the null.

For zero-detection variant screens, "direction" is not meaningful: report callable counts
and an uncertainty bound instead.

No specification is promoted to primary after the fact. S0 stays primary.

### Uncertainty, by outcome

**O4 (spatial).** A permutation p-value is not a confidence interval. The interval is a
**genome-level bootstrap**: genomes are resampled with replacement within season, the
same-minus-different-village median difference recomputed, 999 replicates, percentile
interval. Pairs formed between a genome and its own duplicate are excluded from each
replicate, since they are not observations. Pairs are **never** resampled as if
independent — the unit of independence is the genome, not the pair. The bootstrap
interval and the permutation p-value are reported side by side and described as
answering different questions.

**O5 (variant screen), zero detections.** The upper bound is an **independent-binomial
approximation** and is labelled as such wherever it appears: these genomes are related
by descent and are not independent draws from a population, so the true bound is wider
than stated. Use the **exact one-sided 97.5% binomial bound** (`binom.test(0, n)`), not
the 3/n rule of thumb, because several denominators here are small. If no genome is
callable at a residue, report **"not assessed"** — never a bound.

### Sparse data, by outcome

The version-1 blanket rule ("fewer than 5 events = not estimable") is **removed**; it
would have suppressed the F screen, where zero detections is the finding.

| Outcome | Rule |
|---|---|
| **O5** | zero detections is an estimable descriptive result. Report the count, the callable denominator, and the exact one-sided binomial upper bound, labelled an independent-binomial approximation. "Not assessed" only when no genome is callable. Never "not estimable" |
| **O3** | an adjusted model is fitted only with ≥10 genomes, ≥5 events **and ≥5 non-events** per group. Separation, convergence and singularity are checked and recorded for every fit; a model that fails any check is reported as not estimable with the reason, not as a number. Below the threshold, report descriptive proportions with exact binomial intervals rather than nothing |
| **O4** | a season is analysed only with ≥8 genomes, ≥5 villages, ≥5 same-village pairs and ≥20 different-village pairs after pair exclusions. All four counts are reported whether or not the season qualifies |
| **O1, O2** | report all cells; no suppression |

---

## 6. Declared limits

- Read-level data are unavailable, so no specification can validate a variant call or
  resolve the G repeat-length question.
- The O4 permutation does not condition on collection date or lineage. Time- and
  lineage-conditioned permutation is separate work, **not** part of this panel.
- Nothing here addresses severity-dependent testing selection. Restricting to tested
  episodes avoids miscoding untested children as negative; it does not make the tested
  set representative, and no specification may be described as correcting selection bias.
- Ct is **not** used as a covariate in this panel. Sequencing was attempted on low-Ct
  specimens (17–28), so Ct is a selection variable for which genomes exist at all, and
  adjusting for it would answer a different question from the one asked. Whether Ct
  adjustment would induce collider bias depends on the causal structure assumed and is
  not asserted here; it is out of scope.
- Only 2 sequenced episodes are deaths; no mortality outcome is estimable in any cell.
- The F evidence table is curated, not exhaustive, and was assembled after inspecting
  these genomes. O5 tests robustness of counting, not of variant selection.

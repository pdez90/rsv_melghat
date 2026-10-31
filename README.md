# Melghat RSV genomics — sequencing analysis

RSV whole-genome sequencing (241 genomes, NIV Pune) linked to the MAHAN Melghat
under-2 LRTI cohort. Run everything from this folder.

The pipeline is mostly R (29 scripts: cohort linkage, all statistics, all
figures) with 7 Python scripts and 3 shell drivers. Python is used only where it
is working on sequence alignments character by character, or where the job must
not depend on R being installed: alignment and coordinate handling
(`15_contextual_prep.py`, `17_repeat_region.py`, `20_gapfrac.py`,
`21_prepare_context.py`), comparison of residue states against the contextual
panel (`24_f_context_comparison.py`, `25_g_context_comparison.py`), and the
reproducibility freeze (`99_freeze.py`). `scripts/run_all.R` drives the R
pipeline; the shell scripts in `contextual/` drive MAFFT and FastTree.

## What you need to run it

Check everything at once — this changes nothing, it only reports:

    cd ~/Downloads/RSV_India/Sequencing
    Rscript scripts/check_requirements.R

**R >= 4.1**, and these packages:

    install.packages(c("readr","dplyr","tidyr","stringr","purrr","ggplot2",
                       "jsonlite","ape","haven","lme4","broom.mixed"))

Run that line **inside R**, or from the shell as
`Rscript -e 'install.packages(...)'` — it is not a shell command.

`dplyr` must be >= 1.1.0 (the linkage uses `relationship = "many-to-many"`).

**Nothing else needs installing.** Nextclade was run in the browser at
clades.nextstrain.org, which analyses sequences locally and needs no software;
its output is already in `nextclade/`. MAFFT is optional and only used by
`02b_clade_assign_crosscheck.R`, which `run_all.R` does not call.

Then:

    Rscript scripts/run_all.R

## Layout

    data/        sequencing inputs (NP results CSV, Nextclade input FASTAs)
    nextclade/   Nextclade Web output, RSV-A/ and RSV-B/ (+ run1_superseded/)
    scripts/     numbered R pipeline
    output/      results written by the scripts
    archive_python/   superseded first-pass Python

Cohort data is NOT duplicated here. `scripts/00_setup.R::data_path()` resolves
input files across `data/`, `../Final/Data/` (ari, lrtideaths, demographics),
`../GISproject/GISProject/capstone project Mattie/` (household GPS) and
`../` (village_characteristics.dta). If a file is missing it fails with a message
naming every folder it searched.

## Pipeline

    00_verify_data.R      data identity, counts vs manuscripts, known-bad variables
    01_link_sequencing.R  link 241 swabs to LRTI episodes (studyid + date)
    01b_export_fasta.R    write Nextclade input FASTAs
    02_nextclade_results.R  ingest Nextclade output, merge onto cohort
    03_lineage_turnover.R   subgroup/lineage composition by season, seasonality
    04_fgene_antigenic.R    F-protein variation at antibody target sites
    05_phylogeography.R     genetic distance vs village geography
    06_severity_reinfection.R  severity models, reinfections
    07_data_gaps.R          per-sample disposition, exclusions, quality losses
    08_robustness.R         permutation tests, cluster checks
    09_ggene.R              G-protein variation, length polymorphism, glycosylation
    10_validation.R         translation boundaries, normalised G-vs-F measures

    02b_clade_assign_crosscheck.R  optional; independent R clade assignment (needs MAFFT)

## Traps in this dataset — read before analysing

**Testing denominator.** `rsva_episode` / `rsvb_episode` in ari.csv are **0, not
NA, for episodes that were never swabbed**. Only 4,081 of 9,236 episodes were
tested, and swabbing was heavily selected on severity (40.1% very severe among
swabbed vs 4.6% unswabbed). Restrict every RSV-vs-non-RSV comparison to swabbed
episodes, or the RSV effect is badly inflated.

**Two severity definitions.** `lri_cat` and the WHO flags (`v_severe_lri`,
`severe_lri_who`) disagree on 276 and 377 episodes. The published child-level
counts are only reproducible by mixing them. This pipeline uses `lri_cat`.

**Ambiguity codes.** The consensus genomes contain IUPAC codes (R, Y, W, S, K, M, B)
— 3,772 bases across 102 genomes. `01b_export_fasta.R` converts them to N. They
must never be deleted: that shifts the reading frame and produces spurious
frameshifts.

**Translation boundaries.** Nextclade CDS translations include the terminal stop,
so "575 codons" is 574 amino acids plus stop. Half the RSV-B genomes have lost
the G stop codon (TAA->CAA) and encode a longer protein than the reference window
shows.

**Linkage key.** Link on `studyid`. The migration uid columns were never updated
correctly (Rowena Crow, 3 Nov 2025) and are not usable.

## Linkage summary

203 of 241 sequenced swabs (84.2%) link to an LRTI episode: 190 exact
studyid + episode-date matches, 13 to the child's most recent earlier episode of
the same RSV subgroup within 45 days (flagged in `link_status`). Of the 38
unlinked, 28 are 2021/22 swabs from children already aged 2 or over, 6 have no
matching episode, and 4 study IDs are absent from demographics.csv. The cohort
has 644 RSV-positive episodes, so sequencing covers about 31.5%.
See `output/Melghat_RSV_genotypic_data_gaps.docx` for the full accounting.

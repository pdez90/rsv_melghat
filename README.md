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

Python 3 is required for stages D, E and F of the pipeline (alignment
coordinate handling, the G-repeat helper, and the F/G contextual comparisons)
and for the freeze. The Python scripts use the standard library only -- no
third-party packages, no Biopython -- so any python3 will do.
Nextclade itself needs no installation: it was run in the browser at
clades.nextstrain.org, which analyses sequences locally.

**This repository contains code only.** No participant data, no genome
sequences, no Nextclade output and no result tables are tracked here, because
the cohort data are identifiable and the sequences are not ours to redistribute.
Everything the pipeline reads must be supplied separately -- see "Inputs you must
obtain" below. A clone of this repository therefore will NOT run end to end as
cloned.

External binaries, needed only to rebuild the contextual alignments and trees:

    mafft       contextual/run_add.sh     builds combined_RSV-*.fasta
    FastTree    contextual/run_tree.sh    builds tree_RSV-*.nwk
                contextual/run_masked.sh  builds tree_masked_RSV-*.nwk

Then:

    Rscript scripts/run_all.R

The driver runs in stages and reports which ones ran. If a stage's inputs are
absent it is SKIPPED, named in the closing summary, and the run exits non-zero,
so a partial run cannot be mistaken for a complete one. Pass `--allow-skips` if
a partial run is what you intend, or `--core` to skip the slow resampling panel.
Before any tree-dependent analysis, `28_dependency_check.R` runs two checks of
different strength, and stops the run on a failure of either.

**Taxon-set consistency.** Each tree's tip labels must equal the sequence names in
the alignment currently on disk. Failing this disproves correspondence; passing it
does **not** establish that the tree was built from these inputs, because two
alignments can share sequence names while differing in nucleotides, columns or
mask. File existence alone is weaker still.

**Build provenance.** `contextual/run_tree.sh` and `run_masked.sh` write
`<tree>.provenance` beside each tree, recording the SHA-256 of the alignment (and,
for masked trees, the mask) actually used. Where that file exists the digests are
re-computed from the inputs now on disk and compared, which does establish
correspondence. Trees built before provenance recording was added report
`UNVERIFIED` rather than passing; rebuild with the two shell scripts to close it.
Rebuilding the trees needs **FastTree only** — `run_tree.sh` and `run_masked.sh`
work from the combined alignments already on disk and do not re-align, so MAFFT is
not required and neither the alignments nor the mask definitions change. The tree
files do change, so every tree-dependent result must be re-checked afterwards.

Timestamps are advisory only — here the decompressed contextual input is seconds
newer than the alignment built from it, so an mtime rule would raise a false
alarm. Neither check verifies branch lengths or model settings.

Freeze afterwards with `python3 scripts/99_freeze.py`, and verify with
`python3 scripts/99_freeze.py --verify`. Verification reports two independent
guarantees and fails on either: **live files** re-hashed against the manifest
(changed, missing or unverifiable entries all fail), and **snapshot integrity** —
every archived copy re-hashed, gzipped copies after decompression, so a corrupt or
absent archive cannot pass just because the live file still agrees.

## Inputs you must obtain

None of these are in this repository.

| Input | Where the pipeline looks | How to obtain |
|---|---|---|
| `ari.csv`, `lrtideaths.csv`, `demographics.csv` | `../Final/Data/` | MAHAN study team; controlled access, identifiable cohort records |
| `NP_results_with_clade_sequencing.csv` | `data/` | MAHAN study team / NIV Pune |
| Household GPS | `../GISproject/GISProject/capstone project Mattie/` | MAHAN study team; household coordinates |
| `village_characteristics.dta`, `93 village list.csv` | `../` | MAHAN study team |
| Melghat genome FASTAs | `data/` | NIV Pune; deposition pending |
| Nextclade output (`nextclade.tsv`, `.csv`, `.aligned.fasta`) | `nextclade/RSV-A/`, `nextclade/RSV-B/` | regenerate: upload the FASTAs from `01b_export_fasta.R` to clades.nextstrain.org, datasets `EPI_ISL_412866` (RSV-A) and `EPI_ISL_1653999` (RSV-B) |
| INFORM contextual panel | `contextual/inform/` | Langedijk et al. 2024, Nat Commun 15:3083 |
| Nextclade reference files | `contextual/reference/` | the same Nextclade datasets |

Contact the corresponding author for the controlled items. Checksums for every
input the analysis actually used are recorded in the freeze manifest, so a
recipient can confirm they hold the same files.

## Layout

    data/        sequencing inputs (NP results CSV, Nextclade input FASTAs) -- NOT tracked
    nextclade/   Nextclade Web output, RSV-A/ and RSV-B/                     -- NOT tracked
    contextual/  alignments, trees, INFORM panel, reference  -- only the .sh drivers tracked
    scripts/     numbered R and Python pipeline                              -- tracked
    output/      results written by the scripts                             -- NOT tracked
    freeze/      dated snapshots + SHA-256 manifest                         -- NOT tracked
    manuscript/  Methods and Results                                        -- NOT tracked

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
has 644 RSV-positive episodes. Coverage is reported on the declared definition --
an episode is RSV-positive if the cohort flags say so -- which gives
**201/644 (31.2%)**: two linked genomes are cohort-negative against a positive
laboratory result and are excluded from the coverage numerator while being
retained in the genomic analyses. See `scripts/27_cohort_status_audit.R`. The
wider all-linked numerator, 203/644 (31.5%), is reported alongside in
`output/03_cohort_rsv_by_season.csv` so the two can be read against each other.
See `output/Melghat_RSV_genotypic_data_gaps.docx` for the full accounting.

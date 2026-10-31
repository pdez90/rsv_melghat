# Adding contemporaneous (2021–22) contextual genomes

The INFORM panel ends April 2020 with zero later sequences. Until 2021–22
genomes are added, the 2021/22 Melghat clade cannot be interpreted
geographically and introduction counts cannot be estimated.

NCBI, Zenodo and GISAID are **not reachable** from either the cloud container or
the workspace VM, so this download has to be done in a normal browser.

## 1. Download from NCBI Virus

Go to <https://www.ncbi.nlm.nih.gov/labs/virus/vssi/>

Run the query **twice**, once per subgroup:

| | taxid | |
|---|---|---|
| RSV-A | 208893 | `Human orthopneumovirus` / RSV-A |
| RSV-B | 208895 | RSV-B |

Filters:

- **Nucleotide Completeness**: `complete`
- **Sequence Length**: `14000` to `16000`
- **Collection Date**: `2020-01-01` to `2023-12-31`
  (overlapping 2020 deliberately — it gives a join to the existing panel)
- **Geographic Region**: run once unrestricted, then once restricted to
  `India` and neighbouring South Asian countries, so the Indian sequences are
  not swamped. Keep both downloads.

Download → **Nucleotide** → FASTA. In *Customize fields*, set the defline to:

    Accession | Geo Location | Collection Date

Save as `contextual/incoming/ncbi_RSV-A_2020_2023.fasta` and the RSV-B
equivalent.

## 2. Normalise the headers

The panel's headers are `accession_Country_Date_Subgroup`. Convert with:

    python3 scripts/21_prepare_context.py

It rewrites headers, drops sequences shorter than 14,000 nt or with >5% N,
removes accessions already in the INFORM panel, reports how many sequences
survive per year and country, and writes
`contextual/inform/RSV.{A,B}.aln.wGB.noRC.knownDate.plus2021.fasta`.

Note these are **unaligned** GenBank records, so the merged panel has to be
realigned, not concatenated. The script writes a combined unaligned FASTA and
`run_add.sh` picks it up when `INFORM_PLUS=1` is set:

    INFORM_PLUS=1 ./run_add.sh      # needs mafft; ~10-20 min with the extra sequences
    ./run_tree.sh
    cd .. && python3 scripts/15_contextual_prep.py
    Rscript scripts/16_tree_annotation.R

## 3. What to check once rebuilt

- Do the large Melghat clades stay Melghat-only, or do 2021–22 global sequences
  fall inside them? If sequences from elsewhere interleave, the single-clade
  reading is wrong.
- Does placement support for the 2021/22 RSV-A clade rise above its current
  0.900?
- Do any contemporaneous sequences share the duplication-absent state? That
  would be the strongest available evidence that the G finding is biological
  rather than an assembly artefact — short of the reads themselves.

Until that is done: **no introduction counts, no geographic-origin claims.**

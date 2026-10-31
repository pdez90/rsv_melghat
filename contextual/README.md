# Melghat contextual-alignment package

Unpack into `RSV_India/Sequencing/`. Everything below is runnable from the
project directory; no absolute paths.

    cd RSV_India/Sequencing/contextual
    gunzip -k *.gz inform/*.gz          # alignments ship compressed
    sha256sum -c SHA256SUMS.txt         # verify before use
    ./run_add.sh                        # rebuild combined alignments (~4 min)
    ./run_tree.sh                       # exploratory FastTree trees (slow)
    cd .. && python3 scripts/15_contextual_prep.py

## Contents

    contextual/run_add.sh            builds combined alignments; writes SHA256SUMS.txt
    contextual/run_tree.sh           FastTree; SH-like local support, NOT bootstrap
    contextual/combined_RSV-*.fasta.gz   INFORM + Melghat + reference, mafft --add
    contextual/query_RSV-*.fasta.gz      what was added (reference + Melghat)
    contextual/inform/*.gz           Langedijk et al. 2024 contextual alignments
    contextual/reference/            Nextclade reference + annotation, both subgroups
    contextual/mafft_RSV-*.log       alignment logs
    contextual/SHA256SUMS.txt        checksums over every input and output
    scripts/15_contextual_prep.py    coordinate map, validation, insertion audit
    output/15_coordmap_*.csv         ref_pos -> combined_col, gene, ref_base
    output/15_insertions_*.csv       per-genome Nextclade insertion accounting
    output/15_provenance.txt         SHA-256 over all of the above

## Two things to keep in mind

FastTree node values are **SH-like local support**, not whole-tree bootstrap.
Do not report them as bootstrap. A publication tree needs IQ-TREE or RAxML.

Residue numbering for F and G goes through the per-gene CDS translations
(`nextclade.cds_translation.{F,G}.fasta`), never a genome offset: RSV-A G has
31 distinct coordinate sets across 130 genomes because of length polymorphism.

#!/bin/bash
# Run BLAST against local PDB database
blastp -query data/query.fasta \
       -db /path/to/pdb_seqres.txt \
       -outfmt "6 sseqid pident length evalue bitscore" \
       -max_target_seqs 20 \
       -out results/blast_hits.tsv

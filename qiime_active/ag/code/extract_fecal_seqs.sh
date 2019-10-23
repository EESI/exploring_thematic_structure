#!/bin/bash

DATA_DIR=~/AG/filtered_fasta_combined
SAMPIDS=~/AG/samps_fecal.txt
FASTA=$DATA_DIR/combined_seqs.fna
OUT=$DATA_DIR/combined_seqs_fecal.fna

filter_fasta.py -f $FASTA -o $OUT --sample_id_fp $SAMPIDS

exit 0

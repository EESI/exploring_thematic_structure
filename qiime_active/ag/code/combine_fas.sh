#!/bin/bash

DATA_DIR=~/AG
MAP=$DATA_DIR/map.txt
FASTAS=$DATA_DIR/filtered_fasta
OUT=$DATA_DIR/filtered_fasta_combined

add_qiime_labels.py -m $MAP -i $FASTAS -c InputFileName -o $OUT

exit 0

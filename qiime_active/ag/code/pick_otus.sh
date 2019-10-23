#!/bin/bash

DATA_DIR=~/AG
OUT_DIR=~/Dropbox/stm_microbiome/qiime_active/ag
CODE_DIR=$OUT_DIR/code

REF_SEQS=~/gg_13_5_otus/rep_set/97_otus.fasta
REF_TAX=~/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt

OTUS=$OUT_DIR/picked_otus
SEQS=$DATA_DIR/seqs_filt_rmbloom.fna

PARAMS=$CODE_DIR/sortmerna_params_otus.txt
printf "pick_otus:otu_picking_method sortmerna\npick_otus:threads 50\npick_otus:similarity 0.97" > $PARAMS


pick_closed_reference_otus.py -i $SEQS                  \
                              -o $OTUS                  \
                              -r $REF_SEQS              \
                              -t $REF_TAX               \
                              -p $PARAMS                \
                              -f

exit 0

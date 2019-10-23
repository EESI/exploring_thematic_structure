#!/bin/bash

DATA_DIR=~/AG
SEQS_DIR=$DATA_DIR/filtered_fasta_combined
CODE_DIR=~/Dropbox/stm_microbiome/qiime_active/ag/code

seqs_filt=$SEQS_DIR/combined_seqs.fna
seqs_fecal=$SEQS_DIR/combined_seqs_fecal.fna
seqs_bloom=$DATA_DIR/seqs_bloom.fna
observed_bloom=$DATA_DIR/observed_bloom
seqs_filt_bloom=$DATA_DIR/seqs_filt_rmbloom.fna

params=$CODE_DIR/sortmerna_params_bloom.txt
printf "pick_otus:otu_picking_method sortmerna\npick_otus:threads 50" > $params

echo Picking bloom OTUs from $seqs_fecal

pick_closed_reference_otus.py -i $seqs_fecal            \
                              -o $observed_bloom        \
                              -r $seqs_bloom            \
                              -p $params                \
                              -f

echo Filtering blooms from $seqs_filt

filter_fasta.py -f $seqs_filt         \
                -o $seqs_filt_bloom   \
                -m $observed_bloom/sortmerna_picked_otus/combined_seqs_fecal_otus.txt \
                -n      

exit 0

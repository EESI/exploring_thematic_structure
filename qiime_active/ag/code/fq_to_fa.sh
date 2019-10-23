#!/bin/bash

fqs_dir=~/AG/filtered
fas_dir=~/AG/filtered_fasta

mkdir $fas_dir

fqs=($(ls $fqs_dir))

for fq in ${fqs[@]}
do
  echo converting ${fq} to ${fq%.fastq.gz}.fna.
  zcat ${fqs_dir}/${fq} | awk 'NR%4==1{printf ">%s\n", substr($0,2)}NR%4==2{print}' > ${fas_dir}/${fq%.fastq.gz}.fna
done

exit 0

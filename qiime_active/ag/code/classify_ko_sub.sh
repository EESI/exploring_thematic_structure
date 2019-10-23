#!/bin/bash

for K in 15 25 50 75 100 150
do

  for j in 0 1
  do

    echo running stm $K topics fit $j for kos on training set
    
    rm ~/MiscOut/classify_ko_dump_${K}_${j}.dat
    ~/Dropbox/stm_microbiome/qiime_active/gevers/code/classify_ko_1.R $K $j >> ~/MiscOut/classify_ko_dump_${K}_${j}.dat &
  
  done

done
exit

#!/bin/bash

for i in 1 2 3 4 5 6
do

  for j in 1 2
  do
  
    echo running stan for stm $i fit $j
    
    rm ~/MiscOut/stan_dump_${i}_${j}.dat
    ~/Dropbox/stm_microbiome/qiime_active/gevers/code/stan_1.R $i $j >> ~/MiscOut/stan_dump_${i}_${j}.dat &
    
  done

done
exit

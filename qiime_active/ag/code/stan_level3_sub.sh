#!/bin/bash

for i in 15 25
do

  for j in 0 1
  do
  
    echo running stan for stm $i fit $j
    
    rm ~/MiscOut/stan_dump_${i}_${j}.dat
    ~/Dropbox/stm_microbiome/qiime_active/ag/code/stan_level3_1.R $i $j >> ~/MiscOut/stan_dump_${i}_${j}.dat &
    
  done

done
exit

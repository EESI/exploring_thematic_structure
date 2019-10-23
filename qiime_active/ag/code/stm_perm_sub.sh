#!/bin/bash

for i in 1 2 3 4 5 6
do

  for j in 1 2
  do
  
    echo running stm $i $j topics perm
    
    rm ~/MiscOut/stm_perm_dump_${i}_${j}.dat
    ~/Dropbox/stm_microbiome/qiime_active/ag/code/stm_perm_1.R $i $j >> ~/MiscOut/stm_perm_dump_${i}_${j}.dat &

  done

done
exit

#!/bin/bash

for K in 15 25 50 75 100 150
do
  
  echo running stm $K topics fit for freezing and posterior exploration
  
  rm ~/MiscOut/stm_froze_fit_dump_$K.dat
  ~/Dropbox/stm_microbiome/qiime_active/ag/code/stm_fit_froz_1.R $K >> ~/MiscOut/stm_froze_fit_dump_$K.dat &
    
done
exit

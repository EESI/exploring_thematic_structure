#!/bin/bash

for K in 15 25 50 75 100 150
do

  echo running stm $K topics fit for training set

  rm ~/MiscOut/stm_train_fit_dump_$K.dat
  ~/Dropbox/stm_microbiome/qiime_active/ag/code/stm_fit_train_1.R $K >> ~/MiscOut/stm_train_fit_dump_$K.dat &
  
  done
exit

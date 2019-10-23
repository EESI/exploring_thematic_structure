#!/bin/bash

for K in 15 25 50 75 100 150
do

  echo performing classification for $K topics
  
  rm ~/MiscOut/classify_dump_$K.dat
  ~/Dropbox/stm_microbiome/qiime_active/ag/code/classify_topics_1.R $K >> ~/MiscOut/classify_dump_$K.dat &
  
done
exit

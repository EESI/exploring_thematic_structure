#!/bin/bash

for i in {25,35,50,75,100,150}
do
   ~/Dropbox/stm_microbiome/code_active/gevers_model_prep_rarefied_models.R $i 1000 deseq &
done

exit

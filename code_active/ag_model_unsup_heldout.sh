#!/bin/bash

for i in {25,35,50,75,100,125,150}
do
   ~/Dropbox/stm_microbiome/code_active/ag_model_prep_unsup_models.R ${i} 10000 7 2
done

exit

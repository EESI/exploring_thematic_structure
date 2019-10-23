#!/bin/bash

for i in 1 12 53 74 105
do
   ~/Dropbox/stm_microbiome/code_active/ag_makemodels_traintest_betafeatures.R 35 10000 1 $i &
   ~/Dropbox/stm_microbiome/code_active/ag_makemodels_traintest_betafeatures.R 35 10000 2 $i &
done

exit

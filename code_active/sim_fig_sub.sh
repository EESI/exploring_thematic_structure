#!/bin/bash

rm ~/MiscOut/simfigdump.dat

for IDXcov in 1 2
do
   IDX1=1
   IDX2=0
   while [ $IDX1 -lt 216 ];
   do
      let IDX2=IDX1+9
      echo $IDX1 to $IDX2
      rm ~/MiscOut/simfigdump_$IDXcov_$IDX1.dat
      ~/Dropbox/stm_microbiome/code_active/sim_analysis_looping_info_4.R $IDX1 $IDX2 $IDXcov >> ~/MiscOut/simfigdump_$IDXcov_$IDX1.dat &
      let IDX1=IDX1+10
   done
   echo $IDX1 to 216
   rm ~/MiscOut/simfigdump_$IDXcov_$IDX1.dat
   ~/Dropbox/stm_microbiome/code_active/sim_analysis_looping_info_4.R $IDX1 216 $IDXcov >> ~/MiscOut/simfigdump_$IDXcov_$IDX1.dat &
done

exit

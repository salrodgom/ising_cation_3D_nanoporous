#!/bin/bash
cc zeol_SiGe_distribution_environment_NMR_05.c -o zeol_SiGe_distribution_environment_NMR_05 -lm
for nGe in $(seq 1 1) ; do
 # Analysis:
 InputFile=configuration.${nGe}.txt
 sed -i "s/\.cif//g" $InputFile
 ReportFile=report_${nGe}.txt
 nlines=$(wc -l ${InputFile} | awk '{print $1}')
 echo "${InputFile} ${nlines} ${nGe}
0_0.cif 5
$ReportFile" > main_${nGe}.txt
 ./zeol_SiGe_distribution_environment_NMR_05 main_${nGe}.txt
done

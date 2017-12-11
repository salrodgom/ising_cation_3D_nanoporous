#0 -7719.30894218 2 16 19 31 STW_SiGe_distribution_May2017/stw_Ge_04/stw_scan_04_0009/1d000/conp/stw_scan_04_0009_1d000 .23835468626406336789
#.00000011 -7719.30894207 12 18 20 30 STW_SiGe_distribution_May2017/stw_Ge_04/stw_scan_04_0008/0d985/conp/stw_scan_04_0008_0d985 .23835400711591826055
#.00840390 -7719.30053828 2 16 19 31 STW_SiGe_distribution_May2017/stw_Ge_04/stw_scan_04_0009/1d015/conp/stw_scan_04_0009_1d015 .19172734347508192718
#.07028364 -7719.23865854 2 16 19 31 STW_SiGe_distribution_May2017/stw_Ge_04/stw_scan_04_0009/0d985/conp/stw_scan_04_0009_0d985 .03859829865416564500
#!/bin/bash
for nGe in $(seq 1 59) ; do
 file=configuration.${nGe}.txt
 suma=0.0
 while read line ; do
  energy=$(echo $line | awk '{print $2}')
  boltz=$(echo $line | awk '{print $NF}')
  suma=$(echo "${suma} + ${energy}*${boltz}" | bc -l)
 done < $file
 echo $nGe $suma
done

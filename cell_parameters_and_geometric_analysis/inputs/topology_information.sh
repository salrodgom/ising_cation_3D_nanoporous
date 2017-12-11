#!/bin/bash
cc SiGe_TO_dist_angles_03.c -o SiGe_TO_dist_angles_03 -lm
for nGe in $(seq 1 59) ; do
 file=configuration.${nGe}.txt
 file_topology=configuration_topology.${nGe}.txt
 if [ -f ${file_topology} ] ; then rm ${file_topology} ; fi
 touch ${file_topology}
 while read line ; do
  GeO_distances='#'
  OGeO_angles='#'
  GeOT_angles='#'
  GeGe_distances='#'
  TO_distances='#'
  OTO_angles='#'
  TOT_angles='#'
  TT_distances='#'
  TXTFile=tmp_geom_data.txt
  name=$(echo $line | awk '{print $(NF-1)}')
  name_cif=${name}.cif
  if [ -f $name_cif ] ; then
   echo "$name_cif"
   sed '/F /d' $name_cif > tmp.cif
   sed 's/NameCIFFile/tmp\.cif/g' main_SiGe_TO_dist_angles_00.txt > main
   sed -i 's/PURESILICE/SiO2_STW_c0000\.cif/g' main
   ./SiGe_TO_dist_angles_03 main > salida
   rm tmp.cif main salida
   test1=$(grep -A3 "Overall T (Si or Ge) Analysis" $TXTFile | grep "TO_distances" | tail -n1 | awk '{print $2}')
   test2=$(grep -A3 "Overall T (Si or Ge) Analysis" $TXTFile | grep "TO_distances" | tail -n1 | awk '{print $4}')
   if [ $(echo "$test1 < 1.0" | bc -l) == 1 ] || [ $(echo "$test2 > 4.5" | bc -l) == 1 ] ; then
    TO_distances='-nan'
    echo "[ERROR] with NaN"
    echo "$test1"
    echo "$test2"
   else
#  Ge
    GeO_distances=$(grep -A15 "Overall T (Si or Ge) Analysis" $TXTFile | grep -A8 "Ge Analysis" | grep "GeO_distances" | awk '{print $3}')
    OGeO_angles=$(grep -A15 "Overall T (Si or Ge) Analysis" $TXTFile | grep -A8 "Ge Analysis" | grep "OGeO_angles" | awk '{print $3}')
    GeOT_angles=$(grep -A15 "Overall T (Si or Ge) Analysis" $TXTFile  | grep -A8 "Ge Analysis" | grep "TOT_angles" | awk '{print $3}')
    GeGe_distances=$(grep -A15 "Overall T (Si or Ge) Analysis" $TXTFile | grep -A8 "Ge Analysis" | grep "TT_distances" | awk '{print $3}')
#   Overall
    TO_distances=$(grep -A8 "Overall T (Si or Ge) Analysis" $TXTFile | grep "TO_distances" | awk '{print $3}' )
    OTO_angles=$(grep -A8 "Overall T (Si or Ge) Analysis" $TXTFile | grep "OTO_angles" | awk '{print $3}')
    TOT_angles=$(grep -A8 "Overall T (Si or Ge) Analysis" $TXTFile | grep "TOT_angles" | awk '{print $3}')
    TT_distances=$(grep -A8 "Overall T (Si or Ge) Analysis" $TXTFile | grep "TT_distances" | awk '{print $3}')
#   More:
    TO_T1=$(grep "TO_distances T 1 specific" $TXTFile | awk '{print $6}')
    TO_T2=$(grep "TO_distances T 2 specific" $TXTFile | awk '{print $6}')
    TO_T3=$(grep "TO_distances T 3 specific" $TXTFile | awk '{print $6}')
    TO_T4=$(grep "TO_distances T 4 specific" $TXTFile | awk '{print $6}')
    TO_T5=$(grep "TO_distances T 5 specific" $TXTFile | awk '{print $6}')
    TOT_T1=$(grep "TOT_angles T 1 specific" $TXTFile | awk '{print $6}')
    TOT_T2=$(grep "TOT_angles T 2 specific" $TXTFile | awk '{print $6}')
    TOT_T3=$(grep "TOT_angles T 3 specific" $TXTFile | awk '{print $6}')
    TOT_T4=$(grep "TOT_angles T 4 specific" $TXTFile | awk '{print $6}')
    TOT_T5=$(grep "TOT_angles T 5 specific" $TXTFile | awk '{print $6}')
   fi
   if [ "$OGeO_angles" == "-nan" ] || [ "$GeOT_angles" == "-nan" ] || [ "$TOT_angles" == "-nan" ] || [ "$OTO_angles" == "-nan" ] ; then
    echo "[ERROR] with NaN"
    #exit 0
   fi
   echo "$TO_distances   $GeO_distances   $TT_distances   $GeGe_distances   $OTO_angles   $OGeO_angles   $TOT_angles   $GeOT_angles " >> ${file_topology}
  else
   echo $name_cif NO EXISTE
  fi
 done < $file
done
rm SiGe_TO_dist_angles_03
exit 0

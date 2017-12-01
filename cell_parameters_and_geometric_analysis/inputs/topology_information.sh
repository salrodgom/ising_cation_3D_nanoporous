#!/bin/bash
cc SiGe_TO_dist_angles_01.c -o SiGe_TO_dist_angles_01 -lm
for nGe in $(seq 1 59) ; do
 file=configuration.${nGe}.txt
 file_topology=configuration_topology.${nGe}.txt
 if [ -f ${file_topology} ] ; then rm ${file_topology} ; fi
 echo "# GeO_distances OGeO_angles GeGe_distances TO_distances TT_distances TOT_angles TO_T1 TO_T2 TO_T3 TO_T4 TO_T5 OTO_T1 OTO_T2 OTO_T3 OTO_T4 OTO_T5 name_cif" > ${file_topology}
 while read line ; do
  TXTFile=tmp_geom_data.txt
  name=$(echo $line | awk '{print $4}')
  name_cif=${name}.cif
  cp $name_cif tmp.cif
  echo "$name_cif"
  sed 's/NameCIFFile/tmp\.cif/g' main_SiGe_TO_dist_angles_00.txt > main
  sed -i 's/PURESILICE/SiO2_STW_c0000\.cif/g' main
  ./SiGe_TO_dist_angles_01 main > salida
  rm tmp.cif main salida
# Ge
  GeO_distances=$(grep -A15 "Overall T (Si or Ge) Analysis" $TXTFile | grep -A5 "Ge Analysis" | grep "GeO_distances" | awk '{print $3}')
  OGeO_angles=$(grep -A15 "Overall T (Si or Ge) Analysis" $TXTFile | grep -A5 "Ge Analysis" | grep "OGeO_angles" | awk '{print $3}')
  GeGe_distances=$(grep -A15 "Overall T (Si or Ge) Analysis" $TXTFile | grep -A5 "Ge Analysis" | grep "TT_distances" | awk '{print $3}')
# Overall
  TO_distances=$(grep -A15 "Overall T (Si or Ge) Analysis" $TXTFile | grep "TO_distances" | awk '{print $3}' )
  TT_distances=$(grep -A15 "Overall T (Si or Ge) Analysis" $TXTFile | grep -A5 "Ge Analysis" | grep "TT_distances" | awk '{print $3}')
  TOT_angles=$(grep -A15 "Overall T (Si or Ge) Analysis" $TXTFile | grep -A5 "Ge Analysis" | grep "TOT_angles" | awk '{print $3}')
# More:
  TO_T1=$(grep "TO_distances T 1 specific" $TXTFile | awk '{print $6}')
  TO_T2=$(grep "TO_distances T 2 specific" $TXTFile | awk '{print $6}')
  TO_T3=$(grep "TO_distances T 3 specific" $TXTFile | awk '{print $6}')
  TO_T4=$(grep "TO_distances T 4 specific" $TXTFile | awk '{print $6}')
  TO_T5=$(grep "TO_distances T 5 specific" $TXTFile | awk '{print $6}')
  OTO_T1=$(grep "OTO_angles T 1 specific" $TXTFile | awk '{print $6}')
  OTO_T2=$(grep "OTO_angles T 2 specific" $TXTFile | awk '{print $6}')
  OTO_T3=$(grep "OTO_angles T 3 specific" $TXTFile | awk '{print $6}')
  OTO_T4=$(grep "OTO_angles T 4 specific" $TXTFile | awk '{print $6}')
  OTO_T5=$(grep "OTO_angles T 5 specific" $TXTFile | awk '{print $6}')
  echo $GeO_distances $OGeO_angles $GeGe_distances $TO_distances $TT_distances $TOT_angles $TO_T1 $TO_T2 $TO_T3 $TO_T4 $TO_T5 $OTO_T1 $OTO_T2 $OTO_T3 $OTO_T4 $OTO_T5 >> ${file_topology}
 done < $file
done

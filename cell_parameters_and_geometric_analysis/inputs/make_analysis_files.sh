#!/bin/bash 
kB=8.61734E-5     # eV/K
temperature=448.0 # K
kT=$(echo "${kB}*${temperature}" | sed 's/[eE]+*/*10^/g' | bc -l)
rm analysis_*.txt
for nGe in $(seq 1 59) ; do
 nGeName=$(echo $nGe | awk '{ printf("%02d\n", $1) }')
 ReportFile=report_${nGe}.txt
 if [ -f $ReportFile ] ; then
  grep -A32 "configuration summary information" $ReportFile | grep "relative lattice energy" | awk '{print $1}'  > energy
  Z_partition=0.0
  while read line ; do
   Z_partition=$(echo "${Z_partition} + e(-(${line})/(${kT}))" | bc -l)
  done < energy
  if [ -f BoltzmannWeight ]; then rm -rf BoltzmannWeight ; touch BoltzmannWeight ; fi
  while read line ; do
   if [ ! -z $line ] ; then
    echo $(echo "scale=40; (1.0/(${Z_partition}))*e(-(${line})/(${kT}))" | bc -l) >> BoltzmannWeight
   fi
  done < energy
  echo "$Z_partition"
  grep -A38 "configuration summary information" $ReportFile | grep "0  number of Ge 0 by D4R" | awk '{print $1}' > 0
  grep -A38 "configuration summary information" $ReportFile | grep "1  number of Ge 1 by D4R" | awk '{print $1}' > 1
  grep -A38 "configuration summary information" $ReportFile | grep "2  number of Ge 2 by D4R" | awk '{print $1}' > 2
  grep -A38 "configuration summary information" $ReportFile | grep "3  number of Ge 2 by D4R" | awk '{print $1}' > 3
  grep -A38 "configuration summary information" $ReportFile | grep "4  number of Ge 2 by D4R" | awk '{print $1}' > 4
  grep -A38 "configuration summary information" $ReportFile | grep "5  number of Ge 2 by D4R" | awk '{print $1}' > 5
  grep -A38 "configuration summary information" $ReportFile | grep "6  number of Ge 3 by D4R" | awk '{print $1}' > 6
  grep -A38 "configuration summary information" $ReportFile | grep "7  number of Ge 3 by D4R" | awk '{print $1}' > 7
  grep -A38 "configuration summary information" $ReportFile | grep "8  number of Ge 3 by D4R" | awk '{print $1}' > 8
  grep -A38 "configuration summary information" $ReportFile | grep "9  number of Ge 3 by D4R" | awk '{print $1}' > 9
  grep -A38 "configuration summary information" $ReportFile | grep "10 number of Ge 4 by D4R" | awk '{print $1}' > 10
  grep -A38 "configuration summary information" $ReportFile | grep "11 number of Ge 4 by D4R" | awk '{print $1}' > 11
  grep -A38 "configuration summary information" $ReportFile | grep "12 number of Ge 4 by D4R" | awk '{print $1}' > 12
  grep -A38 "configuration summary information" $ReportFile | grep "13 number of Ge 4 by D4R" | awk '{print $1}' > 13
  grep -A38 "configuration summary information" $ReportFile | grep "14 number of Ge 4 by D4R" | awk '{print $1}' > 14
  grep -A38 "configuration summary information" $ReportFile | grep "15 number of Ge 4 by D4R" | awk '{print $1}' > 15
  grep -A38 "configuration summary information" $ReportFile | grep "16 number of Ge 5 by D4R" | awk '{print $1}' > 16
  grep -A38 "configuration summary information" $ReportFile | grep "17 number of Ge 5 by D4R" | awk '{print $1}'   > 17
  grep -A38 "configuration summary information" $ReportFile | grep "18 number of Ge 5 by D4R" | awk '{print $1}'   > 18
  grep -A38 "configuration summary information" $ReportFile | grep "19 number of 5 by D4R as 5" | awk '{print $1}' > 19
  grep -A38 "configuration summary information" $ReportFile | grep "19 number of Ge 6 by D4R" | awk '{print $1}'   > 20
  grep -A38 "configuration summary information" $ReportFile | grep "20 number of Ge 6 by D4R as forming a square and 2 direct neighbours"   | awk '{print $1}' > 21
  grep -A38 "configuration summary information" $ReportFile | grep "20 number of Ge 6 by D4R as forming a square and 2 Ge in face diagonal" | awk '{print $1}' > 22
  grep -A38 "configuration summary information" $ReportFile | grep "21 numbre of Ge 6 by D4R" | awk '{print $1}' > 23
  grep -A38 "configuration summary information" $ReportFile | grep "22 number of Ge 7 by D4R" | awk '{print $1}' > 24
  grep -A38 "configuration summary information" $ReportFile | grep "23 number of Ge 8 by D4R" | awk '{print $1}' > 25
  grep -A38 "configuration summary information" $ReportFile | grep "27 number of Ge at crystallographic site 1" | awk '{print $1}' > 26
  grep -A38 "configuration summary information" $ReportFile | grep "28 number of Ge at crystallographic site 2" | awk '{print $1}' > 27
  grep -A38 "configuration summary information" $ReportFile | grep "29 number of Ge at crystallographic site 3" | awk '{print $1}' > 28
  grep -A38 "configuration summary information" $ReportFile | grep "30 number of Ge at crystallographic site 4" | awk '{print $1}' > 29
  grep -A38 "configuration summary information" $ReportFile | grep "Current File Name" | awk '{print $1}' > name
  #if [ -f NameUpdated ] ; then rm NameUpdated ; touch NameUpdated ; fi
  #while read line ; do
  # grep $line dictionary.txt | awk '{print $3}' >> NameUpdated
  #done < name
  nGeD4R=0
  paste 26 27 28 29 > tmp.All
  if [ -f tmp.T5 ] ; then rm tmp.T5 ; touch tmp.T5 ; fi
  while read line ; do
   nT1=$(echo $line | awk '{print $1}')
   nT2=$(echo $line | awk '{print $2}')
   nT3=$(echo $line | awk '{print $3}')
   nT4=$(echo $line | awk '{print $4}')
   let "nT5 = $nGe - $nT1 - $nT2 - $nT3 - $nT4"
   echo $nT5 >> tmp.T5
  done < tmp.All
  mv tmp.T5 30
  rm tmp.All
  paste BoltzmannWeight 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 name > analysis_${nGe}.txt
  rm    BoltzmannWeight 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 name
  #
  for i in $(seq 0 30) ; do n[${i}]=0 ; done
  while read line ; do
   p=$(echo $line | awk '{print $1}')
   for i in $(seq 0 30) ; do
    n[${i}]=$(echo "${n[${i}]} + $(echo $line | awk -v col=$(($i+2)) '{print $col}') * $p" | bc -l)
   done
  done < analysis_${nGe}.txt
  echo "Average:" $(echo "${n[@]}") >> analysis_${nGe}.txt
  echo "#Average:" > analysis_${nGe}_col.txt
  for i in $(seq 0 25) ; do
   echo ${i} ${n[${i}]} >> analysis_${nGe}_col.txt
  done
  echo $nGe ${n[26]} >> analysis_T1.txt
  echo $nGe ${n[27]} >> analysis_T2.txt
  echo $nGe ${n[28]} >> analysis_T3.txt
  echo $nGe ${n[29]} >> analysis_T4.txt
  echo $nGe ${n[30]} >> analysis_T5.txt
  #
 fi
done

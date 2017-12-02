#!/bin/bash
PI=$(echo "scale=20; 4*a(1)" | bc -l)
kB=8.61734E-5     # eV/K
temperature=448.0 # K
kT=$(echo "${kB}*${temperature}" | sed 's/[eE]+*/*10^/g' | bc -l)
#
function abs () {
 if [ $(echo "$1 < 0.0" | bc -l) == 1 ] ; then
  echo "0.0 - $1" | bc -l
 else
  echo $1
 fi
}
if [ -f list_CIFFiles ] ; then rm list_CIFFiles ; fi
ls STW_SiGe_distribution_May2017/stw_Ge_*/stw_scan_*_*/*/conp/stw_scan_*_*_*.cif >  list_CIFFiles
ls SiGe_distribution_January_2017/*.*/*/conp/*.cif                               >> list_CIFFiles
if [ -f configurations.txt ] ; then rm configurations.txt ; fi
touch configurations.txt
while read line ; do
 name=$(echo $line | sed 's/\.cif//g')
 name_gout=${name}.gout
 name_cif=${name}.cif
 if [ ! -f $name_gout ] ; then
  echo "WARNING: $name_gout doesn't exits" 
 else
  for i in $(seq 1 59) ; do Ge[${i}]="" ; done
  energy=$(grep "l e" $name_gout | awk '{print $4}')
  sed -n -e '/occ/,$p' $name_cif | sed '/occ/d' | grep -n 'Ge' | sed 's/://g' | awk '{print $1}' > configuration
  n_Ge=$(wc -l configuration | awk '{print $1}')
  i=0
  while read list ; do
   let i++
   Ge[${i}]=$(echo $list | awk '{print $1}')
  done < configuration
  rm configuration
  echo $n_Ge $energy ${Ge[@]} $name >> configurations.txt
 fi
done < list_CIFFiles
for nGe in $(seq 1 59) ; do
 file=configuration.${nGe}.txt
 file_revised=configuration_revised.${nGe}.txt
 awk -v nGe=${nGe} '{if ($1==nGe) print $0}' configurations.txt | awk '{out=""; for(i=2;i<=NF;i++){out=out" "$i}; print out}' | sort -gk1 | uniq > $file
 if [ -f ${file_revised} ] ; then rm ${file_revised} ; fi
 touch $file_revised
 previous_energy=99999.9999
 while read line ; do
  energy=$(echo $line | awk '{print $1}')
  if [ $( echo "$( abs `echo "(${energy}) - (${previous_energy})" | bc -l` ) > 0.0000001 " | bc -l) == 1  ] ; then
   echo $line >> $file_revised
   previous_energy=${energy}
  fi
 done < $file
 #
 energy_min=$(head -n1 ${file_revised} | awk '{print $1}') 
 energy_min=`abs ${energy_min}`
 Z_partition=0.0
 while read line ; do
  energy=$(echo $line | awk '{print $1}')
  Z_partition=$(echo "$Z_partition + e(-(${energy}+${energy_min})/${kT})" | bc -l)
 done < ${file_revised}
 if [ -f tmp ] ; then rm tmp ; fi
 touch tmp
 while read line ; do
  energy=$(echo $line | awk '{print $1}')
  let n_configurations++
  echo "$(echo "scale=6; $energy + ${energy_min}" | bc -l) $line $(echo "(1.0/(${Z_partition}))*e(-(${energy}+${energy_min})/${kT})" | bc -l)" >> tmp
 done < ${file_revised}
 #
 mv tmp ${file_revised}
done

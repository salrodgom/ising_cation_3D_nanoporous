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
declare -a Ge
ls *.gout > list
while read line ; do
 name=$(echo $line | sed 's/\.gout//g')
 cifname=${name}"_optimised.cif"
 if [ -f $cifname ] ; then
  energy=$(grep "l e" $line | awk '{print $4}')
  sed -n -e '/occ/,$p' $cifname | sed '/occ/d' | grep -n 'Ge' | sed 's/://g' | awk '{print $1}' > configuration
  n_Ge=0
  while read list ; do
   let n_Ge++
   Ge[${n_Ge}]=$(echo $list)
  done < configuration
  rm configuration
  echo $energy ${Ge[@]} $name >> configurations.txt
 fi
done < list
sort  -gk1 configurations.txt | uniq > c
mv c  configurations.txt
touch configurations_revised.txt
previous_energy=99999.9999
while read line ; do
 energy=$(echo $line | awk '{print $1}')
 if [ $( echo "$( abs `echo "(${energy}) - (${previous_energy})" | bc -l` ) > 0.000001 " | bc -l) == 1  ] ; then
  echo $line >> configurations_revised.txt
  previous_energy=${energy}
 fi
done < configurations.txt
mv configurations_revised.txt configurations.txt
sort -gk1 configurations.txt | uniq > c
energy_min=$(head -n1 configurations.txt | awk '{print $1}')
energy_min=`abs ${energy_min}`
Z_partition=0.0
while read line ; do
 energy=$(echo $line | awk '{print $1}')
 Z_partition=$(echo "$Z_partition + e(-(${energy}+${energy_min})/${kT})" | bc -l)
done < configurations.txt 
#Z_partition=$(awk -v kT=${kB} -v energy_min=${energy_min} '{Z+=exp(-($1-$energy_min)/$kT)} END {print Z}' configurations.txt )
n_configurations=0
while read line ; do
 energy=$(echo $line | awk '{print $1}')
 let n_configurations++
 echo "$(echo "scale=6; $energy + ${energy_min}" | bc -l) $line $(echo "(1.0/(${Z_partition}))*e(-(${energy}+${energy_min})/${kT})" | bc -l)" >> configurations_revised.txt
done < configurations.txt
exit 0

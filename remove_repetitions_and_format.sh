#!/bin/bash
PI=$(echo "scale=20; 4*a(1)" | bc -l)
kB=8.61734E-5     # eV/K
temperature=448.0 # K
kT=$(echo "${kB}*${temperature}" | sed 's/[eE]+*/*10^/g' | bc -l)
MC_steps=3000
#
function abs () {
 if [ $(echo "$1 < 0.0" | bc -l) == 1 ] ; then
  echo "0.0 - $1" | bc -l
 else
  echo $1
 fi
}
for molar_fraction in $(seq 3 59) ; do
  sed -i 's/\.cif//g' configuration.${molar_fraction}.txt
  # elimino repeticiones
  sort -gk1 configuration.${molar_fraction}.txt | uniq > c
  touch configurations.txt
  previous_energy=99999.9999
  while read line ; do
   energy=$(echo $line | awk '{print $2}')
   if [ $( echo "$( abs `echo "(${energy}) - (${previous_energy})" | bc -l` ) > 0.000001 " | bc -l) == 1  ] ; then
    echo $line >> configurations.txt
    previous_energy=${energy}
    echo $line
   fi
  done < c
  #
  mv configurations.txt configuration.${molar_fraction}.txt
  touch zero.configuration.${molar_fraction}.txt
  rm *.gin 
  energy_min=$(head -n1 configuration.${molar_fraction}.txt | awk '{print $1}')
  Z_partition=0
  while read line ; do
   energy=$(echo $line | awk '{print $1}')
   Z_partition=$(echo "$Z_partition + e(-($energy)/$kT)" | bc -l)
  done < configuration.${molar_fraction}.txt
  #Z_partition=$(awk -v kT=${kB} -v energy_min=${energy_min} '{Z+=exp(-($1-$energy_min)/$kT)} END {print Z}' configuration.${molar_fraction}.txt )
  n_configurations=0
  while read line ; do
   energy=$(echo $line | awk '{print $1}')
   let n_configurations++
   echo "$line  $(echo "(1.0/(${Z_partition}))*e(-(${energy})/${kT})" | bc -l)" >> zero.configuration.${molar_fraction}.txt
  done < configuration.${molar_fraction}.txt
  #mv zero.configuration.${molar_fraction}.txt configuration.${molar_fraction}.txt
done

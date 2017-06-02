#!/bin/bash -x
PI=$(echo "scale=20; 4*a(1)" | bc -l)
kB=8.61734E-5     # eV/K
temperature=448.0 # K
kT=$(echo "${kB}*${temperature}" | sed 's/[eE]+*/*10^/g' | bc -l)
function abs () {
 if [ $(echo "$1 < 0.0" | bc -l) == 1 ] ; then
  echo "0.0 - $1" | bc -l
 else
  echo $1
 fi
}
function go {
ni=$(ps aux | grep 'ising' | sed '/grep/d' | wc -l)
while [ $(echo "$ni >= 4" | bc -l) == 1 ] ; do
 sleep 10
 ni=$(ps aux | grep 'ising' | sed '/grep/d' | wc -l)
done
./ising_frameworks < input > ${molar_fraction}.output
}
make install
if [ ! -d CALCS ] ; then mkdir CALCS ; fi
for molar_fraction in $(seq 1 2) ; do
  echo ${molar_fraction}
  echo "${molar_fraction}" > input
  go
  name=$(grep 'filename:' ${molar_fraction}.output | awk '{print $2}')
  cat ${name} src/oxygen.gin > ${molar_fraction}.gin
  rm ${name}
  # elimino repeticiones
  cp configuration.2.txt configurations.txt
  sort  configurations.txt -gk1 | uniq > c
  rm -rf configurations.txt
  touch configurations.txt
  while read line ; do
   energy=$(echo $line | awk '{print $1}')
   n=0
   while read line2 ; do
     energy_reference=$(echo ${line2} | awk '{print $1}')
     delta=$(echo "(${energy}) - (${energy_reference})" | bc -l)
     delta=$(abs $delta)
     if [ $(echo "$delta <= 0.0001" | bc -l) == 1 ] ; then
      let n++
     fi
   done < configurations.txt
   if [ $(echo "$n == 0" | bc -l) == 1 ] ; then
    echo $line
    echo $line >> configurations.txt
   fi
  done < c
  #
  mv configurations.txt configuration.${molar_fraction}.txt
  touch zero.configuration.${molar_fraction}.txt
  rm ${molar_fraction}.gin 

  energy_min=$(head -n1 configuration.${molar_fraction}.txt | awk '{print $1}')
  Z_partition=0.0
  while read line ; do
   energy=$(echo $line | awk '{print $1}')
   Z_partition=$(echo "${Z_partition} + e(-((${energy})-(${energy_min}))/(${kT}))" | bc -l)
  done < configuration.${molar_fraction}.txt
  while read line ; do
   energy=$(echo $line | awk '{print $1}')
   echo "$line $(echo "(1.0/(${Z_partition}))*e(-((${energy})-(${energy_min}))/(${kT}))" | bc -l)" >> zero.configuration.${molar_fraction}.txt
  done < configuration.${molar_fraction}.txt

  mv zero.configuration.${molar_fraction}.txt CALCS/.
done


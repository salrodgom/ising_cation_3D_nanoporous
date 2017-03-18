#!/bin/bash -x
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
for molar_fraction in $(seq 1 20) ; do
  echo ${molar_fraction}
  echo "${molar_fraction}" > input
  go
  name=$(grep 'filename:' ${molar_fraction}.output | awk '{print $2}')
  cat ${name} src/oxygen.gin > ${molar_fraction}.gin
  rm ${name}
  sort configurations.txt -nk1 | uniq > c
  mv c configuration.${molar_fraction}.txt
  mv configuration.${molar_fraction}.txt ${molar_fraction}.output ${molar_fraction}.gin CALCS
done

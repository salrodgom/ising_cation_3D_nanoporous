#!/bin/bash
mkdir 1st
for nGe in $(seq 1 59) ; do
 file=configuration.${nGe}.txt
 head -n1 $file > c
 while read line ; do
  name=$(echo $line | awk '{print $(NF-1)}')
  name_cif=${name}.res
  if [ -f $name_cif ] ; then
   cp $name_cif 1st/.
  else
   echo $name_cif NO EXISTE
  fi
 done < c
done

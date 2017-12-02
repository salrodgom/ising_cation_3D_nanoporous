#!/bin/bash
PI=$(echo "scale=20; 4*a(1)" | bc -l)
kB=8.61734E-5     # eV/K
temperature=448.0 # K
kT=$(echo "${kB}*${temperature}" | sed 's/[eE]+*/*10^/g' | bc -l)
# pure silica:
cp  c0000.gin                   tmp.gin
sed -i 's/CONFIGURATION/0_0/g'  tmp.gin
gulp < tmp.gin > /dev/null
rm tmp.gin
# 
for nGe in $(seq 1 59) ; do
 file=configuration.${nGe}.txt
 if [ -f $file ] ; then
  sed -i "s/\.cif//g" $file
  if [ -f reduced.$file ] ; then rm reduced.$file ; fi
  touch reduced.$file
  while read line ; do
   #echo $line
   Relative=$(echo $line | awk -v col=1 '{print $col}')
   Absolute=$(echo $line | awk -v col=2 '{print $col}')
   GeList=$(echo $line | awk '{for (i = 3; i <= NF-2; i++) printf $i " "; print ""}')
   let "n_col = $nGe + 3"
   CIFFile=$(echo $line | awk -v col=$n_col '{print $col}')
   BoltzmannWeight=$(echo $line | awk '{print $NF}')
   if [ $( echo "$BoltzmannWeight > 0.00000000001" | bc -l ) == 1 ] ; then
    echo $line
    cp  c0000.gin                          tmp.gin
    sed -i "s/CONFIGURATION/${CIFFile}/g"  tmp.gin
    for i in $(seq 1 $nGe) ; do
     GePosition=$(echo $GeList | awk -v col=$i '{print $col+4}')
     sed -i "${GePosition}s/Si   core/Ge   core/g" tmp.gin
    done
    gulp < tmp.gin > /dev/null
    mv ${CIFFile}.cif CIFFiles/.
    echo $Relative $Absolute $GeList $CIFFile >> reduced.$file
   fi
  done < $file
  Z=0.0
  Z_partition=0
  while read line ; do
   Relative=$(echo $line | awk '{print $1}')
   Z_partition=$(echo "$Z_partition + e(-($Relative)/$kT)" | bc -l)
  done < reduced.$file
  echo $Z_partition
  if [ -f tmp ] ; then rm tmp ; fi
  touch tmp
  while read line ; do
   Absolute=$(echo $line | awk '{print $1}')
   echo "$line  $(echo "(1.0/(${Z_partition}))*e(-(${Absolute})/${kT})" | bc -l)" >> tmp
  done < reduced.$file
  mv tmp reduced.$file
 fi
done
cp  c0000.gin                    tmp.gin
sed -i 's/CONFIGURATION/60_1/g'  tmp.gin
sed -i "s/Si   core/Ge   core/g" tmp.gin
gulp < tmp.gin > /dev/null
rm tmp.gin

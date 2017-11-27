#!/bin/bash
# pure silica:
cp  c0000.gin                   tmp.gin
sed -i 's/CONFIGURATION/0_0/g'  tmp.gin
gulp < tmp.gin > /dev/null
rm tmp.gin
# 
for nGe in $(seq 4 4) ; do
 file=configuration.${nGe}.txt
 if [ -f $file ] ; then
  sed -i "s/\.cif//g" $file
  while read line ; do
   echo $line
   echo ==================
   Relative=$(echo $line | awk -v col=1 '{print $col}')
   Absolute=$(echo $line | awk -v col=2 '{print $col}')
   GeList=$(echo $line | awk '{for (i = 3; i <= NF-3; i++) printf $i " "; print ""}')
   let "n_col = $nGe + 4"
   CIFFile=$(echo $line | awk -v col=$n_col '{print $col}')
   BoltzmannWeight=$(echo $line | awk '{print $NF}')
   # check:
   echo $Relative
   echo $Absolute
   echo $GeList
   echo $CIFFile
   echo $BoltzmannWeight
   #
   cp  c0000.gin                          tmp.gin
   sed -i "s/CONFIGURATION/${CIFFile}/g"  tmp.gin
   for i in $(seq 1 $nGe) ; do
    GePosition=$(echo $GeList | awk -v col=$i '{print $col+4}')
    sed -i "${GePosition}s/Si   core/Ge   core/g" tmp.gin
   done
   gulp < tmp.gin > /dev/null
  done < $file
 fi
done
cp  c0000.gin                   tmp.gin
sed -i 's/CONFIGURATION/60_1/g'  tmp.gin
sed -i "s/Si   core/Ge   core/g" tmp.gin
gulp < tmp.gin > /dev/null
rm tmp.gin

#!/bin/bash
awk '{printf "%3s %3s %3s %3s %6s\n", $3,$4,$5,$6,$1}' OUTSOD > strings
awk '{printf "%3s %3s %3s %3s %6s\n", $3,$4,$5,$6,$1}' OUTSOD_modified >> strings
awk '{if( $1 <= -7718.0 ) print $1,$2,$3,$4,$5}' cosa > list 
if [ -f found_list ]; then rm found_list ; fi
touch found_list
while read line ; do
 a=$(echo $line | awk '{print $2}')
 b=$(echo $line | awk '{print $3}')
 c=$(echo $line | awk '{print $4}')
 d=$(echo $line | awk '{print $5}')
 energy=$(echo $line | awk '{print $1}')
 found=$(grep " $a " strings | grep " $b " | grep " $c " | grep " $d ")
 echo "$energy $found " >> found_list
done < list
uniq found_list | sort -gk1 > c
mv c found_list
rm list strings
i=0
while read line ; do
 let i++
 echo "0.0 $line 4_${i} 0"
done < found_list

#!/bin/bash
if [ -f gp ] ; then rm gp ; fi
if [ -f all.txt ] ; then rm all.txt ; fi
echo "set term postscript eps color enhanced blacktext 'Helvetica,24'
set output 'view_config.eps'
set yrange [0:1]
set xrange [0:20]
set ylabel 'energy / eV'
set xlabel 'configuration / -'
set key outside" > gp
echo "plot 0.0 w l lc rgb 'black' notitle, \\" >> gp
for i in $(seq 1 12) ; do
 if [ -f  zero.configuration.${i}.txt ] ; then
  rm     zero.configuration.${i}.txt
  touch  zero.configuration.${i}.txt
 fi
 energy_min=$(head -n1 configuration.${i}.txt | awk '{print $1}')
 cat configuration.${i}.txt | while read line ; do
  energy=$(echo $line | awk '{print $1}')
  echo "$(echo "${energy} - ${energy_min}" | bc -l) $line" >> zero.configuration.${i}.txt
  echo "$(echo "${energy} - ${energy_min}" | bc -l) $line" >> all.txt
 done
 echo " 'zero.configuration.${i}.txt' u 0:1 w lp pt 6 title '${i}', \\" >> gp
done

echo "0 lc 8 notitle">> gp

gnuplot gp

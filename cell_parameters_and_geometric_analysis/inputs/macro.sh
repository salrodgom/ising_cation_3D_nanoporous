#!/bin/bash
for i in $(seq 1 60) ; do
 echo $i $(head -n5 configuration_topology.${i}.txt | tail -n1 | awk '{print $2}')
done

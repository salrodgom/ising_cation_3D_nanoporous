#!/bin/bash
for i in $(seq 1 59) ; do
 ln -s ./"done"/configuration_revised.${i}.txt configuration.${i}.txt
done

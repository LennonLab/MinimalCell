#!/bin/bash
# Tested using bash version 4.1.5
for ((i=1;i<=96;i++)); 
do 
   cp *B_${i}_S*fastq 3B_${i}/
   echo ${i}
done

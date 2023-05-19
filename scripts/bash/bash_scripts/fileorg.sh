
#!/bin/bash
# Tested using bash version 4.1.5
for ((i=2;i<=96;i++)); 
do 
   cp *1_${i}_S*fastq s1_${i}/
   echo ${i}
done

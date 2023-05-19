
#!/bin/bash
# Tested using bash version 4.1.5
cd /N/dc2/projects/muri2/roy_20200302/MA_s1_GD2VCF
for ((i=1;i<=96;i++)); 
do 
   bgzip < s1_${i}_clean.VCF > gz_versions/s1_${i}_clean.VCF.gz
done
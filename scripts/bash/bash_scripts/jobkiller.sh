
#!/bin/bash
# Tested using bash version 4.1.5
for ((i=1628415;i<=1628510;i++)); 
do 
   qdel ${i}
   echo ${i}
done

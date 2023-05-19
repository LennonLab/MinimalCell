#!/bin/bash

cd /N/dc2/projects/muri2/roy/MA_3B/

#AR=(3B_U 3B_V 3B_W 3B_X 3B_Y)

#for i in "${AR[@]}"
for ((i=1;i<=96;i++));
do
	module load breseq
	echo "cd /N/dc2/projects/muri2/roy/MA_3B/" > 3B_BreSeq_${i}
	echo "" >> 3B_BreSeq_${i}
	echo "module load breseq" >> 3B_BreSeq_${i}
	echo "" >> 3B_BreSeq_${i}
	#echo "breseq -j 8 -p -o 3B_${i}_breseq -r /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.gb Sample_${i}/Sample_${i}_R1_trimmed.fastq" >> 3B_BreSeq_${i}
	echo "breseq -j 8 -p -o 3B_${i}_breseq -r /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.gb 3B_${i}/3B_${i}_R1_trimmed.fastq 3B_${i}/3B_${i}_R2_trimmed.fastq" >> 3B_BreSeq_${i}
	qsub -l walltime=20:00:00,vmem=100gb,nodes=1:ppn=8 ./3B_BreSeq_${i}
	done
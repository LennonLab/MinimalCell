#!/bin/bash

cd /N/dc2/projects/muri2/roy/Mm_NSE

AR=(syn1_1 syn1_2 syn1_3 syn1_4 syn1_5 syn1_6 syn1_7_t.160.clone syn1_8_anc)

for i in "${AR[@]}"

do
	module load breseq
	echo "cd /N/dc2/projects/muri2/roy/Mm_NSE" > Mm_BreSeq_${i}
	echo "" >> Mm_BreSeq_${i}
	echo "module load breseq" >> Mm_BreSeq_${i}
	echo "" >> Mm_BreSeq_${i}
	echo "breseq -j 8 -p -o Sample${i}_breseq -r /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.gb Sample_${i}/Sample_${i}_R1_trimmed.fastq Sample_${i}/Sample_${i}_R2_trimmed.fastq" >> Mm_BreSeq_${i}
	qsub -l walltime=20:00:00,vmem=100gb,nodes=1:ppn=8 ./Mm_BreSeq_${i}
	done
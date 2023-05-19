#!/bin/bash

cd /N/dc2/projects/muri2/roy/MA_3B/
#mkdir MA_3B_outsub/
#AR=(3B_U 3B_V 3B_W 3B_X 3B_Y)

#for i in "${AR[@]}"
for ((i=1;i<=96;i++));
do
	module load breseq
	gdtools subtract -o ../MA_3B_outsub/3B_${i}_sub_annotated.gd 3B_${i}_breseq/output/evidence/annotated.gd MA_3B_anc_annotated.gd
	gdtools GD2VCF -o ../MA_3B_GD2VCF/3B_${i}_clean.VCF ../MA_3B_outsub/3B_${i}_sub_annotated.gd -r /N/dc2/projects/muri2/roy_20200302/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.gb
	done
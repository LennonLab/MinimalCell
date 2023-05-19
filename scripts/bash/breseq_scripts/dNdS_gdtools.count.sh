cd /N/project/MinimalCell/codons
module load breseq

gdtools count -v -r Synthetic.bacterium_JCVI-Syn3A.gb -o syn3_count.csv -s -b -p /N/project/MinimalCell/codons/mm9/output.gd

gdtools count -v -r Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.gb -o syn1_count.csv -s -b -p /N/project/MinimalCell/codons/mm1/output.gd
#!/bin/bash

#$ -q rcc-30d

# SamMuri_2Y
date
cd /N/dc2/projects/muri2/roy/Mm_NSE
AR=(3B_10)

for i in "${AR[@]}"
do
    module load gatk/3.8; module load samtools; module load cutadapt; module load fastqc; module load bwa; module load picard
	cd Sample_${i}
	cat *R1_001.fastq > Sample_${i}_R1.fastq
	cat *R2_001.fastq > Sample_${i}_R2.fastq	
	
	# clean and trim
	echo "#!/bin/bash" > Mm_NSE_${i}_qc_clean.sh
	echo "" >> Mm_NSE_${i}_qc_clean.sh
	echo "#$ -l vmem=50gb walltime=48:00:00" >> Mm_NSE_${i}_qc_clean.sh
	echo "" >> Mm_NSE_${i}_qc_clean.sh
	echo "module load samtools; module load cutadapt; module load fastqc; module load bwa; module load java; module load picard" >> Mm_NSE_${i}_qc_clean.sh
	echo ""
	echo "cd /N/dc2/projects/muri2/roy/Mm_NSE/Sample_${i}" >> Mm_NSE_${i}_qc_clean.sh
	echo "" >> Mm_NSE_${i}_qc_clean.sh
	echo "time cutadapt -u 5 -o Sample_${i}_R1_trimmed.fastq Sample_${i}_R1.fastq" >> Mm_NSE_${i}_qc_clean.sh
	echo "" >> Mm_NSE_${i}_qc_clean.sh
	echo "time cutadapt -u 5 -o Sample_${i}_R2_trimmed.fastq Sample_${i}_R2.fastq" >> Mm_NSE_${i}_qc_clean.sh
	echo "" >> Mm_NSE_${i}_qc_clean.sh
	echo "qsub -l walltime=20:00:00,vmem=64gb,nodes=1:ppn=4 Mm_NSE_${i}_Assemble.sh" >> Mm_NSE_${i}_qc_clean.sh
	echo "" >> Mm_NSE_${i}_qc_clean.sh
	echo "mkdir fastqc" >> Mm_NSE_${i}_qc_clean.sh
	echo "" >> Mm_NSE_${i}_qc_clean.sh
	echo "fastqc -o fastqc/ Sample_${i}_R1_trimmed.fastq" >> Mm_NSE_${i}_qc_clean.sh
	echo "" >> Mm_NSE_${i}_qc_clean.sh
	echo "fastqc -o fastqc/ Sample_${i}_R2_trimmed.fastq" >> Mm_NSE_${i}_qc_clean.sh
	echo "" >> Mm_NSE_${i}_qc_clean.sh
	echo "exit" >> Mm_NSE_${i}_qc_clean.sh
	
	
#### Assembly ####

	echo "#!/bin/bash" > Mm_NSE_${i}_Assemble.sh
	echo "" >> Mm_NSE_${i}_Assemble.sh
	echo "#$ -l vmem=50gb walltime=48:00:00 " >> Mm_NSE_${i}_Assemble.sh
	echo "" >> Mm_NSE_${i}_Assemble.sh
	echo "module load samtools; module load cutadapt; module load fastqc; module load bwa" >> Mm_NSE_${i}_Assemble.sh
	echo ""
	echo "cd /N/dc2/projects/muri2/roy/Mm_NSE/Sample_${i}/" >> Mm_NSE_${i}_Assemble.sh
	echo "" >> Mm_NSE_${i}_Assemble.sh
	echo "echo \"BWA\" >&2" >> Mm_NSE_${i}_Assemble.sh
	echo "time bwa aln /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta Sample_${i}_R1_trimmed.fastq > Sample_${i}_R1_aln.sai" >> Mm_NSE_${i}_Assemble.sh
	echo ""
	echo "time bwa aln /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta Sample_${i}_R2_trimmed.fastq > Sample_${i}_R2_aln.sai" >> Mm_NSE_${i}_Assemble.sh
    echo ""
	echo "time bwa sampe /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta Sample_${i}_R1_aln.sai Sample_${i}_R2_aln.sai Sample_${i}_R1_trimmed.fastq Sample_${i}_R2_trimmed.fastq > Sample_${i}.sam" >> Mm_NSE_${i}_Assemble.sh
   	echo "" >> Mm_NSE_${i}_Assemble.sh
	echo "qsub -l walltime=2:00:00,vmem=20gb Mm_NSE_${i}_compress.fastq.sh" >> Mm_NSE_${i}_Assemble.sh
	echo "" >> Mm_NSE_${i}_Assemble.sh
	echo "echo \"SAMTOOLS\" >&2" >> Mm_NSE_${i}_Assemble.sh
	echo "samtools view -bS -T /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta Sample_${i}.sam > Sample_${i}.bam" >> Mm_NSE_${i}_Assemble.sh
	echo "" >> Mm_NSE_${i}_Assemble.sh
	echo "samtools sort -o Sample_${i}.sorted.bam Sample_${i}.bam" >> Mm_NSE_${i}_Assemble.sh
	echo "" >> Mm_NSE_${i}_Assemble.sh
	echo "samtools index Sample_${i}.sorted.bam" >> Mm_NSE_${i}_Assemble.sh
	echo "" >> Mm_NSE_${i}_Assemble.sh
	echo "" >> Mm_NSE_${i}_Assemble.sh
	
	        ### GATK ###	
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/picard/picard-tools-1.107/AddOrReplaceReadGroups.jar I=Sample_${i}.sorted.bam O=Sample_${i}.sorted.fixed.bam SORT_ORDER=coordinate RGID=GeneConv RGLB=bar RGPL=illumina RGSM=Sample_${i} RGPU=6 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT" >> Mm_NSE_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/picard/picard-tools-1.107/MarkDuplicates.jar I=Sample_${i}.sorted.fixed.bam O=Sample_${i}.sorted.fixed.marked.bam M=Sample_${i}.metrics CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT" >> Mm_NSE_${i}_Assemble.sh
	echo "GATK" >&2
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta -I Sample_${i}.sorted.fixed.marked.bam -T RealignerTargetCreator -o Sample_${i}.intervals" >> Mm_NSE_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -R /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta -I Sample_${i}.sorted.fixed.marked.bam -T IndelRealigner -targetIntervals Sample_${i}.intervals --filter_bases_not_stored -o Sample_${i}.sorted.fixed.marked.realigned.bam" >> Mm_NSE_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 -T UnifiedGenotyper -R /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta -I Sample_${i}.sorted.fixed.marked.realigned.bam -glm BOTH -rf BadCigar -o Sample_${i}.sorted.fixed.marked.realigned.vcf" >> Mm_NSE_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -T BaseRecalibrator -I Sample_${i}.sorted.fixed.marked.realigned.bam -R /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta -rf BadCigar --filter_bases_not_stored -knownSites Sample_${i}.sorted.fixed.marked.realigned.vcf -o Sample_${i}.recal_data.grp" >> Mm_NSE_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta -I Sample_${i}.sorted.fixed.marked.realigned.bam -T PrintReads -rf BadCigar  -o Sample_${i}.mapped.bam -BQSR Sample_${i}.recal_data.grp" >> Mm_NSE_${i}_Assemble.sh
	echo "" >> Mm_NSE_${i}_Assemble.sh
	echo "exit" >> Mm_NSE_${i}_Assemble.sh	
	
 	# compress.fastq again
	echo "#!/bin/bash" > Mm_NSE_${i}_compress.fastq.sh
	echo "" >> Mm_NSE_${i}_compress.fastq.sh
	echo "#$ -l vmem=50gb walltime=48:00:00 " >> Mm_NSE_${i}_compress.fastq.sh
	echo "module load samtools; module load cutadapt; module load fastqc; module load bwa; module load java; module load picard" >> Mm_NSE_${i}_compress.fastq.sh
	echo ""
	echo "" >> Mm_NSE_${i}_compress.fastq.sh
	echo "cd /N/dc2/projects/muri2/roy/Mm_NSE/Sample_${i}/" >> Mm_NSE_${i}_compress.fastq.sh
	echo "rm Sample_${i}_R*_filtered.fastq" >> Mm_NSE_${i}_compress.fastq.sh
	echo "rm Sample_${i}_R*_rmadapter.fastq" >> Mm_NSE_${i}_compress.fastq.sh
	echo "" >> Mm_NSE_${i}_compress.fastq.sh
	echo "exit" >> Mm_NSE_${i}_compress.fastq.sh
	
	
		
	chmod u+x Mm_NSE_${i}_qc_clean.sh
	chmod u+x Mm_NSE_${i}_Assemble.sh
	chmod u+x Mm_NSE_${i}_compress.fastq.sh
		
	qsub -l walltime=2:00:00,vmem=20gb Mm_NSE_${i}_qc_clean.sh
	cd ..
	done
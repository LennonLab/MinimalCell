#!/bin/bash

#$ -q rcc-30d

# SamMuri_2Y
date
cd /N/dc2/projects/muri2/roy/MA_3B/
#AR=(94 95)

for ((i=1;i<=96;i++));
do
    module load gatk/3.8; module load samtools; module load cutadapt; module load fastqc; module load bwa; module load picard
	cd 3B_${i}
	cat *R1_001.fastq > 3B_${i}_R1.fastq
	cat *R2_001.fastq > 3B_${i}_R2.fastq
	
	# clean and trim
	echo "#!/bin/bash" > MA_3B_${i}_qc_clean.sh
	echo "" >> MA_3B_${i}_qc_clean.sh
	echo "#$ -l vmem=50gb walltime=2:00:00" >> MA_3B_${i}_qc_clean.sh
	echo "" >> MA_3B_${i}_qc_clean.sh
	echo "module load samtools; module load cutadapt; module load fastqc; module load bwa; module load java; module load picard" >> MA_3B_${i}_qc_clean.sh
	echo "" >> MA_3B_${i}_qc_clean.sh
	echo "cd /N/dc2/projects/muri2/roy/MA_3B/3B_${i}" >> MA_3B_${i}_qc_clean.sh
	echo "" >> MA_3B_${i}_qc_clean.sh
	echo "time cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o 3B_${i}_R1_rmadapter.fastq -p 3B_${i}_R2_rmadapter.fastq 3B_${i}_R1.fastq 3B_${i}_R2.fastq" >> MA_3B_${i}_qc_clean.sh
	echo "" >> MA_3B_${i}_qc_clean.sh
	echo "time cutadapt -u 5 -o 3B_${i}_R1_trimmed.fastq 3B_${i}_R1_rmadapter.fastq" >> MA_3B_${i}_qc_clean.sh
	echo "" >> MA_3B_${i}_qc_clean.sh
	echo "time cutadapt -u 5 -o 3B_${i}_R2_trimmed.fastq 3B_${i}_R2_rmadapter.fastq" >> MA_3B_${i}_qc_clean.sh
	echo "" >> MA_3B_${i}_qc_clean.sh
	echo "qsub -l walltime=2:00:00,vmem=64gb,nodes=1:ppn=4 MA_3B_${i}_Assemble.sh" >> MA_3B_${i}_qc_clean.sh
	echo "" >> MA_3B_${i}_qc_clean.sh
	echo "mkdir fastqc" >> MA_3B_${i}_qc_clean.sh
	echo "" >> MA_3B_${i}_qc_clean.sh
	echo "fastqc -o fastqc/ 3B_${i}_R1_trimmed.fastq" >> MA_3B_${i}_qc_clean.sh
	echo "" >> MA_3B_${i}_qc_clean.sh
	echo "fastqc -o fastqc/ 3B_${i}_R2_trimmed.fastq" >> MA_3B_${i}_qc_clean.sh
	echo "" >> MA_3B_${i}_qc_clean.sh
	echo "exit" >> MA_3B_${i}_qc_clean.sh
	
#### Assembly ####

	echo "#!/bin/bash" > MA_3B_${i}_Assemble.sh
	echo "" >> MA_3B_${i}_Assemble.sh
	echo "#$ -l vmem=50gb walltime=24:00:00 " >> MA_3B_${i}_Assemble.sh
	echo "" >> MA_3B_${i}_Assemble.sh
	echo "module load samtools; module load cutadapt; module load fastqc; module load bwa" >> MA_3B_${i}_Assemble.sh
	echo "" >> MA_3B_${i}_Assemble.sh
	echo "cd /N/dc2/projects/muri2/roy/MA_3B/3B_${i}/" >> MA_3B_${i}_Assemble.sh
	echo "" >> MA_3B_${i}_Assemble.sh
	echo "echo \"BWA\" >&2" >> MA_3B_${i}_Assemble.sh
	echo "" >> MA_3B_${i}_Assemble.sh
	echo "time bwa aln /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta 3B_${i}_R1_trimmed.fastq > 3B_${i}_R1_aln.sai" >> MA_3B_${i}_Assemble.sh
	echo "" >> MA_3B_${i}_Assemble.sh
	echo "time bwa aln /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta 3B_${i}_R2_trimmed.fastq > 3B_${i}_R2_aln.sai" >> MA_3B_${i}_Assemble.sh
    echo "" >> MA_3B_${i}_Assemble.sh
	echo "time bwa sampe /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta 3B_${i}_R1_aln.sai 3B_${i}_R2_aln.sai 3B_${i}_R1_trimmed.fastq 3B_${i}_R2_trimmed.fastq > 3B_${i}.sam" >> MA_3B_${i}_Assemble.sh
   	echo "" >> MA_3B_${i}_Assemble.sh
	echo "qsub -l walltime=1:00:00,vmem=20gb MA_3B_${i}_compress.fastq.sh" >> MA_3B_${i}_Assemble.sh
	echo "" >> MA_3B_${i}_Assemble.sh
	echo "echo \"SAMTOOLS\" >&2" >> MA_3B_${i}_Assemble.sh
	echo "samtools view -bS -T /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta 3B_${i}.sam > 3B_${i}.bam" >> MA_3B_${i}_Assemble.sh
	echo "" >> MA_3B_${i}_Assemble.sh
	echo "samtools sort -o 3B_${i}.sorted.bam 3B_${i}.bam" >> MA_3B_${i}_Assemble.sh
	echo "" >> MA_3B_${i}_Assemble.sh
	echo "samtools index 3B_${i}.sorted.bam" >> MA_3B_${i}_Assemble.sh
	echo "" >> MA_3B_${i}_Assemble.sh
	echo "" >> MA_3B_${i}_Assemble.sh
	
	        ### GATK ###	
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/picard/picard-tools-1.107/AddOrReplaceReadGroups.jar I=3B_${i}.sorted.bam O=3B_${i}.sorted.fixed.bam SORT_ORDER=coordinate RGID=GeneConv RGLB=bar RGPL=illumina RGSM=3B_${i} RGPU=6 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT" >> MA_3B_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/picard/picard-tools-1.107/MarkDuplicates.jar I=3B_${i}.sorted.fixed.bam O=3B_${i}.sorted.fixed.marked.bam M=3B_${i}.metrics CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT" >> MA_3B_${i}_Assemble.sh
	echo "GATK" >&2
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta -I 3B_${i}.sorted.fixed.marked.bam -T RealignerTargetCreator -o 3B_${i}.intervals" >> MA_3B_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -R /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta -I 3B_${i}.sorted.fixed.marked.bam -T IndelRealigner -targetIntervals 3B_${i}.intervals --filter_bases_not_stored -o 3B_${i}.sorted.fixed.marked.realigned.bam" >> MA_3B_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 -T UnifiedGenotyper -R /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta -I 3B_${i}.sorted.fixed.marked.realigned.bam -glm BOTH -rf BadCigar -o 3B_${i}.sorted.fixed.marked.realigned.vcf" >> MA_3B_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -T BaseRecalibrator -I 3B_${i}.sorted.fixed.marked.realigned.bam -R /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta -rf BadCigar --filter_bases_not_stored -knownSites 3B_${i}.sorted.fixed.marked.realigned.vcf -o 3B_${i}.recal_data.grp" >> MA_3B_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.fasta -I 3B_${i}.sorted.fixed.marked.realigned.bam -T PrintReads -rf BadCigar  -o 3B_${i}.mapped.bam -BQSR 3B_${i}.recal_data.grp" >> MA_3B_${i}_Assemble.sh
	echo "" >> MA_3B_${i}_Assemble.sh
	echo "exit" >> MA_3B_${i}_Assemble.sh	
	
 	# compress.fastq again
	echo "#!/bin/bash" > MA_3B_${i}_compress.fastq.sh
	echo "" >> MA_3B_${i}_compress.fastq.sh
	echo "#$ -l vmem=50gb walltime=4:00:00 " >> MA_3B_${i}_compress.fastq.sh
	echo "module load samtools; module load cutadapt; module load fastqc; module load bwa; module load java; module load picard" >> MA_3B_${i}_compress.fastq.sh
	echo ""
	echo "" >> MA_3B_${i}_compress.fastq.sh
	echo "cd /N/dc2/projects/muri2/roy/MA_3B/3B_${i}/" >> MA_3B_${i}_compress.fastq.sh
	#echo "rm 3B_${i}_R*_filtered.fastq" >> MA_3B_${i}_compress.fastq.sh
	echo "rm 3B_${i}_R*_rmadapter.fastq" >> MA_3B_${i}_compress.fastq.sh
	echo "" >> MA_3B_${i}_compress.fastq.sh
	echo "exit" >> MA_3B_${i}_compress.fastq.sh
	
	
		
	chmod u+x MA_3B_${i}_qc_clean.sh
	chmod u+x MA_3B_${i}_Assemble.sh
	chmod u+x MA_3B_${i}_compress.fastq.sh

		
	qsub -l walltime=2:00:00,vmem=20gb MA_3B_${i}_qc_clean.sh
	cd ..
	done
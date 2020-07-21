#!/bin/bash

DATA_PATH=/home/yudiz/data/peanut_simu

conda activate /home/yudiz/bin

for p in 0.005 0.01 0.02
do
	PARENT_FOLDER="$DATA_PATH/homr${p}"
	cat $PARENT_FOLDER/refA.fa $PARENT_FOLDER/refB.fa > $PARENT_FOLDER/ref_all.fsa
	bwa index $PARENT_FOLDER/ref_all.fsa
		
	# picard does not like .fsa extension (weird)
	mv $PARENT_FOLDER/refA.fsa $PARENT_FOLDER/refA.fa
	mv $PARENT_FOLDER/refB.fsa $PARENT_FOLDER/refB.fa

	# preparing index files for alignment
	bwa index $PARENT_FOLDER/refA.fa
	bwa index $PARENT_FOLDER/refB.fa
	samtools faidx $PARENT_FOLDER/refA.fa
	samtools faidx $PARENT_FOLDER/refB.fa

	# prpeparing dictionary file needed by GATK
	java -jar /usr/local/jar/picard.jar CreateSequenceDictionary REFERENCE=$PARENT_FOLDER/refA.fa OUTPUT=$PARENT_FOLDER/refA.dict
	java -jar /usr/local/jar/picard.jar CreateSequenceDictionary REFERENCE=$PARENT_FOLDER/refB.fa OUTPUT=$PARENT_FOLDER/refB.dict
	
	for q in 3 4 8 12 16
	do
		DATA_FOLDER="$PARENT_FOLDER/cov${q}"

		for m in {0..49}
		do
			# split the reads to A and B
			bwa mem $PARENT_FOLDER/ref_all.fsa $DATA_FOLDER/sim${m}.fq > $DATA_FOLDER/sim${m}test.sam
			samtools view -S -b $DATA_FOLDER/sim${m}test.sam > $DATA_FOLDER/sim${m}test.bam
			samtools sort $DATA_FOLDER/sim${m}test.bam -o $DATA_FOLDER/sim${m}index.test.bam
			samtools index $DATA_FOLDER/sim${m}index.test.bam
			samtools view -b $DATA_FOLDER/sim${m}index.test.bam Genome_B > $DATA_FOLDER/sim${m}in_chr1.bam
			samtools view -b $DATA_FOLDER/sim${m}index.test.bam Genome_A > $DATA_FOLDER/sim${m}in_chr0.bam

			bamToFastq -i $DATA_FOLDER/sim${m}in_chr1.bam -fq $DATA_FOLDER/sim${m}B.fastq
			bamToFastq -i $DATA_FOLDER/sim${m}in_chr0.bam -fq $DATA_FOLDER/sim${m}A.fastq

			bwa mem $PARENT_FOLDER/refA.fa $DATA_FOLDER/sim${m}A.fastq > $DATA_FOLDER/sim${m}A.sam
			bwa mem $PARENT_FOLDER/refB.fa $DATA_FOLDER/sim${m}B.fastq > $DATA_FOLDER/sim${m}B.sam
			
			# sam to bam conversion
			samtools view -bS $DATA_FOLDER/sim${m}A.sam > $DATA_FOLDER/sim${m}A.bam
			samtools view -bS $DATA_FOLDER/sim${m}B.sam > $DATA_FOLDER/sim${m}B.bam

			# Adding @RG to header
			java -jar /usr/local/jar/picard.jar AddOrReplaceReadGroups I=$DATA_FOLDER/sim${m}A.bam O=$DATA_FOLDER/sim${m}A-RG.bam RGID=sim0 RGPL=illumina RGPU=uni1 RGSM=sim0 RGLB=lib1
			java -jar /usr/local/jar/picard.jar AddOrReplaceReadGroups I=$DATA_FOLDER/sim${m}B.bam O=$DATA_FOLDER/sim${m}B-RG.bam RGID=sim0 RGPL=illumina RGPU=uni1 RGSM=sim0 RGLB=lib1

			# Sort bam files
			samtools sort $DATA_FOLDER/sim${m}A-RG.bam -o $DATA_FOLDER/sim${m}A.sorted.bam
			samtools sort $DATA_FOLDER/sim${m}B-RG.bam -o $DATA_FOLDER/sim${m}B.sorted.bam

			# index sorted bam files
			samtools index $DATA_FOLDER/sim${m}A.sorted.bam
			samtools index $DATA_FOLDER/sim${m}B.sorted.bam

			# Call SNPs using GATK
			/usr/local/gatk-4.1.7.0-40-gc65be43-SNAPSHOT/gatk HaplotypeCaller -R $PARENT_FOLDER/refA.fa -I $DATA_FOLDER/sim${m}A.sorted.bam -O $DATA_FOLDER/sim${m}A.vcf
			/usr/local/gatk-4.1.7.0-40-gc65be43-SNAPSHOT/gatk HaplotypeCaller -R $PARENT_FOLDER/refB.fa -I $DATA_FOLDER/sim${m}B.sorted.bam -O $DATA_FOLDER/sim${m}B.vcf

			# haplotyping using hapcut2
			extractHAIRS --bam $DATA_FOLDER/sim${m}A.sorted.bam --VCF $DATA_FOLDER/sim${m}A.vcf --out $DATA_FOLDER/sim${m}A_fragment
			HAPCUT2 --fragments $DATA_FOLDER/sim${m}A_fragment --VCF $DATA_FOLDER/sim${m}A.vcf --output $DATA_FOLDER/sim${m}Ahap_file --outvcf 1
			extractHAIRS --bam $DATA_FOLDER/sim${m}B.sorted.bam --VCF $DATA_FOLDER/sim${m}B.vcf --out $DATA_FOLDER/sim${m}B_fragment
			HAPCUT2 --fragments $DATA_FOLDER/sim${m}B_fragment --VCF $DATA_FOLDER/sim${m}B.vcf --output $DATA_FOLDER/sim${m}Bhap_file --outvcf 1
		done
	done
done

conda deactivate

#bgzip haplotype_output_file.phased.VCF
#tabix haplotype_output_file.phased.vcf.gz
#cat refA.fa | bcftools consensus haplotype_output_file.phased.vcf.gz > outA.fa


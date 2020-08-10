#!/bin/bash

REFERENCE=~/Documents/Karin/alleotetraploid/data/hmm/WGS/reference/ref9.fsa
OUTPUT_SIMU=~/Documents/peanut_simu
ERR_FILE1=~/Documents/Karin/alleotetraploid/data/roshan/simu/miseq250R1.txt
ERR_FILE2=~/Documents/Karin/alleotetraploid/data/roshan/simu/miseq250R2.txt

for p in 0.003 0.005 0.007
do
	PARENT_FOLDER="$OUTPUT_SIMU/ale${p}"
#echo "$PARENT_FOLDER"
	mkdir "$PARENT_FOLDER"
	SEED=$RANDOM
#echo "$SEED"
	~/Documents/practice/run_simu -e $REFERENCE -j .003 -g ${p} -a 100 -b 100 -s $SEED -o $PARENT_FOLDER/indiv -f $PARENT_FOLDER/ref -r $PARENT_FOLDER/ref.sam -n 50 &> $PARENT_FOLDER/truth.txt
	bwa index $PARENT_FOLDER/refA.fsa
	bwa index $PARENT_FOLDER/refB.fsa
	mkdir "$PARENT_FOLDER/pair"
	mkdir "$PARENT_FOLDER/single"
	for q in 4 8 12 16
	do
		DATA_FOLDER="$PARENT_FOLDER/pair/cov${q}"
		mkdir "$DATA_FOLDER"
		DATA_FOLDER_SG="$PARENT_FOLDER/single/cov${q}"
		mkdir "$DATA_FOLDER_SG"
		mkdir "$DATA_FOLDER/hmm_res"
		mkdir "$DATA_FOLDER_SG/hmm_res"
		for m in {0..49}
		do
			#PAIR-END
			art_illumina -1 $ERR_FILE1 -2 $ERR_FILE2 -l 150 -i $PARENT_FOLDER/indiv${m}.fsa -o $DATA_FOLDER/sim${m} -f ${q} -p -m 200 -s 10
		 	bwa mem $PARENT_FOLDER/refA.fsa $DATA_FOLDER/sim${m}1.fq $DATA_FOLDER/sim${m}2.fq > $DATA_FOLDER/aln${m}A.sam
			bwa mem $PARENT_FOLDER/refB.fsa $DATA_FOLDER/sim${m}1.fq $DATA_FOLDER/sim${m}2.fq > $DATA_FOLDER/aln${m}B.sam
			#SINGLE-END
			art_illumina -1 $ERR_FILE1 -l 150 -i $PARENT_FOLDER/indiv${m}.fsa -o $DATA_FOLDER_SG/sim${m} -f ${q}
			bwa mem $PARENT_FOLDER/refA.fsa $DATA_FOLDER_SG/sim${m}.fq > $DATA_FOLDER_SG/aln${m}A.sam
			bwa mem $PARENT_FOLDER/refB.fsa $DATA_FOLDER_SG/sim${m}.fq > $DATA_FOLDER_SG/aln${m}B.sam
		done
	done
done



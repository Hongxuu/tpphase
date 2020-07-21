#!/bin/bash

REFERENCE=~/Documents/Karin/alleotetraploid/data/hmm/WGS/reference/ref9.fsa
OUTPUT_SIMU=~/Documents/Karin/alleotetraploid/data/hmm/WGS
ERR_FILE=~/Documents/Karin/alleotetraploid/data/roshan/simu/miseq250R2.txt

for p in 0.005 0.01 0.02
do
	PARENT_FOLDER="$OUTPUT_SIMU/homr${p}"
#echo "$PARENT_FOLDER"
	mkdir "$PARENT_FOLDER"
	SEED=$RANDOM
#echo "$SEED"
	~/Documents/practice/run_simu -e $REFERENCE -j ${p} -g .003 -a 100 -b 100 -s $SEED -o $PARENT_FOLDER/indiv -f $PARENT_FOLDER/ref -r $PARENT_FOLDER/ref.sam -n 50 &> $PARENT_FOLDER/truth.txt
	bwa index $PARENT_FOLDER/refA.fsa
	bwa index $PARENT_FOLDER/refB.fsa
	for q in 3 4 8 12 16
	do
		DATA_FOLDER="$PARENT_FOLDER/cov${q}"
		mkdir "$DATA_FOLDER"
		for m in {0..49}
		do
			echo "$PARENT_FOLDER/indiv${m}.fsa"
			/Users/yudi/art_bin_MountRainier/art_illumina -1 $ERR_FILE -l 150 -i $PARENT_FOLDER/indiv${m}.fsa -o $DATA_FOLDER/sim${m} -f ${q}
		 	bwa mem $PARENT_FOLDER/refA.fsa $DATA_FOLDER/sim${m}.fq > $DATA_FOLDER/aln${m}A.sam
			bwa mem $PARENT_FOLDER/refB.fsa $DATA_FOLDER/sim${m}.fq > $DATA_FOLDER/aln${m}B.sam
		done
	done
done



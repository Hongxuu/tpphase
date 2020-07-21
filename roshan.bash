#!/bin/bash

EXE_PATH=/home/peanut/peanut2/WGS/simulation_V4/ksd-roshan_wgs
DATA_PATH=/home/yudiz/data/peanut_simu

for p in 0.005 0.01 0.02
do
	PARENT_FOLDER="$DATA_PATH/homr${p}"
	for q in 3 4 8 12 16
	do
		DATA_FOLDER="$PARENT_FOLDER/cov${q}"
		for m in {0..49}
		do
			$EXE_PATH/roshan --geno $PARENT_FOLDER/ref.sam --fsa_files $PARENT_FOLDER/refA.fsa $PARENT_FOLDER/refB.fsa --sam_files $DATA_FOLDER/aln${m}A $DATA_FOLDER/aln${m}B -j $DATA_FOLDER/extracted${m} --ref_names Genome_A:0-5000 Genome_B:0-5000 &> $DATA_FOLDER/roshan_res${m}
		done
	done
done

#!/bin/bash

fastqdir=$1
oligos=$2


for f in `ls $fastqdir/*_R1_001.fastq`
do 
	r=`echo $f|sed -e 's/_R1_/_R2_/'`
	echo sbatch -D $PWD --output $PWD/make_contigs_%j.out --export=ALL -J make_contigs -c 4 --mem 1G -p high --wrap="sh make_contigs.sh $f $r $oligos"
	sbatch -D $PWD --output $PWD/make_contigs_%j.out --export=ALL -J make_contigs -c 4 --mem 1G -p high --wrap="sh make_contigs.sh $f $r $oligos"
done

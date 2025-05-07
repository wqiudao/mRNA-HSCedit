#!/bin/bash
#reads_mapping_bwa_pe fq1 fq2 /path/to
[ "$#" -lt 2 ] && echo "usage: reads_clean_pe fq1[.gz] fq2[.gz] [/path/to]" && exit 0
echo $1;echo $2;echo $3;file_1=$1;file_2=$2;
suffix='.fq';
if [ "${file_1##*.}" = "gz" ]; then
	suffix='.fq.gz';
	file_1=${file_1%.*};file_2=${file_2%.*};
fi
file_1=${file_1%.*};file_2=${file_2%.*};
data_1=$1;data_2=$2;
if [ "$3" != "" ]; then
	data_1=$3/$1;
	data_2=$3/$2;
fi
echo $data_1;echo $data_2;
echo $file_1;echo $file_2;
trimmomatic PE -threads 3 -phred33 $data_1 $data_2  ${file_1}_clean_paired$suffix ${file_1}_unpaired$suffix ${file_2}_clean_paired$suffix ${file_2}_unpaired$suffix ILLUMINACLIP:/home/miniconda3/miniconda3/envs/metaspades/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:0:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:50:20 MINLEN:50 > ${file_1}_${file_2}_clean.log 2>&1

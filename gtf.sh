#!/bin/sh
############  Notice  #################
#Step1: Assembly of novel transcripts by Reference Annotation Based Transcript (RABT) method using Cufflinks in each sample based on RNA-seq reads and reference annotation
echo -e "
====================================================
[$(date +%R:%S)] Assembly of novel transcripts by Reference Annotation Based Transcript (RABT) method\n";

############  Parameters  #################
#Reading defined parameters
TSSfolder=$0;TSSfolder=${TSSfolder:0:${#TSSfolder}-6};
original_data=$1;outpath=$2;
public_gtf=$3
MultiProcessor=$4

func=gtf
codepath=${2}/code/
outpath=${2}/tmp/${func}/
mkdir $codepath
mkdir $outpath

############  De novo assembly step  #################
#Loading all replicates in Group A and B
list=($(echo ${original_data} | tr "," "\n"))

#For each replicate using RABT method to assemble
j=1
for i in ${list[*]}
do
mkdir ${outpath}${j}/
cufflinks -N -p ${MultiProcessor} -g ${public_gtf} -o ${outpath}${j} ${i}
j=$[ $j + 1 ]
done
#Output: GTF format; Exons of transcripts annotation in each row





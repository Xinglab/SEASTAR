#!/bin/sh
############  Notice  #################
#Step1: Assembly of novel transcripts by Reference Annotation Based Transcript (RABT) method using Cufflinks in each sample based on RNA-seq reads and reference annotation
echo "Assembly of novel transcripts by Reference Annotation Based Transcript (RABT) method";

############  Parameters  #################
#Reading defined parameters
TSSfolder=$0;TSSfolder=${TSSfolder:0:${#TSSfolder}-12};
original_data=$1;outpath=$2;
public_gtf=$3
MultiProcessor=$4

func=gtf
codepath=${2}/code/
outpath=${2}/tmp/${func}/
mkdir $codepath
mkdir $outpath
#Current_path record
Current_path=($(pwd))

############  De novo assembly step  #################
#Loading all replicates in Group A and B
list=($(echo ${original_data} | tr "," "\n"))

#For each replicate using RABT method to assemble
j=1
for i in ${list[*]}
do
mkdir ${outpath}${j}/

#Auto detection and change the relative path into absolute path for GTF
tmp=${i%/*};
if [ x${tmp} != "x" ]
then
gtffolder=($(cd ${tmp} ; pwd)/)
else
gtffolder=($(pwd)/)
fi
i=${gtffolder}${i##*/}
cd ${Current_path}


cat ${TSSfolder}header.txt > ${codepath}a00_${func}_${j}.sh
cat>>${codepath}a00_${func}_${j}.sh<<EOF
cufflinks -N -p ${MultiProcessor} -g ${public_gtf} -o ${outpath}${j} ${i}
EOF
j=$[ $j + 1 ]
done
#Output: GTF format; Exons of transcripts annotation in each row






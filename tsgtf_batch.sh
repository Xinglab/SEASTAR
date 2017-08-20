#!/bin/sh
############  Notice  #################
#Step2-2: Reformat the annotation of transcriptome and an interface for more FPKM analysis
echo "Reformat the annotation of transcriptome";

############  Parameters  #################
#Reading defined parameters
TSSfolder=$0;TSSfolder=${TSSfolder:0:${#TSSfolder}-14};
original_data=$1;outpath=$2;
public_gtf=$3
MultiProcessor=$4

func=tsgtf
codepath=${2}/code/
outpath=${2}/tmp/${func}/
inpath=${2}/tmp/mrg/
mkdir $outpath

############  Reformating step  #################
tmp=($(echo ${original_data} | tr "," "\n"))
samtools view -H ${tmp[0]} > ${outpath}tmp.sam
samtools view ${tmp[0]} | head >> ${outpath}tmp.sam
samtools view -Sb ${outpath}tmp.sam > ${outpath}tmp.bam
cat ${TSSfolder}header.txt > ${codepath}a02_${func}.sh
cat>>${codepath}a02_${func}.sh<<EOF

cufflinks -N -p ${MultiProcessor} -G ${public_gtf} -o ${outpath} ${outpath}tmp.bam

cut -d" " -f1-6 ${outpath}transcripts.gtf | tr ' ' '\t' | tr -d ';' |tr -d '"'| cut -f1,3,4,5,7,10,12,14 > ${outpath}transcripts_cut8.gtf
EOF

#Output: GTF format; Exons of transcripts annotation in each row





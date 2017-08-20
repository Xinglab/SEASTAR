#!/bin/sh
############  Notice  #################
#Step4: Count the reads coverage in each replicate; including the exon body region, its downstream splice junction, downstream intron region and upstream intergenic region, as well as the coverage used for Dapars
echo -e "
====================================================
[$(date +%R:%S)] Count the reads coverage in each replicate\n";

############  Parameters  #################
#Reading defined parameters
TSSfolder=$0;TSSfolder=${TSSfolder:0:${#TSSfolder}-10};
original_data=$1;outpath=$2;genome=$3;strand=$4;

func=cov_utr
codepath=${2}/code/
outpath=${2}/tmp/${func}/
inpath=${2}/
mkdir $outpath

if [ $strand != "U" ]
then
  strand=-${strand};
fi
if [ $strand = "U" ]
then
  strand=;
fi


############  Counting step  #################
#Loading all replicates in Group A and B
list=($(echo ${original_data} | tr "," "\n"))

#For each replicate counting the reads covered the exon body region, its downstream splice junction, downstream intron region and upstream intergenic region
j=1
for i in ${list[*]}
do
mkdir ${outpath}${j}/

coverageBed ${strand} -abam ${i} -b ${inpath}utr_split.bed -split | tr -d "\r" > ${outpath}${j}/tsscov.cov
#coverageBed ${strand} -abam ${i} -b ${inpath}utr_split.bed -split | tr -d "\r" > ${outpath}${j}/tsscov_u.cov
#coverageBed ${strand} -abam ${i} -b ${inpath}utr_split.bed -split | tr -d "\r" > ${outpath}${j}/tsscov_d.cov
samtools view -H ${i} > ${outpath}${j}/junction
samtools view ${i} | awk -F"\t" '$6~/N/' >> ${outpath}${j}/junction
samtools view -Sb ${outpath}${j}/junction > ${outpath}${j}/junction.b
rm ${outpath}${j}/junction
coverageBed ${strand} -abam ${outpath}${j}/junction.b -b ${inpath}utr_split.bed -split | tr -d "\r" > ${outpath}${j}/tsscov_j.cov
#genomeCoverageBed -bg -ibam ${i} -g ${genome} -split > ${outpath}${j}/coverage.bedgraph


j=$[ $j + 1 ]

done

#Output: tsscov.cov, tsscov_j.cov, tsscov_u.cov, tsscov_d.cov and coverage.bedgraph
#Reads counts covered the exon body region, its downstream splice junction, downstream intron region and upstream intergenic region, as well as the coverage used for Dapars





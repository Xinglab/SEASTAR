#!/bin/sh
############  Notice  #################
#Step4: Count the reads coverage in each replicate; including the exon body region, its downstream splice junction, downstream intron region and upstream intergenic region, as well as the coverage used for Dapars
echo "Count the reads coverage in each replicate";

############  Parameters  #################
#Reading defined parameters
TSSfolder=$0;TSSfolder=${TSSfolder:0:${#TSSfolder}-12};
original_data=$1;outpath=$2;genome=$3;strand=$4;

func=cov
codepath=${2}/code/
outpath=${2}/tmp/${func}/
inpath=${2}/tmp/tsgtf/
#Current_path record
Current_path=($(pwd))
mkdir $outpath

if [ $strand != "U" ]
then
  strand=-${strand};
fi
if [ $strand = "U" ]
then
  strand=;
fi

################  bedtools version detection  ################
## bedtools has a important change of coverageBed as version 2.24.0
vercomp () {
    if [[ $1 == $2 ]]
    then
        return 0
    fi
    local IFS=.
    local i ver1=($1) ver2=($2)
    # fill empty fields in ver1 with zeros
    for ((i=${#ver1[@]}; i<${#ver2[@]}; i++))
    do
        ver1[i]=0
    done
    for ((i=0; i<${#ver1[@]}; i++))
    do
        if [[ -z ${ver2[i]} ]]
        then
            # fill empty fields in ver2 with zeros
            ver2[i]=0
        fi
        if ((10#${ver1[i]}>10#${ver2[i]}))
        then
            return 1
        fi
        if ((10#${ver1[i]}<10#${ver2[i]}))
        then
            return 2
        fi
    done
    return 0
}
bedver=($(bedtools -version |tr -d "bedtools v"))
    vercomp $bedver 2.24.0
    case $? in
        0) bedA='-b'
        bedB='-a';; # "="
        1) bedA='-b'
        bedB='-a';; # ">"
        2) bedA='-abam'
        bedB='-b';; # "<"
    esac


############  Counting step  #################
#Loading all replicates in Group A and B
list=($(echo ${original_data} | tr "," "\n"))

#For each replicate counting the reads covered the exon body region, its downstream splice junction, downstream intron region and upstream intergenic region
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


cat ${TSSfolder}header.txt > ${codepath}a04_${func}_${j}.sh
cat>>${codepath}a04_${func}_${j}.sh<<EOF
coverageBed ${strand} ${bedA} ${i} ${bedB} ${inpath}tsscov.bed -split | tr -d "\r" > ${outpath}${j}/tsscov.cov
#coverageBed ${strand} ${bedA} ${i} ${bedB} ${inpath}tsscov_u.bed -split | tr -d "\r" > ${outpath}${j}/tsscov_u.cov
#coverageBed ${strand} ${bedA} ${i} ${bedB} ${inpath}tsscov_d.bed -split | tr -d "\r" > ${outpath}${j}/tsscov_d.cov
samtools view -H ${i} > ${outpath}${j}/junction
samtools view ${i} | awk -F"\t" '\$6~/N/' >> ${outpath}${j}/junction
samtools view -Sb ${outpath}${j}/junction > ${outpath}${j}/junction.b
rm ${outpath}${j}/junction
coverageBed ${strand} ${bedA} ${outpath}${j}/junction.b ${bedB} ${inpath}tsscov.bed -split | tr -d "\r" > ${outpath}${j}/tsscov_j.cov
genomeCoverageBed -bg -ibam ${i} -g ${genome} -split > ${outpath}${j}/coverage.bedgraph

EOF
j=$[ $j + 1 ]

done

#Output: tsscov.cov, tsscov_j.cov, tsscov_u.cov, tsscov_d.cov and coverage.bedgraph
#Reads counts covered the exon body region, its downstream splice junction, downstream intron region and upstream intergenic region, as well as the coverage used for Dapars






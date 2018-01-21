#!/bin/sh
############  Notice  #################
#Step6-2: To identify differential usage of tandem 5'UTR by using DaPars package
echo "To identify differential usage of tandem 5UTR by using DaPars model";

############  Parameters  #################
#Reading defined parameters
TSSfolder=$0;TSSfolder=${TSSfolder:0:${#TSSfolder}-12};
original_data=$1;outpath=$2;genome=$3;data1=$4;data2=$5;syspath=$6;Strand=$7;length1=$8;Splicing_diff_cutoff=$9;TypeData=${10};MultiProcessor=${11};

func=utr
codepath=${2}/code/
inpath=${2}/tmp/cov/
#Current_path record
Current_path=($(pwd))
mkdir $outpath

############  DaPars running step  #################
#Generating the configure file required parameters in DaPars
count1=0
j=1
group1=x
list=($(echo ${data1} | tr "," "\n"))
for i in ${list[*]}
do
group1=${group1}${inpath}${j}/coverage.bedgraph,
count1=$[ $count1 + 1 ]
j=$[ $j + 1 ]
done
group1=${group1%,*}
group1=${group1#*x}
count2=0
group2=x
list=($(echo ${data2} | tr "," "\n"))
for i in ${list[*]}
do
group2=${group2}${inpath}${j}/coverage.bedgraph,
count2=$[ $count2 + 1 ]
j=$[ $j + 1 ]
done
group2=${group2%,*}
group2=${group2#*x}


##Auto detection and change the relative path into absolute path for GTF
list=($(echo ${original_data} | tr "," "\n"))
for i in ${list[*]}
do
tmp=${i%/*};
if [ x${tmp} != "x" ]
then
gtffolder=($(cd ${tmp} ; pwd)/)
else
gtffolder=($(pwd)/)
fi
i=${gtffolder}${i##*/}
cd ${Current_path}

group3=${group3}${i},
done
original_data=${group3%,*}


#Generating the configure file required by DaPars model
cat>${2}/utr_configure.txt<<EOF
Annotated_3UTR=${2}/tmp/tsgtf/utr_tmp.bed
Group1_Tophat_aligned_Wig=${group1}
Group2_Tophat_aligned_Wig=${group2}
Output_directory=${2}
Output_result_file=utr_result
#Parameters
Num_least_in_group1=${count1}
Num_least_in_group2=${count2}
Coverage_cutoff=30
FDR_cutoff=0.05
PDUI_cutoff=0.3
Fold_change_cutoff=0.288

EOF

#Detecting the tandem 5'UTR events by using DaPars model
#Using split.R to segment the 5'UTR into two regions based on change-point model
#Count the reads coverage in each replicate; including the exon body region, its downstream splice junction, downstream intron region and upstream intergenic region
#The table.cov can be used for differential expression analysis of each first exon as the input of DESeq or edgeR
#Using format.R to estimate the the coverage and effective length of each first exon in one AFE event preparing for the estimation of PSI
#Using RMATS.sh to perform rMATS package. The null hypothesis and alternative hypothesis are proposed to detect whether the change of PSI for one first exon between the case and control groups exceeds the threshold    	
#Using FDR.R to calculate the adjusted P-Value

cat ${TSSfolder}header.txt > ${codepath}a06_${func}_test.sh
cat>>${codepath}a06_${func}_test.sh<<EOF
python ${syspath}dapars/DaPars_main.py ${2}/utr_configure.txt
Rscript ${syspath}split.R ${2}/utr_result_All_Prediction_Results.txt ${2}/
bash ${TSSfolder}cov_utr.sh ${original_data} ${outpath} ${genome} ${Strand}
bash ${TSSfolder}table_utr.sh ${original_data} ${outpath}
Rscript ${TSSfolder}format_utr.R ${outpath} ${length1}
bash ${TSSfolder}RMATS_utr.sh ${outpath} ${Splicing_diff_cutoff} ${TypeData} ${MultiProcessor} ${TSSfolder}
Rscript ${TSSfolder}FDR_utr.R ${outpath}

EOF




#Output: utr_result_All_Prediction_Results.txt, utr.annotation, utr_split.bed
#A standard output format as the DaPars documentation used for the prediction of tandem 5'UTR events
#utr.annotation: Annotation file of utr adding with fit value and changed point predicted by Dapars prediction model
#utr_split.bed: bed region of utr split into two regions
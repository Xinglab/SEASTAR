#!/bin/sh
############  Notice  #################
#Step6-2: To identify differential usage of tandem 5'UTR by using DaPars package
echo "To identify differential usage of tandem 5UTR by using DaPars model";

############  Parameters  #################
#Reading defined parameters
TSSfolder=$0;TSSfolder=${TSSfolder:0:${#TSSfolder}-12};
original_data=$1;outpath=$2;genome=$3;data1=$4;data2=$5;syspath=$6;

func=utr
codepath=${2}/code/
inpath=${2}/tmp/cov/

############  DaPars running step  #################
#Generating the configure file required parameters in DaPars
count1=0
j=1
list=($(echo ${data1} | tr "," "\n"))
for i in ${list[*]}
do
group1=${group1}${inpath}${j}/coverage.bedgraph,
count1=$[ $count1 + 1 ]
j=$[ $j + 1 ]
done
group1=${group1%,*}
count2=0
list=($(echo ${data2} | tr "," "\n"))
for i in ${list[*]}
do
group2=${group2}${inpath}${j}/coverage.bedgraph,
count2=$[ $count2 + 1 ]
j=$[ $j + 1 ]
done
group2=${group2%,*}

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
cat ${TSSfolder}header.txt > ${codepath}a06_${func}_test.sh
cat>>${codepath}a06_${func}_test.sh<<EOF
python ${syspath}dapars/DaPars_main.py ${2}/utr_configure.txt

EOF

#Output: utr_result_All_Prediction_Results.txt
#A standard output format as the DaPars documentation used for the prediction of tandem 5'UTR events

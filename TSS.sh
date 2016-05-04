#!/bin/bash

##Requied:
##	(1) Cufflinks >= 1.4.1
##  (2) samtools >= 0.1.19
##  (3) bedtools >= 2.15.0
##  (4) R >= 3.0.2
##  (5) rMATS >= 3.0.8 (Optional; for differential promoter detection)
##
##
##input example: scriptPath/TSS.sh -A inputdata -o outputfolder -g publicgtf -s sequence -d distance of TSSs cluster -c diff_cutoff(rMATS) -p multiprocessor -t modeltype(rMATS) -b batchprocess



#
#Parameters
TSSfolder=$0;
tmp=${TSSfolder:0:${#TSSfolder}-6};
if [ x${tmp} != "x" ]
then
TSSfolder=$(cd ${tmp} ; pwd)/
else
TSSfolder=$(pwd)/
fi


original_data="0";
output_folder="0";
public_gtf="0";
sequence="0";
distance="100"
Splicing_diff_cutoff=0.1;#Cutoff of >=splicing difference(eg. 0.1 for at least 10% splicing difference)
MultiProcessor=1;#Number of processor
TypeData="U";#Default is running unpaired data ('U' for unpaired data. 'P' for paired data.)
BatchType="U" #Default is generating multiple command files for batch processing. ("B" for batch processing. "U")
current_folder=$(pwd);check_code="replicate";

while getopts 'A:o:g:s:c:p:t:b:' optname
  do
    case $optname in
      A)
        original_data=$OPTARG;;
      o)
        output_folder=$OPTARG;;
      g)
        public_gtf=$OPTARG;;
      s)
        sequence=$OPTARG;;
      d)
        distance=$OPTARG;;
      c)
        Splicing_diff_cutoff=$OPTARG;;
      p)
        MultiProcessor=$OPTARG;;
      t)
        TypeData=$OPTARG;;
      b)
        BatchType=$OPTARG;;  
    esac
  done

if [ ${original_data} = "0" ];then
	echo "Error: No input data.";
else
	if [ ${output_folder} = "0" ];then
		echo "Error: No output folder.";
	else
		#select betch processing or not
    if [ "$BatchType" = "B" ];then
    	echo 'for batch processing manually';
    	rm -r ${output_folder} > /dev/null 2>&1
      mkdir ${output_folder} >/dev/null 2>&1
      mkdir ${output_folder}/tmp >/dev/null 2>&1
      output_folder=$(cd ${output_folder} ; pwd)

    	bash ${TSSfolder}gtf_batch.sh ${original_data} ${output_folder} ${public_gtf} ${MultiProcessor} 
    	bash ${TSSfolder}mrg_batch.sh ${original_data} ${output_folder} ${public_gtf} ${sequence} ${MultiProcessor} 
    	bash ${TSSfolder}tsgtf_batch.sh ${original_data} ${output_folder} ${output_folder}/tmp/mrg/merged_asm/merged.gtf ${MultiProcessor} 
    	
    	cat ${TSSfolder}header.txt > ${output_folder}/code/a03_Nrtss.sh
    	cat>>${output_folder}/code/a03_Nrtss.sh<<EOF
Rscript ${TSSfolder}FirstExon.R ${output_folder}/tmp/tsgtf
coverageBed -s -a ${output_folder}/tmp/tsgtf/nontss.bed -b ${output_folder}/tmp/tsgtf/ant.bed | tr -d "\r" > ${output_folder}/tmp/tsgtf/cmb
Rscript ${TSSfolder}TssMerge.R ${output_folder}/tmp/tsgtf ${distance}
coverageBed -s -a ${output_folder}/tmp/tsgtf/tssant.bed -b ${output_folder}/tmp/tsgtf/tssant.bed | tr -d "\r" > ${output_folder}/tmp/tsgtf/tssant.cov
Rscript ${TSSfolder}Nrtss.R ${output_folder}/tmp/tsgtf
EOF

    	bash ${TSSfolder}cov_batch.sh ${original_data} ${output_folder}
    	bash ${TSSfolder}table_batch.sh ${original_data} ${output_folder}
    	
    	cat ${TSSfolder}header.txt > ${output_folder}/code/a06_Test.sh
     	cat>>${output_folder}/code/a06_Test.sh<<EOF
Rscript ${TSSfolder}PresenceTest.R ${TSSfolder} ${output_folder}
 rm -r ${output_folder}/tmp/
EOF
   	
		else
    	echo 'for processing automated';
    	rm -r ${output_folder} > /dev/null 2>&1
      mkdir ${output_folder} >/dev/null 2>&1
      mkdir ${output_folder}/tmp >/dev/null 2>&1
      output_folder=$(cd ${output_folder} ; pwd)

    	bash ${TSSfolder}gtf.sh ${original_data} ${output_folder} ${public_gtf} ${MultiProcessor} 
    	bash ${TSSfolder}mrg.sh ${original_data} ${output_folder} ${public_gtf} ${sequence} ${MultiProcessor} 
    	bash ${TSSfolder}tsgtf.sh ${original_data} ${output_folder} ${output_folder}/tmp/mrg/merged_asm/merged.gtf ${MultiProcessor} 
    	
      Rscript ${TSSfolder}FirstExon.R ${output_folder}/tmp/tsgtf
    	coverageBed -s -a ${output_folder}/tmp/tsgtf/nontss.bed -b ${output_folder}/tmp/tsgtf/ant.bed | tr -d "\r" > ${output_folder}/tmp/tsgtf/cmb
    	Rscript ${TSSfolder}TssMerge.R ${output_folder}/tmp/tsgtf ${distance}
    	coverageBed -s -a ${output_folder}/tmp/tsgtf/tssant.bed -b ${output_folder}/tmp/tsgtf/tssant.bed | tr -d "\r" > ${output_folder}/tmp/tsgtf/tssant.cov
    	Rscript ${TSSfolder}Nrtss.R ${output_folder}/tmp/tsgtf

    	bash ${TSSfolder}cov.sh ${original_data} ${output_folder}
    	bash ${TSSfolder}table.sh ${original_data} ${output_folder}
    	
    	Rscript ${TSSfolder}PresenceTest.R ${TSSfolder} ${output_folder}
     rm -r ${output_folder}/tmp/
    	
		fi
	fi
fi

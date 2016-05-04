#!/bin/bash
#
#Parameters
MATSfolder=$0;
MATSfolder=${MATSfolder:0:${#MATSfolder}-8};
original_data="0";
output_folder="0";
Splicing_diff_cutoff=0.1;#Cutoff of >=splicing difference(eg. 0.1 for at least 10% splicing difference)
MultiProcessor=1;#Number of processor
TypeData="U";#Default is running unpaired data ('U' for unpaired data. 'P' for paired data.)
current_folder=$(pwd);check_code="replicate";

while getopts 'd:o:c:p:t:' optname
  do
    case $optname in
      d)
        original_data=$OPTARG;;
      o)
        output_folder=$OPTARG;;
      c)
        Splicing_diff_cutoff=$OPTARG;;
      p)
        MultiProcessor=$OPTARG;;
      t)
        TypeData=$OPTARG;;
    esac
  done

if [ "$original_data" = "0" ];then
	echo "Error: No input data.";
else
	if [ "$output_folder" = "0" ];then
		echo "Error: No output folder.";
	else
		#check the input file
		check_code=$(python "$MATSfolder"check_input.py "$original_data")
		if [ "$check_code" = "replicate" ];then
	             	if [ "$TypeData" = "P" ];then
				bash "$MATSfolder"rMATS_Paired.sh "$original_data" "$output_folder" "$Splicing_diff_cutoff" "$MultiProcessor" 
			else
				bash "$MATSfolder"rMATS_Unpaired.sh "$original_data" "$output_folder" "$Splicing_diff_cutoff" "$MultiProcessor" 
			fi
		else
			if [ "$check_code" = "pooled" ];then
				bash "$MATSfolder"MATS_LRT.sh "$original_data" "$output_folder" "$Splicing_diff_cutoff" "$MultiProcessor" 
			else
				echo "$check_code";
			fi
		fi
	fi
fi

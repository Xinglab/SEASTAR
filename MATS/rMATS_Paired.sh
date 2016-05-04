#!/bin/bash
#
#Parameters
MATSfolder=$0;MATSfolder=${MATSfolder:0:${#MATSfolder}-15};
original_data=$1;output_folder=$2;
Splicing_diff_cutoff=$3;#Cutoff of >=splicing difference(eg. 0.1 for at least 10% splicing difference)
MultiProcessor=$4;#1 for |psi1-psi2|>=cutoff; 2 for switch
current_folder=$(pwd)

#Create the output folder
rm -r "$2" > /dev/null 2>&1
mkdir "$2" >/dev/null 2>&1

#Paste the ID column out
awk '{print $1}' "$1" > "$2/ID.txt" 
awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' "$1" > "$2/Data.txt"

#Calculate Inclusion Level
python "$MATSfolder"inclusion_level.py "$2/Data.txt" "$2/Inc.txt" 

#Main
python "$MATSfolder"GLM_MS_paired.py "$1" "$2/" 50 84 "$4" "$3" > /dev/null 2>&1
python "$MATSfolder"FDR.py "$2/rMATS_Result_P.txt" "$2/rMATS_Result_FDR.txt"
paste "$2/rMATS_Result_FDR.txt" "$2/Inc.txt" > "$2/rMATS_Result.txt"
#awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$10"\t"$11"\t"$8"\t"$9}' "$2/rMATS_Result.tmp" > "$2/rMATS_Result.txt"

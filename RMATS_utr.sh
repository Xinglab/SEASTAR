#!/bin/sh
############  Notice  #################
#Using RMATS.sh to perform rMATS package. The null hypothesis and alternative hypothesis are proposed to detect whether the change of PSI for one first exon between the case and control groups exceeds the threshold    	
echo -e "
====================================================
[$(date +%R:%S)] To perform rMATS package\n";

############  Parameters  #################
#Reading defined parameters
TSSfolder=$0;TSSfolder=${TSSfolder:0:${#TSSfolder}-12};
outpath=$1;cutoff=$2;TypeData=$3;CPU=$4;syspath=$5;


func=RMATS_utr
inpath=${outpath}/tmp/RMATS_utr/
outpath=${outpath}/tmp/RMATS_utr/

############  rMATS running step  #################
#Loading all input files
cd $inpath
list=($(ls *.txt | tr "\n" " "))

#for each of input calculating the PSI and difference of PSI by rMATS model
for i in ${list[*]}
do
mkdir ${outpath}${i%%.*}_${cutoff}
cd ${outpath}${i%%.*}_${cutoff}

bash ${syspath}/MATS/rMATS.sh -d ${inpath}${i} -o ${outpath}${i%%.*}_${cutoff} -t ${TypeData} -c ${cutoff} -p ${CPU}

done

#Combining all the output files together
cd ${inpath}
list=($(ls -F |grep /$ | tr "\n" " "))
for i in ${list[*]}
do
  row=`cat ${outpath}rMATS_Result.txt | wc -l`
  if [ ${row} -eq 0 ]; then cat ${inpath}${i}rMATS_Result.txt | sed   -n   "1,   1p" > ${outpath}rMATS_Result.txt; fi

  linenum=`cat ${inpath}${i}rMATS_Result_FDR.txt | wc -l`
  if [ ${linenum} -ne 1 ]; then cat ${inpath}${i}rMATS_Result.txt | sed   -n   "2,   ${linenum}p" >> ${outpath}rMATS_Result.txt; fi
  
done

#rMATS_Result.txt: A standard output format as the rMATS documentation


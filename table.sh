#!/bin/sh
############  Notice  #################
#Step4-2: The table of reads counts for all samples; Each column represent one sample in the order of input and each raw represent one TSS
echo -e "
====================================================
[$(date +%R:%S)] Combine reads counts for all samples into the table\n";

############  Parameters  #################
#Reading defined parameters
TSSfolder=$0;TSSfolder=${TSSfolder:0:${#TSSfolder}-8};
original_data=$1;outpath=$2;

func=table
codepath=${2}/code/
outpath=${2}/tmp/${func}/
inpath=${2}/tmp/cov/
mkdir $outpath

############  Combining step  #################
#Loading all replicates in Group A and B
cd $inpath
list=($(ls -F | grep /$ | sort -n))

#Combining the counts of reads covered the exon body region
i=${list[0]}
cd ${inpath}
sort -k 4,4 ${inpath}${i}tsscov.cov | cut -f4 > ${outpath}table.cov
for i in ${list[*]}
do
sort -k 4,4 ${inpath}${i}tsscov.cov | cut -f7 | paste ${outpath}table.cov - > ${outpath}table.tmp
rm ${outpath}table.cov
mv ${outpath}table.tmp ${outpath}table.cov
done
cp ${outpath}table.cov ${2}/table.cov

##Combining the counts of reads covered the exon upstream intergenic region
#i=${list[0]}
#cd ${inpath}
#sort -k 4,4 ${inpath}${i}tsscov_u.cov | cut -f4 > ${outpath}u_table.cov
#for i in ${list[*]}
#do
#sort -k 4,4 ${inpath}${i}tsscov_u.cov | cut -f7 | paste ${outpath}u_table.cov - > ${outpath}u_table.tmp
#rm ${outpath}u_table.cov
#mv ${outpath}u_table.tmp ${outpath}u_table.cov
#done
#
##Combining the counts of reads covered the exon downstream intron region
#i=${list[0]}
#cd ${inpath}
#sort -k 4,4 ${inpath}${i}tsscov_d.cov | cut -f4 > ${outpath}d_table.cov
#for i in ${list[*]}
#do
#sort -k 4,4 ${inpath}${i}tsscov_d.cov | cut -f7 | paste ${outpath}d_table.cov - > ${outpath}d_table.tmp
#rm ${outpath}d_table.cov
#mv ${outpath}d_table.tmp ${outpath}d_table.cov
#done

#Combining the counts of reads covered the exon downstream splice junction
i=${list[0]}
cd ${inpath}
sort -k 4,4 ${inpath}${i}tsscov_j.cov | cut -f4 > ${outpath}j_table.cov
for i in ${list[*]}
do
sort -k 4,4 ${inpath}${i}tsscov_j.cov | cut -f7 | paste ${outpath}j_table.cov - > ${outpath}j_table.tmp
rm ${outpath}j_table.cov
mv ${outpath}j_table.tmp ${outpath}j_table.cov
done
cp ${outpath}j_table.cov ${2}/j_table.cov

############  FE presence test step  #################
# To identify the real first exons
# Using PresenceTest.R to identify the real first exons using the logistic model

Rscript ${TSSfolder}PresenceTest.R ${TSSfolder} ${2}


#Output: table.cov, u_table.cov, d_table.cov, j_table.cov and FilteredNrtss.annotation
#Counts table of reads covered the exon body region, downstream intron region, upstream intergenic region, downstream splice junction and real first exons



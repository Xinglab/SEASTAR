#!/bin/bash
#
TSSfolder=$0;
TSSfolder=${TSSfolder:0:${#TSSfolder}-11};
sequence=$1;
bash ${TSSfolder}TSS.sh -A ${TSSfolder}test/adipose2.bam -o ${TSSfolder}testresult -g ${TSSfolder}test/Homo_sapiens.Ensembl.GRCh37.72.chr7.gtf -s ${sequence} -p 4 -b U

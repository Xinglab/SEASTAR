#!/bin/bash
#
## used for simply test whether the pipeline works
echo "Used for simply test whether the pipeline works";
TSSfolder=$0;
TSSfolder=${TSSfolder:0:${#TSSfolder}-15};
sequence=$1;

bash ${TSSfolder}SEASTAR.sh -A ${TSSfolder}test/adipose1.bam,${TSSfolder}test/adipose2.bam -B ${TSSfolder}test/adrenal1.bam,${TSSfolder}test/adrenal2.bam -o ${TSSfolder}testresult -i ${TSSfolder}test/hg19.chrom.sizes -g ${TSSfolder}test/Homo_sapiens.Ensembl.GRCh37.72.chr7.gtf -s ${sequence} -c 0.1 -p 1 -b U


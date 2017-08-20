#!/bin/bash

################  Description  ################
##         Name: SEASTAR - Systematic Evaluation of Alternative STArt site in RNA
##      Authors: Zhiyi Qin, Bioinfomatics Laboratory, Department of Automation, Tsinghua University, qzy06@mails.tsinghua.edu.cn
##               Yi Xing, Department of Microbiology, Immunology, & Molecular Genetics, University of California, Los Angeles, yxing@ucla.edu
##      Created: May 5, 2014 for creation
## Last Revised: Apr 05, 2017 for fixed the bug in test command and gtf file with relative path
###############################################


################  Requirements  ################
##	(1) Cufflinks >= 1.4.1
##  (2) samtools >= 0.1.19
##  (3) bedtools >= 2.15.0
##  (4) R >= 3.0.2
##  (5) Python 2.6.x or Python 2.7.x; corresponding versions of NumPy, SciPy and rpy2
################################################


################  Quick Start  ################
## (1)   Installation: Run the install.sh (bash ./install.sh)
## (2)   Test example: bash [scriptPath]/test_SEASTAR.sh [BowtieIndex]
## (3)   De novo mode: bash [scriptPath]/SEASTAR.sh -A [inputdata of sample_1] -B [inputdata of sample_2] -o [outputfolder] -g [publicgtf] -i [genome_sizes] -s [sequence] -d [range of non-redundant TSSs] -p [multiprocessor] -c [diff_cutoff(for detecting AFE)] -t [modeltype(for detecting AFE)] -b [batchprocess] -S [strand-specific data]
## (4) Reference mode: bash [scriptPath]/SEASTAR.sh -A [inputdata of sample_1] -B [inputdata of sample_2] -o [outputfolder] -G [publicgtf] -i [genome_sizes] -s [sequence] -d [range of non-redundant TSSs] -p [multiprocessor] -c [diff_cutoff(for detecting AFE)] -t [modeltype(for detecting AFE)] -b [batchprocess] -S [strand-specific data]
###############################################


############  Default Parameters  #################
#Parameter of script path, used for both relative path and absolute path
TSSfolder=$0;
#Current_path record
Current_path=($(pwd))
#Auto detection and change the relative path into absolute path
tmp=${TSSfolder:0:${#TSSfolder}-10};
if [ x${tmp} != "x" ]
then
TSSfolder=($(cd ${tmp} ; pwd)/)
else
TSSfolder=($(pwd)/)
fi
cd ${Current_path}

#Quick Installation
cd ${TSSfolder}
bash ${TSSfolder}install.sh
cd ${Current_path}

#Input data; Mapping results for the sample in bam format
original_data="0";
original_data_1="0";
original_data_2="0";
#Output path
output_folder="0";
#Reference trancriptome; An annotation of genes and transcripts in GTF format
public_gtf="0";
#Reference trancriptome for Reference mode: skipping assembly step provided by Cufflinks
ref_gtf="0";
#Size of genome; The lengths of all chromosomes
genome_size="0";
#The fasta file of the bowtie indexes (fa files)
sequence="0";
#The distance among TSSs derived from same promoter region
distance="max"
#Cutoff of >=splicing difference(eg. 0.1 for at least 10% splicing difference)
Splicing_diff_cutoff=0.1;
#Number of processor
MultiProcessor=1;
#Default is running unpaired data ('U' for unpaired data. 'P' for paired data.)
TypeData="U";
#Default is running nonstrand specific data ('U' for nonstrand-specific data. 's' for strand-specific data.)
Strand="U";
#Default is generating multiple command files for batch processing. ("B" for batch processing. "U") 
#Batch processing is used for parallel computing on cluster. Non-batch processing is used for directly computing on PC.
BatchType="U" 
###################################################


###############  Reading user defined parameters  ###################
##Read in
while getopts 'A:B:o:g:G:i:s:d:c:p:t:S:b:' optname
  do
    case $optname in
      A)
        original_data_1=$OPTARG;;
      B)
        original_data_2=$OPTARG;;
      o)
        output_folder=$OPTARG;;
      g)
        public_gtf=$OPTARG;;
      G)
        ref_gtf=$OPTARG;;
      i)
        genome_size=$OPTARG;;
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
      S)
        Strand=$OPTARG;;
      b)
        BatchType=$OPTARG;;  
    esac
  done

#Showing the package version
echo -e "
====================================================
- SEASTAR [version 0.9.4]

Please see the webpage on Github about SEASTAR for more details: https://github.com/Xinglab/SEASTAR.git

===================================================="

#Showing help if there is not any input
tmp=$1
if [[ x${tmp} = "x" || ${tmp} = "-h" || ${tmp} = "-help" || ${tmp} = "--help" ]]
then
echo -e "Quick Start:

Test example: bash [scriptPath]/test_SEASTAR.sh [BowtieIndex]

De novo mode: bash [scriptPath]/SEASTAR.sh -A [inputdata of sample_1] -B [inputdata of sample_2] -o [outputfolder] -g [publicgtf] -i [genome_sizes] -s [sequence] -d [range of non-redundant TSSs] -p [multiprocessor] -c [diff_cutoff(for detecting AFE)] -t [modeltype(for detecting AFE)] -b [batchprocess] -S [strand-specific]

\nRequired parameters:

-A A_r1.bam[,A_r2.bam] The mapping results for the sample in bam format for the case group. Multiple alignments must be in a comma separated list (if using bam).

-B B_r1.bam[,B_r2.bam] The mapping results for the sample in bam format for the control group. Multiple alignments must be in a comma separated list (if using bam).

-g gtfFile Annotation of genes and transcripts in GTF format in De novo mode

-G gtfFile Annotation of genes and transcripts in GTF format skipping assembly step provided by Cufflinks

-i genomeSizes The lengths of all chromosomes. The format can be achieved from UCSC. More details can be found from the command genomeCoverageBed in bedtools.

-s bowtieIndexBase The fasta file of the bowtie indexes (fa files). The name should use hg19.fa instead of hg19. (Only used for assembly)

-o outDir The output directory for the generated results

\nOptional parameters:

-p <int> The number of processors to be used. The default value is 1.

-d <int> The distance among the TSSs derived from the same promoter region. The default is max (bps).

-c <float> The splicing difference cutoff. The cutoff is used in the null hypothesis test for differential splicing. The default is 0.1 for a 10% difference. The valid range is 0 ¡Ü cutoff < 1.

-t modeType Mode used in the MATS analysis. The options are 'U' for unpaired data and 'P' for paired data. Default is unpaired data.

-b batchProcess Type of batch processing to be done. The 'B' option generates several scripts to be used for batch processing on a parallel computing cluster. The 'U' option starts the analysis immediately. Default is 'U'.

-S The strand specific type of input data. The 's' option represents strand-specific data. The 'U' option represents non-strand specific data. Default is '`U`'.

\n====================================================
Please see the webpage on Github about SEASTAR for more details: https://github.com/Xinglab/SEASTAR.git \n
"
fi



#All of input data including Group A and Group B
original_data=${original_data_1},${original_data_2}
#Number of replicates in Group A
list1=($(echo ${original_data_1} | tr "," "\n"))
length1=${#list1[@]}
##########################################################


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
        2) bedA='-a'
        bedB='-b';; # "<"
    esac


################  Running the pipeline  ################
#Error without input data
if [ ${original_data_1} = "0" ];then
	echo "Error: No input data for case group.";
else
if [ ${original_data_2} = "0" ];then
	echo "Error: No input data for control group.";
else
#Error without output path
	if [ ${output_folder} = "0" ];then
		echo "Error: No output folder.";
	else
#Select betch processing; Batch processing is used for parallel computing on cluster.
    if [ "$BatchType" = "B" ];then
    	echo 'for batch processing manually';
#Parameter of output path, used for both relative path and absolute path
      rm -r ${output_folder} > /dev/null 2>&1
      mkdir ${output_folder} >/dev/null 2>&1
      mkdir ${output_folder}/tmp >/dev/null 2>&1
      mkdir ${output_folder}/code >/dev/null 2>&1

      output_folder=($(cd ${output_folder} ; pwd))
      cd ${Current_path}

      if [ "$ref_gtf" = "0" ];then
            
#Auto detection and change the relative path into absolute path for GTF
tmp=${public_gtf%/*};
if [ x${tmp} != "x" ]
then
gtffolder=($(cd ${tmp} ; pwd)/)
else
gtffolder=($(pwd)/)
fi
public_gtf=${gtffolder}${public_gtf##*/}
cd ${Current_path}


#Step1: Assembly of novel transcripts by Reference Annotation Based Transcript (RABT) method using Cufflinks in each sample based on RNA-seq reads and reference annotation
    	  bash ${TSSfolder}gtf_batch.sh ${original_data} ${output_folder} ${public_gtf} ${MultiProcessor} 
    	
#Step2-1: Merge each annotation file (GTF file) generated from the assembly step in each sample together using Cuffmerge
    	  bash ${TSSfolder}mrg_batch.sh ${original_data} ${output_folder} ${public_gtf} ${sequence} ${MultiProcessor} 

#Step2-2: Reformat the annotation of transcriptome and an interface for more FPKM analysis
    	  bash ${TSSfolder}tsgtf_batch.sh ${original_data} ${output_folder} ${output_folder}/tmp/mrg/merged_asm/merged.gtf ${MultiProcessor} 
    	else
    	
#Auto detection and change the relative path into absolute path for GTF
tmp=${ref_gtf%/*};
if [ x${tmp} != "x" ]
then
gtffolder=($(cd ${tmp} ; pwd)/)
else
gtffolder=($(pwd)/)
fi
ref_gtf=${gtffolder}${ref_gtf##*/}
cd ${Current_path}


#Step2-2: Reformat the annotation of transcriptome and an interface for more FPKM analysis
    	  bash ${TSSfolder}tsgtf_batch.sh ${original_data} ${output_folder} ${ref_gtf} ${MultiProcessor} 
    	fi
    	
#Step3: If the first exons are overlapped with each other, we merge these first exons together, so called non-redundant first exon
#Write the commands into the script file
    	cat ${TSSfolder}header.txt > ${output_folder}/code/a03_Nrtss.sh
    	cat>>${output_folder}/code/a03_Nrtss.sh<<EOF
Rscript ${TSSfolder}FirstExon.R ${output_folder}/tmp/tsgtf
coverageBed -s ${bedA} ${output_folder}/tmp/tsgtf/nontss.bed ${bedB} ${output_folder}/tmp/tsgtf/ant.bed | tr -d "\r" > ${output_folder}/tmp/tsgtf/cmb
Rscript ${TSSfolder}TssMerge.R ${output_folder}/tmp/tsgtf ${distance}
coverageBed -s ${bedA} ${output_folder}/tmp/tsgtf/tssant.bed ${bedB} ${output_folder}/tmp/tsgtf/tssant.bed | tr -d "\r" > ${output_folder}/tmp/tsgtf/tssant.cov
Rscript ${TSSfolder}Nrtss.R ${output_folder}/tmp/tsgtf
EOF

#Step4: Count the reads coverage in each replicate; including the exon body region, its downstream splice junction, downstream intron region and upstream intergenic region, as well as the coverage used for Dapars
    	bash ${TSSfolder}cov_batch.sh ${original_data} ${output_folder} ${genome_size} ${Strand}
    	
#Step5: The table for all samples; Each column represent one sample in the order of input and each raw represent one TSS
#       To identify the real first exons
#The table.cov can be used for differential expression analysis of each first exon as the input of DESeq or edgeR
#Using PresenceTest.R to identify the real first exons using the logistic model
    	bash ${TSSfolder}table_batch.sh ${original_data} ${output_folder}
    	
#Step6-1: To test whether the differential usage of AFE is significant by using rMATS package
#Using format.R to estimate the the coverage and effective length of each first exon in one AFE event preparing for the estimation of PSI
#Using RMATS.sh to perform rMATS package. The null hypothesis and alternative hypothesis are proposed to detect whether the change of PSI for one first exon between the case and control groups exceeds the threshold    	
#Using FDR.R to calculate the adjusted P-Value
    	cat ${TSSfolder}header.txt > ${output_folder}/code/a06_Test.sh
     	cat>>${output_folder}/code/a06_Test.sh<<EOF
      mkdir ${output_folder}/tmp/RMATS >/dev/null 2>&1
Rscript ${TSSfolder}format.R ${output_folder} ${length1}
    	bash ${TSSfolder}RMATS.sh ${output_folder} ${Splicing_diff_cutoff} ${TypeData} ${MultiProcessor} ${TSSfolder}
Rscript ${TSSfolder}FDR.R ${output_folder}
    	rm -r ${output_folder}/tmp/RMATS/ >/dev/null 2>&1
EOF

#Step6-2: To identify differential usage of tandem 5'UTR by using DaPars package
    	mkdir ${output_folder}/tmp/RMATS_utr >/dev/null 2>&1
    	bash ${TSSfolder}utr_batch.sh ${original_data} ${output_folder} ${genome_size} ${original_data_1} ${original_data_2} ${TSSfolder} ${Strand} ${length1} ${Splicing_diff_cutoff} ${TypeData} ${MultiProcessor}

    	
    	
    	
    	
    	
   	# rm -r ${output_folder}/tmp/


#Select non-betch processing; Non-batch processing is used for directly computing on PC.
		else
    	echo 'for processing automated';
#Parameter of output path, used for both relative path and absolute path
      set -x

    	rm -r ${output_folder} > /dev/null 2>&1
      mkdir ${output_folder} >/dev/null 2>&1
      mkdir ${output_folder}/tmp >/dev/null 2>&1
      mkdir ${output_folder}/code >/dev/null 2>&1

      output_folder=($(cd ${output_folder} ; pwd))
      cd ${Current_path}


#Record and show each command on the screen
      if [ "$ref_gtf" = "0" ];then

#Auto detection and change the relative path into absolute path for GTF
tmp=${public_gtf%/*};
if [ x${tmp} != "x" ]
then
gtffolder=($(cd ${tmp} ; pwd)/)
else
gtffolder=($(pwd)/)
fi
public_gtf=${gtffolder}${public_gtf##*/}
cd ${Current_path}


#Step1: Assembly of novel transcripts by Reference Annotation Based Transcript (RABT) method using Cufflinks in each sample based on RNA-seq reads and reference annotation
    	  bash ${TSSfolder}gtf.sh ${original_data} ${output_folder} ${public_gtf} ${MultiProcessor} 
    	
#Step2-1: Merge each annotation file (GTF file) generated from the assembly step in each sample together using Cuffmerge
    	  bash ${TSSfolder}mrg.sh ${original_data} ${output_folder} ${public_gtf} ${sequence} ${MultiProcessor} 
    	
#Step2-2: Reformat the annotation of transcriptome and an interface for more FPKM analysis
    	  bash ${TSSfolder}tsgtf.sh ${original_data} ${output_folder} ${output_folder}/tmp/mrg/merged_asm/merged.gtf ${MultiProcessor} 
    	else

#Auto detection and change the relative path into absolute path for GTF
tmp=${ref_gtf%/*};
if [ x${tmp} != "x" ]
then
gtffolder=($(cd ${tmp} ; pwd)/)
else
gtffolder=($(pwd)/)
fi
ref_gtf=${gtffolder}${ref_gtf##*/}
cd ${Current_path}


#Step2-2: Reformat the annotation of transcriptome and an interface for more FPKM analysis
    	  bash ${TSSfolder}tsgtf.sh ${original_data} ${output_folder} ${ref_gtf} ${MultiProcessor} 
      fi
    	
#Step3: If the first exons are overlapped with each other, we merge these first exons together, so called non-redundant first exon
      Rscript ${TSSfolder}FirstExon.R ${output_folder}/tmp/tsgtf
    	coverageBed -s ${bedA} ${output_folder}/tmp/tsgtf/nontss.bed ${bedB} ${output_folder}/tmp/tsgtf/ant.bed | tr -d "\r" > ${output_folder}/tmp/tsgtf/cmb
    	Rscript ${TSSfolder}TssMerge.R ${output_folder}/tmp/tsgtf ${distance}
    	coverageBed -s ${bedA} ${output_folder}/tmp/tsgtf/tssant.bed ${bedB} ${output_folder}/tmp/tsgtf/tssant.bed | tr -d "\r" > ${output_folder}/tmp/tsgtf/tssant.cov
    	Rscript ${TSSfolder}Nrtss.R ${output_folder}/tmp/tsgtf
    	
#Step4: Count the reads coverage in each replicate; including the exon body region, its downstream splice junction, downstream intron region and upstream intergenic region
    	bash ${TSSfolder}cov.sh ${original_data} ${output_folder} ${genome_size}  ${Strand}
    	
#Step5: The table for all samples; Each column represent one sample in the order of input and each raw represent one TSS
#       To identify the real first exons
#The table.cov can be used for differential expression analysis of each first exon as the input of DESeq or edgeR
#Using PresenceTest.R to identify the real first exons using the logistic model
    	bash ${TSSfolder}table.sh ${original_data} ${output_folder}

#Step6-1: To test whether the differential usage of AFE is significant by using rMATS package
#Using format.R to estimate the the coverage and effective length of each first exon in one AFE event preparing for the estimation of PSI
#Using RMATS.sh to perform rMATS package. The null hypothesis and alternative hypothesis are proposed to detect whether the change of PSI for one first exon between the case and control groups exceeds the threshold    	
#Using FDR.R to calculate the adjusted P-Value
    	mkdir ${output_folder}/tmp/RMATS >/dev/null 2>&1
    	Rscript ${TSSfolder}format.R ${output_folder} ${length1}
    	bash ${TSSfolder}RMATS.sh ${output_folder} ${Splicing_diff_cutoff} ${TypeData} ${MultiProcessor} ${TSSfolder}
    	Rscript ${TSSfolder}FDR.R ${output_folder}
    	rm -r ${output_folder}/tmp/RMATS/ >/dev/null 2>&1
    	
#Step6-2: To identify differential usage of tandem 5'UTR by using DaPars package
    	mkdir ${output_folder}/tmp/RMATS_utr >/dev/null 2>&1
    	bash ${TSSfolder}utr.sh ${original_data} ${output_folder} ${genome_size} ${original_data_1} ${original_data_2} ${TSSfolder} ${Strand} ${length1} ${Splicing_diff_cutoff} ${TypeData} ${MultiProcessor}


    	
    	
    	# rm -r ${output_folder}/tmp/
		fi
	fi
fi
fi
#######################################
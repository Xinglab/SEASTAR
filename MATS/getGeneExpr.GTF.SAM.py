#This script generate the gene expression based on the SAM files
import re,os,sys,commands

#Use dictionary to speed up the match
chunk=1000;

#build annotation dictionary
#read the GTF files
ifile=open(sys.argv[1]);
ilines=ifile.readlines();
Loc2Gene={};Gene_dic={};Loc_Gene_dic={};
for i in ilines:
	element=re.findall('[^\t\n]+',i);
	if element[2]=='exon':
		#gene=re.findall('ENSG\d+',i);
		gene=re.findall('gene_id[^\t]+',i);
		if len(gene)==0: #No gene name
			continue; 
		gene=re.findall('[^ ";]+',gene[0]);#remove unused marks in the gene names
		if len(gene)<2: #No gene name
			continue; 
		gene=gene[1];
		Gene_dic[gene]=0;
		start=int(element[3]);
		stop=int(element[4]);
		group_start=start/chunk;group_stop=stop/chunk;
		#To speed up, the same exon of the same gene only register once
		#Partically overlap exons won't cause problem here because counts only register once for a gene
		if not((element[0]+'_'+str(start)+'_'+str(stop)+'_'+str(gene)) in Loc_Gene_dic):
			Loc_Gene_dic[element[0]+'_'+str(start)+'_'+str(stop)+'_'+str(gene)]=0;
			if not(element[0] in Loc2Gene):
				Loc2Gene[element[0]]={};
			for this_group in range(group_start,group_stop+1):
				if not(this_group in Loc2Gene[element[0]]):
					Loc2Gene[element[0]][this_group]=[];
				Loc2Gene[element[0]][this_group].append([start,stop,gene]);

gene_list=Gene_dic.keys();gene_list_count={};
ifile.close();
#print(Loc2Gene);

#Read in the Total Number of Counts Provided by Upstream
count_list=re.findall('[^,]+',sys.argv[3]);
for i in range(len(count_list)):
	count_list[i]=int(count_list[i]);

#Read in the SAM files and add counts to genes based on the SAM files
ifile_list=re.findall('[^,]+',sys.argv[2]);
for i in range(len(ifile_list)):
	gene_count={};ofile=open(ifile_list[i]+'.GeneExpr.tmp','w');
	ifile=open(ifile_list[i]);
	for iline in ifile:
                if iline[0]=='@': ### it is the header..
                	continue;
		element=re.findall('[^\t\n]+',iline);
		chr=element[2];
		if chr[0:3]!='chr':
			chr = 'chr'+chr;
		start=int(element[3])-1; #convert to 0 base
		group=start/chunk; 
		mString=element[5]; # mapping string, 50M or aMbNcM format
		if 'D' in mString or 'I' in mString or 'S' in mString or 'H' in mString or 'P' in mString or 'X' in mString or '=' in mString:
			continue;
		# check to see if the line is either exonic read or junction read
		split_mString = mString.split('M');
		if len(split_mString) == 2: #exonic read
			rL=int(split_mString[0]);
			stop=start+rL;
			#Find the location of the exonic read
			if chr in Loc2Gene:
				if group in Loc2Gene[chr]:
					gene_of_read={};#make sure only add 1 to each gene
					for loc in Loc2Gene[chr][group]:
						if ((start>=loc[0]) & (stop <=loc[1])): #compare with the indxed exon location
							if not(loc[2] in gene_of_read):
								gene_of_read[loc[2]]=0;
								if not(loc[2] in gene_count):
									gene_count[loc[2]]=0;
								gene_count[loc[2]]+=1;
	#finish locating SAM reads in the genes, start to write temp files in case the script is aborted
	for j in gene_list:
		if not(j in gene_list_count):
			gene_list_count[j]=[];
		if j in gene_count:
			#ofile.write(j+'\t'+str(round(float(gene_count[j])/count_list[i],3))+'\n');
			#gene_list_count[j].append(str(round(float(gene_count[j])/count_list[i],3)));
			ofile.write(j+'\t'+str(1000000.0*float(gene_count[j])/count_list[i])+'\n');
			gene_list_count[j].append(str(1000000.0*float(gene_count[j])/count_list[i]));
		else:
			ofile.write(j+'\t0\n');
			gene_list_count[j].append('0');
	ofile.close();
	ifile.close();

#Generate the Gene Expression Summary File
ofile=open(sys.argv[4],'w');ofile.write('GeneID');
for i in ifile_list:
	ofile.write('\t'+i);
ofile.write('\n')
for i in gene_list:
	ofile.write(i);
	for j in gene_list_count[i]:
		ofile.write('\t'+j);
	ofile.write('\n');
ofile.close();


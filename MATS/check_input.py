#this script scans the input file for errors

import re,os,sys

res='replicate';

ifile=open(sys.argv[1]);
ifile.readline();
ilines=ifile.readlines();

len_group_1=-1;len_group_2=-1;
for i in ilines:
	element=re.findall('[^\t]+',i);
	inc1=re.findall('[^,]+',element[1]);skp1=re.findall('[^,]+',element[2]);
	inc2=re.findall('[^,]+',element[3]);skp2=re.findall('[^,]+',element[4]);

	if (len(inc1)+len(inc2))==2:
		res='pooled';

	#check length of inclusion / skipping
	if len(inc1)!=len(skp1):
		res='Error: different number of inclusion and skipping counts.\n'+str(i);break;
	if len(inc2)!=len(skp2):
		res='Error: different number of inclusion and skipping counts.\n'+str(i);break;

	#check the length of inclusion / skipping that is the same for all the exons
	if len_group_1==-1:
		len_group_1=len(inc1);
	if len_group_2==-1:
		len_group_2=len(inc2);
	if len_group_1!=len(inc1):
		res='Error: number of inclusion and skipping counts are not the same for different exons.\n'+str(i);break;
	if len_group_2!=len(inc2):
		res='Error: number of inclusion and skipping counts are not the same for different exons.\n'+str(i);break;
	
print(res);




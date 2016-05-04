
#####   FirstExon.R
##Collecting all first exons from transcripts annotation
##Last command: Linux, cut -d" " -f1-6 transcripts.gtf | sed 's/ /\t/g' |sed 's/;//g' |sed 's/"//g'| cut -f1,3,4,5,7,10,12,14 >transcripts_cut8.gtf
argv <- commandArgs(TRUE)
path <- paste(argv[1],"",sep='')
#path <-"E:\\Placenta\\Data"
gtf <- as.matrix(read.table(paste(path,"/transcripts_cut8.gtf",sep=''),header = FALSE, sep = "\t",quote = ""))
iso <- as.matrix(read.table(paste(path,"/isoforms.fpkm_tracking",sep=''),header = TRUE, sep = "\t",quote = ""))
iso <- iso[,1:8]
tmp <- rep(0,nrow(gtf))

for (i in nrow(gtf):1) {
 if (gtf[i,5]=="-") {
  if (i==nrow(gtf) || (i<nrow(gtf) && gtf[i+1,2]=="transcript")) {   tmp[i] <- 1  }
 } else {
  if (as.numeric(gtf[i,8])==1) {   tmp[i] <- 1  }
 }
}
ant <- gtf[tmp==1,]  ##all first exons of each transcript
nontss <- gtf[((tmp==0)+(gtf[,2]=="exon"))==2,]  ##all other exons
ant <- ant[order(ant[,7]),]
nontss <- nontss[order(nontss[,7]),]
iso <- iso[order(iso[,1]),]  ##all isoforms information

#ant[,6]!=iso[,4]
ant <- cbind(iso,ant[,c(1,3,4,5,8,8,8)])
colnames(ant) <- c(colnames(ant)[1:8],"tss_chr","tss_ss","tss_es","strand","exon_num","in/out","over/non")

## annotate out=1,in=0
ant <- ant[order(ant[,4],as.numeric(ant[,10])),]
for (i in 1:nrow(ant)) {
 if (ant[i,12]=="+"){
  if (i==1 || (i!=1 && ant[i,4]!=ant[i-1,4])){ first <- as.numeric(ant[i,10]) }
  if (as.numeric(ant[i,10])==first) {   ant[i,14] <- 1
  } else {   ant[i,14] <- 0  }
 }
}
ant <- ant[order(ant[,4],-as.numeric(ant[,11])),]
for (i in 1:nrow(ant)) {
 if (ant[i,12]=="-"){
  if (i==1 || (i!=1 && ant[i,4]!=ant[i-1,4])){ first <- as.numeric(ant[i,11]) }
  if (as.numeric(ant[i,11])==first) {   ant[i,14] <- 1
  } else {   ant[i,14] <- 0  }
 }
}
ant <- ant[order(ant[,1]),]

colnames(nontss) <- c("chr","type","ss","es","strand","gene_id","tracking_id","exon_num")
write.table(ant, paste(path,"/ant.annotation",sep=''),sep="\t", row.names=FALSE, col.name=TRUE, quote = FALSE)

#ant is all first exons based on transcripts record, 
ant <- gsub(pattern=" ", replacement="", ant)
nontss <- gsub(pattern=" ", replacement="", nontss)
write.table(ant, paste(path,"/ant.annotation",sep=''),sep="\t", row.names=FALSE, col.name=TRUE, quote = FALSE)
write.table(ant[,c(9:11,1,14,12)], paste(path,"/ant.bed",sep=''),sep="\t", row.names=FALSE, col.name=FALSE, quote = FALSE)
write.table(nontss[,c(1,3,4,7,8,5)], paste(path,"/nontss.bed",sep=''),sep="\t", row.names=FALSE, col.name=FALSE, quote = FALSE)

#ant.annotation: All the first exons abstracted from transcript annotation
#ant.bed and nontss.bed are used for detecting whether the first exons are overlapped with other non-first exons
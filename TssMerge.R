#####  TssMerge.R
##If the first exons are overlapped with each other, we merge these first exons together, so called non-redundant first exon
##User given threshold "distss" for merging overlapped first exons together (if "distss"> distance between TSSs for the overlapped first exons)
argv <- commandArgs(TRUE)
path <- argv[1]
distss <- argv[2]
#path <-"E:\\Placenta\\Data"
#distss <- 100
#### in each transcript the first exon' overlapped fraction with non-tss exons ####
## need ant.bed, nontss.bed


frc <- as.matrix(read.table(paste(path,"/cmb",sep=''),header = FALSE, sep = "\t",quote = ""))
frc <- frc[order(frc[,4]),]
ant <- as.matrix(read.table(paste(path,"/ant.annotation",sep=''),header = TRUE, sep = "\t",quote = ""))

ant <- cbind(ant[,1:14],frc[,7:10])
colnames(ant) <- c(colnames(ant)[1:14],"over_cnt","over_len","tss_len","over_frc")
write.table(ant, paste(path,"/ant.annotation",sep=''),sep="\t", row.names=FALSE, col.name=TRUE, quote = FALSE)
####first exon of transcrpts, ant.annotation(ant)



#### each non-redundant tss' overlapped fraction with non-tss exons, using the max result ####
#### non-redundant tss, need ant.annotation(ant)
#path <-"E:\\Placenta\\Data\\"
ant <- read.delim(paste(path,"/ant.annotation",sep=''), header=TRUE)

#### annotate TSS_ID by user given margin (distance among TSS clusters)
ant <- ant[order(ant[,"tss_chr"],ant[,"strand"],ant[,"tss_ss"]),]
ss <- ant[,"tss_ss"]
ss[ant[,"strand"]=="-"] <- -1*as.numeric(as.vector(ant[ant[,"strand"]=="-","tss_es"]))
ant <- cbind(ant,ss)
ant <- ant[order(ant[,"tss_chr"],ant[,"strand"],ant[,"ss"]),]
ant[,"tss_id"] <- "TSS"
tmp <- ant

#identify the shared tss_id among overlapped first exons, positive strand
ant <- tmp[tmp[,"strand"]!="-",]
j <- 1
ant[1,"tss_id"] <- paste("TSS",j,sep="")
mergedr <- as.numeric(ant[1,"tss_es"])
for (i in 2:nrow(ant)) {
 if ((as.vector(ant[i-1,"tss_chr"])==as.vector(ant[i,"tss_chr"])) && (as.vector(ant[i-1,"strand"])==as.vector(ant[i,"strand"])) && (mergedr >= as.numeric(ant[i,"tss_ss"])) && (distss=="max" || abs(as.numeric(ant[i-1,"ss"])-as.numeric(ant[i,"ss"]))<=as.numeric(distss))) {
  ant[i,"tss_id"] <- as.vector(ant[i-1,"tss_id"])
  
  #demand for the shorest end site (no varied 5'ss region; min)
  mergedr <- min(mergedr,as.numeric(ant[i,"tss_es"]))

 } else {
  j <- j+1
  ant[i,"tss_id"] <- paste("TSS",j,sep="")
  mergedr <- as.numeric(ant[i,"tss_es"])
 }
}
tmpp <- ant

#identify the shared tss_id among overlapped first exons, negative strand
ant <- tmp[tmp[,"strand"]=="-",]
j <- j+1
ant[1,"tss_id"] <- paste("TSS",j,sep="")
mergedr <- as.numeric(ant[1,"tss_ss"])
for (i in 2:nrow(ant)) {
 if ((as.vector(ant[i-1,"tss_chr"])==as.vector(ant[i,"tss_chr"])) && (as.vector(ant[i-1,"strand"])==as.vector(ant[i,"strand"])) && (mergedr <= as.numeric(ant[i,"tss_es"])) && (distss=="max" || abs(as.numeric(ant[i-1,"ss"])-as.numeric(ant[i,"ss"]))<=as.numeric(distss))) {
  ant[i,"tss_id"] <- as.vector(ant[i-1,"tss_id"])
  
  #demand for the shorest end site (no varied 5'ss region; max)
  mergedr <- max(mergedr,as.numeric(ant[i,"tss_ss"]))

 } else {
  j <- j+1
  ant[i,"tss_id"] <- paste("TSS",j,sep="")
  mergedr <- as.numeric(ant[i,"tss_ss"])
 }
}
ant <- rbind(tmpp,ant)

#generate the non-redundant first exons' ss and es (without the varied 5'ss region)


ant <- ant[order(ant[,"tss_id"],-ant[,"length"]),]
tss <- ant[!duplicated(ant$tss_id),c(6,2,4,5,7:12,14:18)]
tss[tss$strand!="-",8:9] <- aggregate(ant[ant[,12]!="-",10:11], list(ant[ant[,12]!="-","tss_id"]), min)[,2:3]
tss[tss$strand=="-",8:9] <- aggregate(ant[ant[,12]=="-",10:11], list(ant[ant[,12]=="-","tss_id"]), max)[,2:3]
ant <- ant[order(ant[,"tss_id"],-ant[,"over_frc"]),]
tss[,12:15] <- ant[!duplicated(ant$tss_id),15:18]
tss[,"tss_len"] <- tss[,"tss_es"] - tss[,"tss_ss"]
tss <- tss[tss[,"tss_len"]!="0",]


write.table(tss, paste(path,"/tss_short.annotation",sep=""), row.names=FALSE, col.name=TRUE, quote = FALSE, sep="\t")
####tss_id(old nrtss), tss_short.annotation(tss)

#### the non-redundant tss' overlapped fraction with themself ####
## need to generate tssant.bed
tssant <- as.matrix(read.table(paste(path,"/tss_short.annotation",sep=''),header = TRUE, sep = "\t",quote = ""))
tssant <- gsub(pattern=" ", replacement="", tssant)
write.table(tssant[,c(7,8,9,1,11,10)], paste(path,"/tssant.bed",sep=''),sep="\t", row.names=FALSE, col.name=FALSE, quote = FALSE)

#tss_short.annotation: the non-redundant first exons annotation
#tssant.bed: used for calculate the overlapped fractions between non-redundant first exons
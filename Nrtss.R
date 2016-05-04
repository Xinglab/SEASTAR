#####  Nrtss.R
##Generating the final annotation of non-redundant first exons

argv <- commandArgs(TRUE)
path <- argv[1]


tssant_frc <- as.matrix(read.table(paste(path,"/tssant.cov",sep=''),header = FALSE, sep = "\t",quote = ""))
tssant <- as.matrix(read.table(paste(path,"/tss_short.annotation",sep=''),header = TRUE, sep = "\t",quote = ""))
tssant_frc <- tssant_frc[order(tssant_frc[,4]),]
tssant <- tssant[order(tssant[,1]),]
colnames(tssant)[9] <- "tss_es"
tssant <- cbind(tssant,tssant_frc[,7])
colnames(tssant)[16] <- "tss_ovcnt"
tssant <- gsub(pattern=" ", replacement="", tssant)

write.table(tssant, paste(path,"/tss.annotation",sep=''),sep="\t", row.names=FALSE, col.name=TRUE, quote = FALSE)
####tss_id(old nrtss) with overlap_exon and overlap_nrtss, tss.annotation(tssant)

##generate the non-redundant tss genome annotation, tssant[,15]==0&&tssant[,16]==1
tssant <- as.matrix(read.table(paste(path,"/tss.annotation",sep=''),header = TRUE, sep = "\t",quote = ""))
tssant <- gsub(pattern=" ", replacement="", tssant)
write.table(tssant[((as.numeric(tssant[,15])==0)+(as.numeric(tssant[,16])==1))==2,c(7,8,9,1,11,10)], paste(path,"/tsscov.bed",sep=''),sep="\t", row.names=FALSE, col.name=FALSE, quote = FALSE)
write.table(tssant[((as.numeric(tssant[,15])==0)+(as.numeric(tssant[,16])==1))==2,], paste(path,"/nrtss.annotation",sep=''),sep="\t", row.names=FALSE, col.name=TRUE, quote = FALSE)
####non-overlapped tss(final nrtss), nrtss.annotation(tssant after filtering)


############################  generate up down region for test  #######################################
##generate the non-redundant tss genome annotation, tssant[,15]==0&&tssant[,16]==1
#upstream and downstream
tssant <- as.matrix(read.table(paste(path,"/tss.annotation",sep=''),header = TRUE, sep = "\t",quote = ""))
tssant[,"tss_es"] <- as.numeric(tssant[,"tss_es"])
tssant[,"tss_ss"] <- as.numeric(tssant[,"tss_ss"])
up <- tssant
up[,"tss_ss"] <- 2*as.numeric(up[,"tss_ss"])-as.numeric(up[,"tss_es"])-1
up[as.numeric(up[,"tss_ss"])<=0,"tss_ss"] <- 1
up[,"tss_es"] <- as.numeric(tssant[,"tss_ss"])-1
up[as.numeric(up[,"tss_es"])<=0,"tss_es"] <- 1
#tssant[as.numeric(up[,"tss_ss"])>=as.numeric(up[,"tss_es"]),]

down <- tssant
down[,"tss_es"] <- 2*as.numeric(down[,"tss_es"])-as.numeric(down[,"tss_ss"])+1
down[,"tss_ss"] <- as.numeric(tssant[,"tss_es"])+1
down[as.numeric(down[,"tss_ss"])<=0,"tss_ss"] <- 1
down[as.numeric(down[,"tss_es"])<=0,"tss_es"] <- 1
write.table(tssant[((as.numeric(tssant[,15])==0)+(as.numeric(tssant[,16])==1))==2,c(7,8,9,1,11,10)], paste(path,"/tsscov.bed",sep=''),sep="\t", row.names=FALSE, col.name=FALSE, quote = FALSE)
write.table(up[((as.numeric(tssant[,15])==0)+(as.numeric(tssant[,16])==1))==2,c(7,8,9,1,11,10)], paste(path,"/tsscov_u.bed",sep=''),sep="\t", row.names=FALSE, col.name=FALSE, quote = FALSE)
write.table(down[((as.numeric(tssant[,15])==0)+(as.numeric(tssant[,16])==1))==2,c(7,8,9,1,11,10)], paste(path,"/tsscov_d.bed",sep=''),sep="\t", row.names=FALSE, col.name=FALSE, quote = FALSE)

##tsscov.bed, tsscov_u.bed and tsscov_d.bed: the genomic regions for counting reads used for identifing the real first exons
##### PresenceTest.R
#Using PresenceTest.R to identify the real first exons using the logistic model

argv <- commandArgs(TRUE)
syspath <- argv[1]
path <- argv[2]
#path <-"E:\\Placenta\\Data"
#######################  exon presence test  #####################################
#### check assembled TSS of RNA-seq with CAGE

##check TSS after filter after overlap (as just count tss after overlap)
##tssant tss.annotation is all tss before overlap check and filter check
HBM <- read.delim(paste(path,"/tmp/table/table.cov",sep=''), header = FALSE)
PLC <- read.delim(paste(path,"/tmp/table/table.cov",sep=''), header = FALSE)
table <- cbind(PLC[,-1],HBM[,-1])
rownames(table) <- PLC[,1]
colnames(table) <- 1:ncol(table)
table <- table[order(as.vector(PLC[,1])),]
e <- apply(table,1,max)
HBM <- read.delim(paste(path,"/tmp/table/j_table.cov",sep=''), header = FALSE)
PLC <- read.delim(paste(path,"/tmp/table/j_table.cov",sep=''), header = FALSE)
table <- cbind(PLC[,-1],HBM[,-1])
rownames(table) <- PLC[,1]
colnames(table) <- 1:ncol(table)
table <- table[order(as.vector(PLC[,1])),]
j <- apply(table,1,max)
nrtss <- read.delim(paste(path,"/tmp/tsgtf/nrtss.annotation",sep=''),header = TRUE, sep = "\t",quote = "")
tss <- cbind(nrtss,e,j,matrix(0,nrow(nrtss),1))
colnames(tss)[(ncol(tss)-2):ncol(tss)] <- c("e_max","j_max", "result")

#### load in the pca and logit parameters
##pca
tmp <- read.delim(paste(syspath,"/model/tss_nuc_4.csv",sep=''),sep="\t", header=TRUE)
tp <- read.delim(paste(syspath,"/model/cage_0.1_tsscov.ns.cov",sep=''),sep="\t", header=FALSE)
fp <- read.delim(paste(syspath,"/model/cage_0.1-5_tsscov.ns.cov",sep=''),sep="\t", header=FALSE)
pca <- prcomp(tmp[,c("j_max","e_max")])
tmp <- cbind(tmp,pca$x)
a <- cbind(tss[,"j_max"]-pca$center["j_max"],tss[,"e_max"]-pca$center["e_max"]) %*% pca$rotation
tss <- cbind(tss,a)

##logit
tp <- tp[order(tp[,4]),]
fp <- fp[order(fp[,4]),]
tmp[,"result"] <- -1
tmp[which(table(c(as.matrix(tmp[,1]),as.matrix(fp[,4])))==2),"result"] <- 0
tmp[which(table(c(as.matrix(tmp[,1]),as.matrix(tp[,4])))==2),"result"] <- 1
tmp <- tmp[tmp[,"result"]> -1,]

a=glm(result~PC1+PC2,family=binomial(link=logit),data=tmp)



## apply logistic
pre=apply(cbind(predict(a,tss),rep(100,nrow(tss))),1,min)
p=exp(pre)/(1+exp(pre))
tss <- cbind(tss,p)
colnames(tss)[ncol(tss)] <- "logit"
tss$result <- tss$logit > 0.916

##check whether novel tss, true represent true TSS after filter, 0 represent novel (not present before assembly)
##first exons of each transcript
ant <- read.delim(paste(path,"/tmp/tsgtf/ant.annotation",sep=''), header=TRUE)
tss <- tss[order(tss[,1]),]
tss <- cbind(tss, rep(0,nrow(tss)))
tss[,ncol(tss)] <- table(factor(unique(ant[ant$class_code=="=","tss_id"]),tss$tss_id)) ##m1 represent known isoform
colnames(tss)[ncol(tss)] <- "known"
table(tss$result,tss$known)
write.table(tss, paste(path,"/tmp/tsgtf/nrtss_filter.annotation",sep=''),sep="\t", row.names=FALSE, col.name=TRUE, quote = FALSE)
write.table(tss[tss$result==1,], paste(path,"/FilteredNrtss.annotation",sep=''),sep="\t", row.names=FALSE, col.name=TRUE, quote = FALSE)

#Generate the bed format of each first exons used for DaPars input to detect the tandem 5'UTR
#Reverse the direction of each first exons to satisfy the requirements of DaPars model
tss <- tss[tss[,10]!=".",]
write.table(tss[tss$result==1,c(7,8,9,1,11,10)], paste(path,"/tmp/tsgtf/utr.bed",sep=''),sep="\t", row.names=FALSE, col.name=F, quote = FALSE)
write.table(tss[tss$result==1,c(7,8,9,1,11,10)], paste(path,"/utr.bed",sep=''),sep="\t", row.names=FALSE, col.name=F, quote = FALSE)
utr <- tss[tss$result==1,c(7,8,9,1,11,10)]
tmpp <- utr[utr[,6]=="+",]
tmpn <- utr[utr[,6]=="-",]
tmpp[,6] <- "-"
tmpn[,6] <- "+"
tmp <- rbind(tmpp,tmpn)
write.table(tmp, paste(path,"/tmp/tsgtf/utr_tmp.bed",sep=''), sep="\t", row.names=FALSE, col.name=F, quote = FALSE)

#nrtss_filter.annotation: All the first exons annotation combining with Logistic model predicted values
#FilteredNrtss.annotation: Real first exons identified by logistic model
#utr_tmp.bed: the regions of real first exons in reverse direction used to satisfy the requirements of DaPars model
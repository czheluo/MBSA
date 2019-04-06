#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
    'input','i',1,'character',
	'output','o',1,'character',
	'bs','b',2,'character',
	'ws','w',2,'character',
	'pta','c',2,'character',
	'ptb','d',2,'character',
	'ptc','e',2,'character',
	'ptd','f',2,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example:

Usage:
    --input input step01 dir
	--output	output dir
    --bs    bulk size
    --ws window size
	--pta pemutation test confidence interval
	--ptb pemutation test confidence interval
	--help		usage
\n")
	q(status=1);
}

if (is.null(opt$input)) { print_usage(spec)}
if (is.null(opt$output)){ print_usage(spec) }
#if (is.null(opt$hb)){ print_usage(spec) }
#if (is.null(opt$lb)){ print_usage(spec) }
#if (is.null(opt$wid)){ print_usage(spec) }
#if (is.null(opt$mid)){ print_usage(spec) }
if (is.null(opt$bs)){ opt$bs=30;print(opt$bs)}
if (is.null(opt$ws)){ opt$ws=1e6;print(opt$ws)}
if (is.null(opt$pta)){ opt$pta=95;print(opt$pta) }
if (is.null(opt$ptb)){ opt$ptb=99;print(opt$ptb)}
if (is.null(opt$ptc)){ opt$ptc=99.9;print(opt$ptb)}
if (is.null(opt$ptd)){ opt$ptd=99.99;print(opt$ptb)}
if (!is.null(opt$pop)){ print(opt$pop) }else{opt$pop="F2";print(opt$popt)}

times<-Sys.time()
library(magrittr)
library(dplyr)
library(tidyr)
library(parallel)
source("/mnt/ilustre/centos7users/dna/Pipeline/05.BSA/bin/R/takagi_sim.R")
source("/mnt/ilustre/centos7users/dna/Pipeline/05.BSA/bin/R/G_functions.R")
source("/mnt/ilustre/centos7users/dna/Pipeline/05.BSA/bin/R/RcppExports.R")
source("/mnt/ilustre/centos7users/dna/Pipeline/05.BSA/bin/R/export_functions.R")
library("Rcpp")
sourceCpp("/mnt/ilustre/centos7users/dna/Pipeline/05.BSA/bin/R/countSNPs.cpp")

setwd(opt$input)
load("step01.RData")

setwd(opt$output)

seq_G <- runGprimeAnalysis(SNPset = ntable,windowSize =opt$ws,outlierFilter = "deltaSNP")
write.table(seq_G,file = "G.analysis.result",sep = "\t",quote=F,row.names = FALSE)


##Gprime RESULT
##

getFDRThreshold <- function(pvalues, alpha = 0.01,method="BH"){
    sortedPvals <- sort(pvalues, decreasing = FALSE)
    pAdj <- p.adjust(sortedPvals, method = method)
    if (!any(pAdj < alpha)) {
        fdrThreshold <- NA
    } else {
    fdrThreshold <- sortedPvals[max(which(pAdj < alpha))]
    }
    return(fdrThreshold)
}
##change the Gprime which the nSNP less than 10
seq_G$nSNPs<-nSNPs
seq_G$Gprime[seq_G$nSNPs<10,]<-0

##three thershold
alpha<-c(0.05,0.01,0.001,0.0001)
##Threshold methods
#c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")

fdrT=NULL;GprimeT=NULL;
for (i in 1:length(alpha)){
	fdrT [i]<- getFDRThreshold(seq_G$pvalue, alpha = alpha[i],method="fdr")#bonferroni
	GprimeT[i]<- seq_G[which(seq_G$pvalue == fdrT[i]), "Gprime"]
	Gp<-matrix(GprimeT[i],length(seq_G$Gprime))
	dat1<-data.frame(seq_G$CHROM,seq_G$POS,seq_G$Gprime,Gp)
	colnames(dat1)<-c("Chr","Pos","Gprime","Gthreshold")
	write.table(dat1,file=paste(i,".Gprime.result",sep=""),sep="\t",quote=F,row.names=F,col.names=T)
}

#SNPset$Gqvalue[SNPset$Gqvalue==0]<-0.0001
#Gp<-matrix(GprimeT,length(seq_G$Gprime))
#colnames(Gp)<-"Gp"
#Gq<--log10(SNPset$Gqvalue)
#colnames(Gq)<-"Gq"
#dat1<-data.frame(seq_G$CHROM,seq_G$POS,seq_G$Gprime,Gp)
#colnames(dat1)<-c("Chr","Pos","Gprime","Gthreshold")
#write.table(dat1,file="Gprime.result",sep="\t",quote=F,row.names=F,col.names=T)

##Gprime PLOT
source("/mnt/ilustre/users/meng.luo/Pipeline/05.BSA/bin/R/CMplot.R")
SNP<-paste("SNP",c(1:length(ntable$CHROM)),sep="")
ch<-matrix(0,length(ntable$CHROM))
CHR<-unique(ntable$CHROM)
#ntable$CHROM<-factor(ntable$CHROM,levels =paste(unique(ntable$CHROM)))

for (i in 1:length(CHR)){
	ch[which(ntable$CHROM %in% CHR[i])]<-rep(paste(i),length(ntable$CHROM==CHR))
	}

Chr<-as.numeric(ch)

fdrT=NULL;GprimeT=NULL;
for (i in 1:4){

	fdrT [i]<- getFDRThreshold(seq_G$pvalue, alpha = alpha[i],method="fdr")#bonferroni
	GprimeT[i]<- seq_G[which(seq_G$pvalue == fdrT[i]), "Gprime"]
	Gp<-matrix(GprimeT[i],length(seq_G$Gprime))
	dat1<-data.frame(seq_G$CHROM,seq_G$POS,seq_G$Gprime,Gp)
	map<-data.frame(SNP,Chr,dat1[,-1])
	png(paste(i,"Gprime.fdr",".png",sep=""),width=1000, height=800)
	CMplot(map[c(1,2,3,4)], plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2, # ylim=c(-1,1),
		amplify=TRUE,main=expression("Gprime"),chr.labels=NULL,threshold=GprimeT[i],threshold.lty=2,threshold.lwd=1, threshold.col="red",
		bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)#,legend=expression(paste(GprimeT[i])),threLE=GprimeT[i])

	dev.off()
	pdf(paste(i,"Gprime.fdr",".pdf",sep=""),height=9,width=16)
	CMplot(map[c(1,2,3,4)], plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2, # ylim=c(-1,1),
		amplify=TRUE,main=expression("Gprime"),chr.labels=NULL,threshold=GprimeT[i],threshold.lty=2,threshold.lwd=1, threshold.col="red",
		bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)#,legend=expression(paste(GprimeT[i])),threLE=GprimeT[i])
	dev.off()

}

fdrT=NULL;GprimeT=NULL;
for (i in 1:4){

	fdrT [i]<- getFDRThreshold(seq_G$pvalue, alpha = alpha[i],method="bonferroni")#bonferroni
	GprimeT[i]<- seq_G[which(seq_G$pvalue == fdrT[i]), "Gprime"]
	Gp<-matrix(GprimeT[i],length(seq_G$Gprime))
	dat1<-data.frame(seq_G$CHROM,seq_G$POS,seq_G$Gprime,Gp)
	map<-data.frame(SNP,Chr,dat1[,-1])
	png(paste(i,"Gprime.bonf",".png",sep=""),width=1000, height=800)
	CMplot(map[c(1,2,3,4)], plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2, # ylim=c(-1,1),
		amplify=TRUE,main=expression("Gprime"),chr.labels=NULL,threshold=GprimeT[i],threshold.lty=2,threshold.lwd=1, threshold.col="red",
		bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)#,legend=expression(paste(GprimeT[i])),threLE=GprimeT[i])

	dev.off()
	pdf(paste(i,"Gprime.bonf",".pdf",sep=""),height=9,width=16)
	CMplot(map[c(1,2,3,4)], plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2, # ylim=c(-1,1),
		amplify=TRUE,main=expression("Gprime"),chr.labels=NULL,threshold=GprimeT[i],threshold.lty=2,threshold.lwd=1, threshold.col="red",
		bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)#,legend=expression(paste(GprimeT[i])),threLE=GprimeT[i])
	dev.off()
}

##DeltaSNP

setwd(paste(opt$input,"/Index/",sep="")

seq <- runQTLseqAnalysis(SNPset =ntable,windowSize = opt$ws,popStruc =paste(opt$pop),bulkSize =opt$bs,replications = 10000,intervals = c(opt$pta, opt$ptb,opt$ptc, opt$ptd))
write.table(seq,file = "qtlseq.analysis.result",sep = "\t",quote=F,row.names = FALSE)


##remove the nSNP less than 10

seq$tricubeDeltaSNP[seq$nSNP<10,]<-0


CI1<-abs(seq[,22])
CI2<-abs(seq[,23])
CI3<-abs(seq[,24])
CI4<-abs(seq[,25])

WindowDeltaSNP<-abs(seq_G$tricubeDeltaSNP)
dat1<-data.frame(seq_G$CHROM,seq_G$POS,WindowDeltaSNP,CI1)
colnames(dat1)<-c("Chr","Pos","WindowDeltaSNP","CI1") #opt$pta
write.table(dat1,file="Index1.result",sep="\t",quote=F,row.names=F,col.names=T)

dat2<-data.frame(seq_G$CHROM,seq_G$POS,WindowDeltaSNP,CI2)
colnames(dat2)<-c("Chr","Pos","WindowDeltaSNP","CI2") #opt$pta
write.table(dat2,file="Index2.result",sep="\t",quote=F,row.names=F,col.names=T)

dat3<-data.frame(seq_G$CHROM,seq_G$POS,WindowDeltaSNP,CI3)
colnames(dat3)<-c("Chr","Pos","WindowDeltaSNP","CI3") #opt$pta
write.table(dat3,file="Index3.result",sep="\t",quote=F,row.names=F,col.names=T)

dat4<-data.frame(seq_G$CHROM,seq_G$POS,WindowDeltaSNP,CI4)
colnames(dat4)<-c("Chr","Pos","WindowDeltaSNP","CI4") #opt$pta
write.table(dat4,file="Index4.result",sep="\t",quote=F,row.names=F,col.names=T)


map1<-data.frame(SNP,Chr,seq_G$POS,seq_G$tricubeDeltaSNP,seq[,22])

png(paste("DeltaSNP1",".png",sep=""),width=1000, height=800)

CMplot(map1[,c(1,2,3,4)], plot.type="m", LOG10=FALSE, line=2.5,cex.lab=2,ylim=c(round(min(map1[,4]),digits=1),round(max(map1[,4]),digits=1)),
       amplify=FALSE,main=expression("WindowDeltaSNP"),chr.labels=NULL,
       threshold=c(min(map1[,5]),-min(map1[,5])),threshold.lty=c(2,2),threshold.lwd=c(2,2), threshold.col=c("red","red"),
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)

dev.off()
pdf(paste(i,"DeltaSNP1",".pdf",sep=""),height=9,width=16)
CMplot(map1[,c(1,2,3,4)], plot.type="m", LOG10=FALSE, line=2.5,cex.lab=2,ylim=c(round(min(map1[,4]),digits=1),round(max(map1[,4]),digits=1)),
       amplify=FALSE,main=expression("WindowDeltaSNP"),chr.labels=NULL,
       threshold=c(min(map1[,5]),-min(map1[,5])),threshold.lty=c(2,2),threshold.lwd=c(2,2), threshold.col=c("red","red"),
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)

dev.off()


map2<-data.frame(SNP,Chr,seq_G$POS,seq_G$tricubeDeltaSNP,seq[,23])
pdf(paste(i,"DeltaSNP2",".pdf",sep=""),height=9,width=16)
CMplot(map2[,c(1,2,3,4)], plot.type="m", LOG10=FALSE, line=2.5,cex.lab=2,ylim=c(round(min(map2[,4]),digits=1),round(max(map2[,4]),digits=1)),
       amplify=FALSE,main=expression("WindowDeltaSNP"),chr.labels=NULL,
       threshold=c(min(map2[,5]),-min(map2[,5])),threshold.lty=c(2,2),threshold.lwd=c(2,2), threshold.col=c("red","red"),
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)
dev.off()


png(paste("DeltaSNP2",".png",sep=""),width=1000, height=800)
CMplot(map2[,c(1,2,3,4)], plot.type="m", LOG10=FALSE, line=2.5,cex.lab=2,ylim=c(round(min(map2[,4]),digits=1),round(max(map2[,4]),digits=1)),
       amplify=FALSE,main=expression("WindowDeltaSNP"),chr.labels=NULL,
       threshold=c(min(map2[,5]),-min(map2[,5])),threshold.lty=c(2,2),threshold.lwd=c(2,2), threshold.col=c("red","red"),
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)
dev.off()

map3<-data.frame(SNP,Chr,seq_G$POS,seq_G$tricubeDeltaSNP,seq[,24])
pdf(paste(i,"DeltaSNP3",".pdf",sep=""),height=9,width=16)
CMplot(map3[,c(1,2,3,4)], plot.type="m", LOG10=FALSE, line=2.5,cex.lab=2,ylim=c(round(min(map3[,4]),digits=1),round(max(map3[,4]),digits=1)),
       amplify=FALSE,main=expression("WindowDeltaSNP"),chr.labels=NULL,
       threshold=c(min(map3[,5]),-min(map3[,5])),threshold.lty=c(2,2),threshold.lwd=c(2,2), threshold.col=c("red","red"),
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)
dev.off()

png(paste("DeltaSNP3",".png",sep=""),width=1000, height=800)
CMplot(map3[,c(1,2,3,4)], plot.type="m", LOG10=FALSE, line=2.5,cex.lab=2,ylim=c(round(min(map3[,4]),digits=1),round(max(map3[,4]),digits=1)),
       amplify=FALSE,main=expression("WindowDeltaSNP"),chr.labels=NULL,
       threshold=c(min(map3[,5]),-min(map3[,5])),threshold.lty=c(2,2),threshold.lwd=c(2,2), threshold.col=c("red","red"),
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)
dev.off()


map4<-data.frame(SNP,Chr,seq_G$POS,seq_G$tricubeDeltaSNP,seq[,25])
pdf(paste(i,"DeltaSNP4",".pdf",sep=""),height=9,width=16)
CMplot(map4[,c(1,2,3,4)], plot.type="m", LOG10=FALSE, line=2.5,cex.lab=2,ylim=c(round(min(map4[,4]),digits=1),round(max(map4[,4]),digits=1)),
       amplify=FALSE,main=expression("WindowDeltaSNP"),chr.labels=NULL,
       threshold=c(min(map4[,5]),-min(map4[,5])),threshold.lty=c(2,2),threshold.lwd=c(2,2), threshold.col=c("red","red"),
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)
dev.off()
png(paste("DeltaSNP4",".png",sep=""),width=1000, height=800)
CMplot(map4[,c(1,2,3,4)], plot.type="m", LOG10=FALSE, line=2.5,cex.lab=2,ylim=c(round(min(map4[,4]),digits=1),round(max(map4[,4]),digits=1)),
       amplify=FALSE,main=expression("WindowDeltaSNP"),chr.labels=NULL,
       threshold=c(min(map4[,5]),-min(map4[,5])),threshold.lty=c(2,2),threshold.lwd=c(2,2), threshold.col=c("red","red"),
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)
dev.off()

escaptime<-Sys.time()-times;
print("Done!");
print(escaptime)


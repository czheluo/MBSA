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
	'pop','p',2,'character',
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
source("/mnt/ilustre/centos7users/meng.luo/Pipeline/05.BSA_new/bin/R/takagi_sim.R")
source("/mnt/ilustre/centos7users/meng.luo/Pipeline/05.BSA_new/bin/R/G_functions.R")
source("/mnt/ilustre/centos7users/meng.luo/Pipeline/05.BSA_new/bin/R/RcppExports.R")
source("/mnt/ilustre/centos7users/meng.luo/Pipeline/05.BSA_new/bin/R/export_functions.R")
library("Rcpp")
sourceCpp("/mnt/ilustre/centos7users/meng.luo/Pipeline/05.BSA_new/bin/R/countSNPs.cpp")

setwd(opt$input)
load("step01.RData")

##DeltaSNP
setwd(opt$output)

seq <- runQTLseqAnalysis(SNPset =ntable,windowSize = opt$ws,popStruc =paste(opt$pop),bulkSize =opt$bs,replications = 10000,intervals = c(opt$pta, opt$ptb,opt$ptc, opt$ptd))
write.table(seq,file = "qtlseq.analysis.result",sep = "\t",quote=F,row.names = FALSE)

##remove the nSNP less than 10

seq$tricubeDeltaSNP[ntable$nSNPs<10,]<-0

CI1<-abs(seq[,paste0("CI_",opt$pta)])
CI2<-abs(seq[,paste0("CI_",opt$ptb)])
CI3<-abs(seq[,paste0("CI_",opt$ptc)])
CI4<-abs(seq[,paste0("CI_",opt$ptd)])


WindowDeltaSNP<-abs(seq$tricubeDeltaSNP)

dir.create(file.path(opt$output, "Index1"), showWarnings = FALSE)
setwd(file.path(opt$output, "Index1"))

dat1<-data.frame(seq_G$CHROM,seq_G$POS,WindowDeltaSNP,CI1)
colnames(dat1)<-c("Chr","Pos","WindowDeltaSNP","CI1") #opt$pta
write.table(dat1,file="Index1.result",sep="\t",quote=F,row.names=F,col.names=T)

dir.create(file.path(opt$output, "Index2"), showWarnings = FALSE)
setwd(file.path(opt$output, "Index2"))

dat2<-data.frame(seq_G$CHROM,seq_G$POS,WindowDeltaSNP,CI2)
colnames(dat2)<-c("Chr","Pos","WindowDeltaSNP","CI2") #opt$pta
write.table(dat2,file="Index2.result",sep="\t",quote=F,row.names=F,col.names=T)

dir.create(file.path(opt$output, "Index3"), showWarnings = FALSE)
setwd(file.path(opt$output, "Index3"))

dat3<-data.frame(seq_G$CHROM,seq_G$POS,WindowDeltaSNP,CI3)
colnames(dat3)<-c("Chr","Pos","WindowDeltaSNP","CI3") #opt$pta
write.table(dat3,file="Index3.result",sep="\t",quote=F,row.names=F,col.names=T)

dir.create(file.path(opt$output, "Index4"), showWarnings = FALSE)
setwd(file.path(opt$output, "Index4"))

dat4<-data.frame(seq_G$CHROM,seq_G$POS,WindowDeltaSNP,CI4)
colnames(dat4)<-c("Chr","Pos","WindowDeltaSNP","CI4") #opt$pta
write.table(dat4,file="Index4.result",sep="\t",quote=F,row.names=F,col.names=T)

##plot
source("/mnt/ilustre/users/meng.luo/Pipeline/05.BSA/bin/R/CMplot.R")
SNP<-paste("SNP",c(1:length(ntable$CHROM)),sep="")
ch<-matrix(0,length(ntable$CHROM))
CHR<-unique(ntable$CHROM)
#ntable$CHROM<-factor(ntable$CHROM,levels =paste(unique(ntable$CHROM)))

for (i in 1:length(CHR)){
	ch[which(ntable$CHROM %in% CHR[i])]<-rep(paste(i),length(ntable$CHROM==CHR))
	}

Chr<-as.numeric(ch)

setwd(file.path(opt$output, "Index1"))
map1<-data.frame(SNP,Chr,seq_G$POS,seq_G$tricubeDeltaSNP,seq[,paste0("CI_",opt$pta)])

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

setwd(file.path(opt$output, "Index2"))
map2<-data.frame(SNP,Chr,seq_G$POS,seq_G$tricubeDeltaSNP,seq[,paste0("CI_",opt$ptb)])
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

setwd(file.path(opt$output, "Index3"))
map3<-data.frame(SNP,Chr,seq_G$POS,seq_G$tricubeDeltaSNP,seq[,paste0("CI_",opt$ptc])
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

setwd(file.path(opt$output, "Index3"))
map4<-data.frame(SNP,Chr,seq_G$POS,seq_G$tricubeDeltaSNP,seq[,paste0("CI_",opt$ptd]])
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
setwd(file.path(opt$output)
save.image("step04.RData")

escaptime<-Sys.time()-times;
print("Done!");
print(escaptime)


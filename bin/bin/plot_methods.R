#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'input','i',1,'character',
	'output','o',1,'character',
	'span','s',2,'character',
	'power','p',2,'character',
	'method','m',2,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example:
	Rscript G.R --input --hb --lb --output --all --opt
Usage:
	--input step01.RData file dir
	--output	output dir
	--span the span size if marker more than 1m,please,setting the span was bigger than 0.01
	--power ED power 1-4
  --method here two method to choose aicc(default) and gcv
	--help		usage
\n")
	q(status=1);
}
if (is.null(opt$input)) { print_usage(spec)}
if (is.null(opt$output)){ print_usage(spec) }
if (is.null(opt$span)){ opt$span=0.001;print(opt$span)}else{print(opt$span)}
if (is.null(opt$power)){ opt$power=4;print(opt$power)}else{print(opt$power)}
if (is.null(opt$method)){ opt$method="aicc";print(opt$method)}else{print(opt$method)}
times<-Sys.time()

setwd(opt$output)

load("/RData/step01.RData")



##plot

source("/mnt/ilustre/users/meng.luo/Pipeline/05.BSA/bin/R/CMplot.R")
SNP<-paste("SNP",c(1:length(SNPset$Chr)),sep="")
ch<-matrix(0,length(SNPset$Chr))
CHR<-unique(SNPset$Chr)
SNPset$Chr<-factor(SNPset$Chr,levels =paste(unique(SNPset$Chr)))

for (i in 1:length(CHR)){
	ch[which(SNPset$Chr %in% CHR[i])]<-rep(paste(i),length(SNPset$Chr==CHR))
	}

Chr<-as.numeric(ch)
map<-data.frame(SNP,Chr,SNPset[,-1])

png(paste("BSA_threemethods",".png",sep=""),width=1000, height=800)

par(mfrow=c(3,1))

CMplot(map[c(1,2,3,16)], plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2,#ylim=(c(round(min(map[,16]),digits=1),round(max(map[,16]),digits=1))),
       amplify=TRUE,main="ED4 losse fit",chr.labels=NULL,threshold=cutoff[1],threshold.lty=2,threshold.lwd=1, threshold.col="red",
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE,legend="cutoff",threLE=cutoff[1])

CMplot(map[c(1,2,3,12)], plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2, # ylim=c(-1,1),
       amplify=TRUE,main=expression("Gprime"),chr.labels=NULL,threshold=GprimeT,threshold.lty=2,threshold.lwd=1, threshold.col="red",
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE,legend=expression("q<-0.01"),threLE=GprimeT)

CMplot(map[c(1,2,3,10)], plot.type="m", LOG10=FALSE, ylim=c(round(min(map[,10]),digits=1),round(max(map[,10]),digits=1)),line=2.5,cex.lab=2,
       amplify=TRUE,main=expression("WindowDeltaSNP"),chr.labels=NULL,
       threshold=-max(map$I95),threshold.lty=2,threshold.lwd=1, threshold.col="red",
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE,legend="CI95",threLE=-0.38)

dev.off()


pdf(paste("BSA_threemethods","pdf",sep="."),height=9,width=16)
par(mfrow=c(3,1))

CMplot(map[c(1,2,3,16)], plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2,ylim=(c(round(min(map[,16]),digits=1),round(max(map[,16]),digits=1))),
       amplify=TRUE,main="ED4 losse fit",chr.labels=NULL,threshold=cutoff[1],threshold.lty=2,threshold.lwd=1, threshold.col="red",
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE,legend="cutoff",threLE=cutoff[1])

CMplot(map[c(1,2,3,12)], plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2, # ylim=c(-1,1),
       amplify=TRUE,main=expression("Gprime"),chr.labels=NULL,threshold=GprimeT,threshold.lty=2,threshold.lwd=1, threshold.col="red",
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE,legend=expression("q<-0.01"),threLE=GprimeT)

CMplot(map[c(1,2,3,10)], plot.type="m", LOG10=FALSE, ylim=c(round(min(map[,10]),digits=1),round(max(map[,10]),digits=1)),line=2.5,cex.lab=2,
       amplify=TRUE,main=expression("WindowDeltaSNP"),chr.labels=NULL,
       threshold=-max(map$I95),threshold.lty=2,threshold.lwd=1, threshold.col="red",
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE,legend="CI95",threLE=-0.38)

dev.off()

png(paste("BSA_index",".png",sep=""),width=1000, height=800)
par(mfrow=c(3,1))

CMplot(map[c(1,2,3,6)], plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2,ylim=c(0,1),
       amplify=TRUE,main="Indel1",chr.labels=NULL,
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)

CMplot(map[c(1,2,3,7)], plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2, ylim=c(0,1),
       amplify=TRUE,main=expression("Indel2"),chr.labels=NULL,
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)

CMplot(map[c(1,2,3,8)], plot.type="m", LOG10=FALSE, line=2.5,cex.lab=2,ylim=c(-1,1),
       amplify=TRUE,main=expression(Delta * '(SNP-index)'),chr.labels=NULL,
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)

dev.off()

pdf(paste("BSA_index","pdf",sep="."),height=9,width=16)
par(mfrow=c(3,1))

CMplot(map[c(1,2,3,6)], plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2,ylim=c(0,1),
       amplify=TRUE,main="Indel1",chr.labels=NULL,
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)

CMplot(map[c(1,2,3,7)], plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2, ylim=c(0,1),
       amplify=TRUE,main=expression("Indel2"),chr.labels=NULL,
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)

CMplot(map[c(1,2,3,8)], plot.type="m", LOG10=FALSE, line=2.5,cex.lab=2,ylim=c(-1,1),
       amplify=TRUE,main=expression(Delta * '(SNP-index)'),chr.labels=NULL,
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)

dev.off()

escaptime=Sys.time()-times;
print("Done!");
print(escaptime)


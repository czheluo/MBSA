#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
    'input','i',1,'character',
	'output','o',1,'character',
	'ws','w',2,'character',
	'method','m',2,'character',
	'alpha','a',2,'character',
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
if (is.null(opt$ws)){ opt$ws=1e6;print(opt$ws)}
if (!is.null(opt$methhod)){ print(opt$method) }else{opt$pop="bnf";print(opt$popt)}
if (!is.null(opt$alpha)){ print(opt$alpha) }else{opt$alpha<-c(0.05,0.01,0.001,0.0001);print(opt$alpha)}

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
seq_G$Gprime[ntable$nSNPs<10,]<-0
##three thershold
##Threshold methods
#c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")

fdrT=NULL;GprimeT=NULL;
if (opt$method=="bnf"){
	for (i in 1:length(opt$alpha)){
		dir.create(file.path(opt$output, paste("Gprime",i,sep="")), showWarnings = FALSE)
		setwd(file.path(opt$output, paste("Gprime",i,sep=""))
		fdrT [i]<- getFDRThreshold(seq_G$pvalue, alpha = opt$alpha[i],method="bonferroni")#bonferroni
		GprimeT[i]<- seq_G[which(seq_G$pvalue == fdrT[i]), "Gprime"]
		Gp<-matrix(GprimeT[i],length(seq_G$Gprime))
		dat1<-data.frame(seq_G$CHROM,seq_G$POS,seq_G$Gprime,Gp)
		colnames(dat1)<-c("Chr","Pos","Gprime","Gthreshold")
		write.table(dat1,file=paste(i,".Gprime.result",sep=""),sep="\t",quote=F,row.names=F,col.names=T)
	}
}else{
	for (i in 1:length(opt$alpha)){
		dir.create(file.path(opt$output, paste("Gprime",i,sep="")), showWarnings = FALSE)
		setwd(file.path(opt$output, paste("Gprime",i,sep=""))
		fdrT [i]<- getFDRThreshold(seq_G$pvalue, alpha = opt$alpha[i],method=opt$method)
		GprimeT[i]<- seq_G[which(seq_G$pvalue == fdrT[i]), "Gprime"]
		Gp<-matrix(GprimeT[i],length(seq_G$Gprime))
		dat1<-data.frame(seq_G$CHROM,seq_G$POS,seq_G$Gprime,Gp)
		colnames(dat1)<-c("Chr","Pos","Gprime","Gthreshold")
		write.table(dat1,file=paste(i,".Gprime.result",sep=""),sep="\t",quote=F,row.names=F,col.names=T)
	}
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
if (opt$method=="bnf"){
	for (i in 1:length(opt$alpha)){
		setwd(file.path(opt$output, paste("Gprime",i,sep=""))
		fdrT [i]<- getFDRThreshold(seq_G$pvalue, alpha = alpha[i],method="bonferroni")#bonferroni
		GprimeT[i]<- seq_G[which(seq_G$pvalue == fdrT[i]), "Gprime"]
		Gp<-matrix(GprimeT[i],length(seq_G$Gprime))
		dat1<-data.frame(seq_G$CHROM,seq_G$POS,seq_G$Gprime,Gp)
		map<-data.frame(SNP,Chr,dat1[,-1])
		png(paste(i,"Gprime.bonferroni",".png",sep=""),width=1000, height=800)
		CMplot(map[c(1,2,3,4)], plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2, # ylim=c(-1,1),
			amplify=TRUE,main=expression("Gprime"),chr.labels=NULL,threshold=GprimeT[i],threshold.lty=2,threshold.lwd=1, threshold.col="red",
			bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)#,legend=expression(paste(GprimeT[i])),threLE=GprimeT[i])

		dev.off()
		pdf(paste(i,"Gprime.bonferroni",".pdf",sep=""),height=9,width=16)
		CMplot(map[c(1,2,3,4)], plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2, # ylim=c(-1,1),
			amplify=TRUE,main=expression("Gprime"),chr.labels=NULL,threshold=GprimeT[i],threshold.lty=2,threshold.lwd=1, threshold.col="red",
			bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)#,legend=expression(paste(GprimeT[i])),threLE=GprimeT[i])
		dev.off()
}else{
	for (i in 1:length(opt$alpha)){
		setwd(file.path(opt$output, paste("Gprime",i,sep=""))
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
		pdf(paste(i,"Gprime.",opt$method,".pdf",sep=""),height=9,width=16)
		CMplot(map[c(1,2,3,4)], plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2, # ylim=c(-1,1),
			amplify=TRUE,main=expression("Gprime"),chr.labels=NULL,threshold=GprimeT[i],threshold.lty=2,threshold.lwd=1, threshold.col="red",
			bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)#,legend=expression(paste(GprimeT[i])),threLE=GprimeT[i])
		dev.off()
	}
}


setwd(file.path(opt$output)
save.image("step05.RData")

escaptime<-Sys.time()-times;
print("Done!");
print(escaptime)


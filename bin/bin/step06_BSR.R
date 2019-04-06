#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
  'input','i',1,'character',
  'output','o',1,'character',
  'tgl','s',2,'character',
  'bs','b',2,'character',
  'threshold',2,'character',
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
      --tgl  total genetic length (default was 2000)
      --bs Set N to be the number of individuals in the mutant pool
      --help		usage
      \n")
  q(status=1);
}

if (is.null(opt$input)) { print_usage(spec)}
if (is.null(opt$output)){ print_usage(spec) }
if (is.null(opt$tgl)){ opt$tgl=2000;print(opt$tgl)}else{print(opt$tgl)}
if (is.null(opt$bs)){ opt$bs=30;print(opt$bs)}else{print(opt$bs)}
if (is.null(opt$threshold)){ opt$threshold=0.999;print(opt$thresholds)}else{print(opt$threshold)}

times<-Sys.time()
setwd(opt$input)
load("step01.RData")
library(magrittr)
library(dplyr)
library(tidyr)
library(parallel)

setwd(opt$output)
### required data
#SNP data counts, 5-8 columuns are read counts for mutant-ref, mutant-alt, wildtype-ref and wildtype-alt
dat<-ltable[,c("Chrom","Pos","Ref","Alt",
               paste(HD,"1",sep = ""),paste(HD,"2",sep = ""),
               paste(LD,"1",sep = ""),paste(LD,"2",sep = ""))]
##Find estimated probability that a randomly selected
##SNP is in complete linkage disequilibrium with mutantgene.
dtor=function(dat)
{
  r = 0.5 * (1 - exp(-2 *dat))
  return(r)
}
cM=c(seq(0,20,by=0.01))
rf=dtor(cM/100)
pnr=(1-rf)^(2*opt$bs)
plot(cM,pnr)

##This is numerical integration of P(no recomb|dist)P(dist)
ptheta0=sum(2*(0.01/2000)*(pnr[-length(pnr)]+pnr[-1])/2)
ptheta0
##Find prior distribution of theta
thetahat=dat[,6]/apply(dat[,5:6],1,sum)
thetahat=c(thetahat,1-thetahat)
getppp=function(x)
{
  nm=sum(x[3:4])
  nw=sum(x[1:2])
  nonma=which.min(x[3:4])
  xm=min(x[3:4])
  xw=(x[1:2])[nonma]
  ##R="some recombinaton between SNP and causal gene in mutant pool"
  if(xm>0){
    ppmt=0
    ppwt=0
  }
  else{
    ##Find P(x1,x2|R)
    ##Integrate P(x1=0 or n|theta)prior(theta)
    p1=mean((1-thetahat)^nm+thetahat^nm)
    ##Find P(R|x1,x2)
    ppmt=(1/(1+p1*((1-ptheta0)/ptheta0)))
    ##Find P(wm,w) by integrating P(wm,w|theta)prior(theta)
    p2=mean(dbinom(xw,nw,thetahat))
    ##Find P(wm,w|R)
    p3=mean(dbinom(xw,nw,thetahat[thetahat>=0.5]))
    ##Find P(R|wm,w)
    ppwt=.5*p3/p2
    ##Find the product of the posterior probabilities.
  }
  ppp=ppmt*ppwt
  c(ppmt,ppwt,ppp)
}
#times<-Sys.time()
o=t(apply(dat[,5:8],1,getppp))
#otime<-Sys.time()-times
sortedPvals <- sort(o[,3], decreasing = FALSE)
ntable$ocutoff<-c(matrix(sortedPvals[round(length(sortedPvals)*opt$threshold)],length(ntable$POS)))


##BSR result

dat6<-data.frame(ntable$CHROM,ntable$POS,o[,3],ntable$ocutoff)
colnames(dat6)<-c("Chr","Pos","BSR","cutoff")
write.table(dat6,file="BSR.result",sep="\t",quote=F,row.names=F,col.names=T)

source("/mnt/ilustre/centos7users/dna/Pipeline/05.BSA/bin/R/CMplot.R")
SNP<-paste("SNP",c(1:length(dat4$Chr)),sep="")
ch<-matrix(0,length(dat4$Chr))
CHR<-unique(dat4$Chr)
dat1$Chr<-factor(dat4$Chr,levels =paste(unique(dat4$Chr)))

for (i in 1:length(CHR)){
	ch[which(dat4$Chr %in% CHR[i])]<-rep(paste(i),length(dat4$Chr==CHR))
	}

Chr<-as.numeric(ch)

map6<-data.frame(SNP,Chr,ntable$POS,o[,3])

png(paste("BSR",".png",sep=""),width=1000, height=500)

CMplot(map6, plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2,#ylim=(c(round(min(map[,16]),digits=1),round(max(map[,16]),digits=1))),amplify=FALSE,
       main="Probability",chr.labels=NULL,threshold=ntable$ocutoff[1],threshold.lty=2,threshold.lwd=1, threshold.col="red",
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)

dev.off()
pdf(paste(i,"BSR",".pdf",sep=""),height=9,width=16)
CMplot(map6, plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2,#ylim=(c(round(min(map[,16]),digits=1),round(max(map[,16]),digits=1))),amplify=FALSE,
       main="Probability",chr.labels=NULL,threshold=ntable$ocutoff[1],threshold.lty=2,threshold.lwd=1, threshold.col="red",
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)

dev.off()

#manhattan(dat4, main = "Manhattan Plot",logp = FALSE,#ylim=0.3,
  #  cex.axis = 0.9, col = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3"),
    #suggestiveline = F, genomewideline = F)

##Columns are
##posterior probability from mutant data,
##posterior probability from wild type data,
##product of the posterior probabilities.
setwd(file.path(opt$output)
save.image("step06.RData")

escaptime=Sys.time()-times;
print("Done!");
print(escaptime)

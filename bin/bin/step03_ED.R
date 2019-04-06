#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'input','i',1,'character',
	'output','o',1,'character',
	'span','s',2,'character',
	#'power','p',2,'character',
	'ws','w',2,'character',
	'fold','f',2,'character',
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
#if (is.null(opt$power)){ opt$power=4;print(opt$power)}else{print(opt$power)}
if (is.null(opt$fold)){ opt$fold=3;print(opt$fold)}else{print(opt$fold)}
if (is.null(opt$method)){ opt$method="aicc";print(opt$method)}else{print(opt$method)}
if (is.null(opt$ws)){ opt$ws=1e6;print(opt$ws)}
times<-Sys.time()
setwd(opt$input)
load("step01.RData")
setwd(opt$output)
library(magrittr)
library(dplyr)
library(tidyr)
library(parallel)
source("/mnt/ilustre/centos7users/meng.luo/Pipeline/05.BSA_new/bin/R/takagi_sim.R")
source("/mnt/ilustre/centos7users/meng.luo/Pipeline/05.BSA_new/bin/R/G_functions.R")
source("/mnt/ilustre/centos7users/meng.luo/Pipeline/05.BSA_new/bin/R/RcppExports.R")
source("/mnt/ilustre/centos7users/meng.luo/Pipeline/05.BSA_new/bin/R/export_functions.R")
library("Rcpp")
#sourceCpp("/mnt/ilustre/centos7users/meng.luo/Pipeline/05.BSA_new/bin/R/countSNPs.cpp")

loessopt <- function (s, e, p, m) {
  x <- try(loess(e ~ p, span=s, degree=1, family="symmetric", surface='direct'), silent=T)
  if(class(x)=="try-error"){return(NA)}
  span <- x$pars$span
  n <- x$n
  traceL <- x$trace.hat
  sigma2 <- sum( x$residuals^2 ) / (n-1)
  delta1 <- x$one.delta
  delta2 <- x$two.delta
  enp <- x$enp
  if(m=="aicc"){return(log(sigma2) + 1 + 2* (2*(traceL+1)) / (n-traceL-2))}
  if(m=="gcv"){return(n*sigma2 / (n-traceL)^2)}
}

#ED[ntable$nSNPs<10,]<-0
## six node
cl <- makeCluster(as.numeric(6))
#ED_INDEX<-sqrt((seq_G$SNPindex.LOW-seq_G$SNPindex.HIGH)^2)
chrs <- unique(ntable$CHROM)
#chrs <- chrs[grep(chrname, chrs)]
#ED<-ED_INDEX

time<-Sys.time()

for(chr in chrs) {
  e<-ED[ntable$CHROM==chr]^4
  p<-ntable$POS[ntable$CHROM==chr]
  spanresults.span <- NULL
  spanresults.aicc <- NULL
  #if(as.integer(length(p)) < 100){
    #cat(chr, " has less than 100 snps (n=", length(p), ")\n", sep="")
   # next
  #}
  ###change AIC for more running times
  if (length(p)<50){
	spanresults.span <- seq(round(1/length(p),digits=3),1, 0.001)#2*round(50/length(p),digits=3)
	spanresults.aicc <- parSapply(cl, spanresults.span, loessopt, e=e, p=p, m="aicc")
	usespan <- spanresults.span[spanresults.aicc==min(spanresults.aicc, na.rm=T)]
	if (length(usespan)>1){
	usespan<-max(usespan,na.rm=T)
	}
	lo <- loess(e ~ p, span = usespan[1], degree=1, family='symmetric', surface='direct')
	ntable$fitted[ntable$CHROM==chr] <- lo$fitted
	ntable$unfitted[ntable$CHROM==chr] <- lo$y
	cat(chr, length(p), round(usespan[1], digits=3), "\n", sep="\t")
  }else if (length(p)>50 && length(p)<5001){
	spanresults.span <- seq(round(50/length(p),digits=3),1, 0.001)#2*round(50/length(p),digits=3)
	spanresults.aicc <- parSapply(cl, spanresults.span, loessopt, e=e, p=p, m="aicc")
	usespan <- spanresults.span[spanresults.aicc==min(spanresults.aicc, na.rm=T)]
	if (length(usespan)>1){
	usespan<-max(usespan,na.rm=T)
	}
	lo <- loess(e ~ p, span = usespan[1], degree=1, family='symmetric', surface='direct')
	ntable$fitted[ntable$CHROM==chr] <- lo$fitted
	ntable$unfitted[ntable$CHROM==chr] <- lo$y
	cat(chr, length(p), round(usespan[1], digits=3), "\n", sep="\t")
  }else{
  	spanresults.span <- seq(round(length(p)/1.2/length(p),digits=3),1, 0.001)#2*round(50/length(p),digits=3)
	spanresults.aicc <- parSapply(cl, spanresults.span, loessopt, e=e, p=p, m="aicc")
	usespan <- spanresults.span[spanresults.aicc==min(spanresults.aicc, na.rm=T)]
	if (length(usespan)>1){
	usespan<-max(usespan,na.rm=T)
	}
	lo <- loess(e ~ p, span = usespan[1], degree=1, family='symmetric', surface='direct')
	ntable$fitted[ntable$CHROM==chr] <- lo$fitted
	ntable$unfitted[ntable$CHROM==chr] <- lo$y
	cat(chr, length(p), round(usespan[1], digits=3), "\n", sep="\t")

  }
}

#stopCluster(cl = NULL)
times<-Sys.time()-time

save.image("step03.RData")

fold<-as.numeric(opt$fold)
power<-as.numeric(opt$power)
ntable$fitted[ntable$fitted<0]<-0
cutoff <- fold*(sd(ntable$fitted)+median(ntable$fitted))
ntable$cutoff<-c(matrix(cutoff,length(ntable$fitted)))
EDpower<-ED^power

dir.create(file.path(opt$output, "ED1"), showWarnings = FALSE)
setwd(file.path(opt$output, "ED1")

##ED result

dat5<-data.frame(ntable$CHROM,ntable$POS,ntable$fitted,ntable$cutoff)
colnames(dat5)<-c("Chr","Pos","ED4_fit","cutoff")
write.table(dat5,file="ED1.result",sep="\t",quote=F,row.names=F,col.names=T)

dir.create(file.path(opt$output, "ED2"), showWarnings = FALSE)
setwd(file.path(opt$output, "ED2")

ntable$fitted[ntable$fitted<0]<-0
cutoff1 <- fold*4*(sd(ntable$fitted)+median(ntable$fitted))
ntable$cutoff1<-c(matrix(cutoff1,length(ntable$fitted)))

dat6<-data.frame(ntable$CHROM,ntable$POS,ntable$fitted,ntable$cutoff1)
colnames(dat6)<-c("Chr","Pos","ED4_fit","cutoff")
write.table(dat6,file="ED2.result",sep="\t",quote=F,row.names=F,col.names=T)



source("/mnt/ilustre/centos7users/meng.luo/Pipeline/05.BSA_new/bin/R/CMplot.R")
SNP<-paste("SNP",c(1:length(dat5$Chr)),sep="")
ch<-matrix(0,length(dat5$Chr))
CHR<-unique(dat5$Chr)
dat5$Chr<-factor(dat5$Chr,levels =paste(unique(dat5$Chr)))

for (i in 1:length(CHR)){
	ch[which(dat5$Chr %in% CHR[i])]<-rep(paste(i),length(dat5$Chr==CHR))
	}

Chr<-as.numeric(ch)

##OR grun

setwd(file.path(opt$output, "ED1")

map3<-data.frame(SNP,Chr,ntable$POS,ntable$fitted)

png(paste("ED_Losse_Fit",".png",sep=""),width=1000, height=500)

CMplot(map3, plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2,#ylim=(c(round(min(map[,16]),digits=1),round(max(map[,16]),digits=1))),amplify=FALSE,
       main="losse fit of ED",chr.labels=NULL,threshold=ntable$cutoff[1],threshold.lty=2,threshold.lwd=1, threshold.col="red",
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)

dev.off()
pdf(paste("ED_Losse_Fit",".pdf",sep=""),height=9,width=16)

CMplot(map3, plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2,#ylim=(c(round(min(map[,16]),digits=1),round(max(map[,16]),digits=1))),amplify=FALSE,
       main="losse fit of ED",chr.labels=NULL,threshold=ntable$cutoff[1],threshold.lty=2,threshold.lwd=1, threshold.col="red",
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)

dev.off()

setwd(file.path(opt$output, "ED2")

map4<-data.frame(SNP,Chr,ntable$POS,ntable$fitted)

png(paste("ED_Losse_Fit",".png",sep=""),width=1000, height=500)

CMplot(map4, plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2,#ylim=(c(round(min(map[,16]),digits=1),round(max(map[,16]),digits=1))),amplify=FALSE,
       main="losse fit of ED",chr.labels=NULL,threshold=ntable$cutoff1[1],threshold.lty=2,threshold.lwd=1, threshold.col="red",
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)

dev.off()
pdf(paste("ED_Losse_Fit",".pdf",sep=""),height=9,width=16)

CMplot(map4, plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2,#ylim=(c(round(min(map[,16]),digits=1),round(max(map[,16]),digits=1))),amplify=FALSE,
       main="losse fit of ED",chr.labels=NULL,threshold=ntable$cutoff1[1],threshold.lty=2,threshold.lwd=1, threshold.col="red",
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)

dev.off()


####FITTING WITH DIFFERENT REGRESSION

tricubeStat <- function(POS, Stat, windowSize = 2e6, ...)
{
    if (windowSize <= 0)
        stop("A positive smoothing window is required")
    stats::predict(locfit::locfit(Stat ~ locfit::lp(POS, h = windowSize, deg = 0),maxk=500, ...), POS)
}

tED <- tricubeStat(POS = ntable$POS,Stat = EDpower,windowSize = opt$ws)

## cutoff value use more than 3 times than default value (was three)
dir.create(file.path(opt$output, "ED3"), showWarnings = FALSE)
setwd(file.path(opt$output, "ED3")
cutofft <- opt$fold*(sd(tED)+median(tED))## 12*(sd(tED)+median(tED))

cutofft<-c(matrix(cutofft,length(tED)))

dat6<-data.frame(ntable$CHROM,ntable$POS,tED,cutofft)
colnames(dat6)<-c("Chr","Pos","ED4_t","cutofft")
write.table(dat6,file="ED3.result",sep="\t",quote=F,row.names=F,col.names=T)

#result<-
#write.table(result,file="three.methods.result",sep="\t",quote=F,row.names=F,col.names=T)

source("/mnt/ilustre/centos7users/meng.luo/Pipeline/05.BSA_new/bin/R/CMplot.R")
SNP<-paste("SNP",c(1:length(dat5$Chr)),sep="")
ch<-matrix(0,length(dat5$Chr))
CHR<-unique(dat5$Chr)
dat1$Chr<-factor(dat5$Chr,levels =paste(unique(dat5$Chr)))

for (i in 1:length(CHR)){
	ch[which(dat5$Chr %in% CHR[i])]<-rep(paste(i),length(dat5$Chr==CHR))
	}

Chr<-as.numeric(ch)
##OR grun
map5<-data.frame(SNP,Chr,ntable$POS,ntable$fitted)

png(paste("ED_Losse_Fit",".png",sep=""),width=1000, height=500)
CMplot(map5, plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2,#ylim=(c(round(min(map[,16]),digits=1),round(max(map[,16]),digits=1))),amplify=FALSE,
       main="losse fit of ED",chr.labels=NULL,threshold=cutofft[1],threshold.lty=2,threshold.lwd=1, threshold.col="red",
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)


dev.off()

pdf(paste("ED_Losse_Fit",".pdf",sep=""),height=9,width=16)
CMplot(map5, plot.type="m", LOG10=FALSE,line=2.5,cex.lab=2,#ylim=(c(round(min(map[,16]),digits=1),round(max(map[,16]),digits=1))),amplify=FALSE,
       main="losse fit of ED",chr.labels=NULL,threshold=cutofft[1],threshold.lty=2,threshold.lwd=1, threshold.col="red",
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)


dev.off()

setwd(file.path(opt$output)
save.image("step03.RData")

escaptime=Sys.time()-times;
print("Done!");
print(escaptime)

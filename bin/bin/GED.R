#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'input','i',1,'character',
	'hb','p',1,'character',
	'lb','l',1,'character',
	'pop','t',2,'character',
	'bs','b',2,'character',
	'ws','w',2,'character',
	'pta','c',2,'character',
	'ptb','d',2,'character',
	'output','o',1,'character',
	'opt','u','2','character',
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
	--input *vcf.table (the vcf to TABLE from GATK table)
	--chr the number of chromosomes
	--hb	HighBulk sampleID
	--lb	LowBulk sampleID
	--pop   all kinds of population
	--bs    bulk size
	--pta pemutation test confidence interval
	--ptb pemutation test confidence interval
	--output	output dir
	--ws window size
	--opt save the file name of the p-value less 0.01
	--help		usage
\n")
	q(status=1);
}

times<-Sys.time()

if (is.null(opt$input)) { print_usage(spec)}
if (is.null(opt$output)){ print_usage(spec) }
if (is.null(opt$hb)){ print_usage(spec) }
if (is.null(opt$lb)){ print_usage(spec) }
if (is.null(opt$bs)){ opt$bs=30 }
print(opt$bs)
if (is.null(opt$ws)){ opt$ws=1e6 }

print(opt$ws)
if (is.null(opt$pta)){ opt$pta=95 }
print(opt$pta)
if (is.null(opt$ptb)){ opt$ptb=99 }
print(opt$ptb)
if (is.null(opt$pop)){ opt$pop="F2" }
print(opt$pop)

times<-Sys.time()

#library(QTLseqr)
library(magrittr)
library(dplyr)
library(tidyr)
source("/mnt/ilustre/users/meng.luo/Pipeline/05.BSA/bin/R/takagi_sim.R")
source("/mnt/ilustre/users/meng.luo/Pipeline/05.BSA/bin/R/G_functions.R")
source("/mnt/ilustre/users/meng.luo/Pipeline/05.BSA/bin/R/RcppExports.R")
library(Rcpp)
sourceCpp("/mnt/ilustre/users/meng.luo/Pipeline/05.BSA/bin/R/countSNPs.cpp")

##ONE WAY
#library(magrittr)
#pwd<-setwd(opt$output)
#source(sprintf("%s%s",pwd,"/importFromGATK.R"))
#source(sprintf("%s%s",pwd,"/importfilter.R")
#source("/mnt/ilustre/users/meng.luo/QTLseqR/BSA_G_ED/bin/importFromGATK.R")
#source("/mnt/ilustre/users/meng.luo/QTLseqR/BSA_G_ED/bin/importfilter.R")
##TWO WAY
setwd(opt$output)
#chr<-opt$chr
#Chroms <- paste0(rep("Chr",chr), 1:chr)
#windowSize<- opt$ws
##Import SNP data from GATKTable file
HighBulk <- opt$hb
LowBulk <- opt$lb
#data <-importFromGATK(file =opt$input,highBulk = opt$hb,lowBulk = opt$lb)
#times=Sys.time()
table<-read.table(file=opt$input,header=T)
#caltime=Sys.time()-times


colnames(table)<-c("CHROM","POS","REF","ALT","HIGH","DP.HIGH","GQ.HIGH","PL.HIGH","LOW","DP.LOW","GQ.LOW","PL.LOW")
table<-table%>%tidyr::separate(col = HIGH,into = "AD_REF.HIGH",sep = ",",extra = "drop",convert = TRUE)
table<-table%>%tidyr::separate(col = LOW,into = "AD_REF.LOW",sep = ",",extra = "drop",convert = TRUE)
table<-table%>%dplyr::mutate(
            AD_ALT.HIGH = DP.HIGH - AD_REF.HIGH,
            AD_ALT.LOW = DP.LOW - AD_REF.LOW,
            SNPindex.HIGH = AD_ALT.HIGH / DP.HIGH,
            SNPindex.LOW = AD_ALT.LOW / DP.LOW,
            REF_FRQ = (AD_REF.HIGH + AD_REF.LOW) / (DP.HIGH + DP.LOW),
            deltaSNP = SNPindex.HIGH - SNPindex.LOW
        )
table<-table%>%dplyr::select(
            -dplyr::contains("HIGH"),
            -dplyr::contains("LOW"),
            -dplyr::one_of("deltaSNP", "REF_FRQ"),
            dplyr::matches("AD.*.LOW"),
            dplyr::contains("LOW"),
            dplyr::matches("AD.*.HIGH"),
            dplyr::contains("HIGH"),
            dplyr::everything()
        )

## if else for another way to filter data
#data<-filterSNPs(SNPset = data,refAlleleFreq = 0.20,minTotalDepth = 100,maxTotalDepth = 400,minSampleDepth = 40,minGQ = 99)
#write.table(data,file = "filter.index.result",sep = "\t",quote=F,row.names = FALSE)
seq_G <- runGprimeAnalysis(SNPset = table,windowSize =opt$ws,outlierFilter = "deltaSNP")
write.table(seq_G,file = "G.analysis.result",sep = "\t",quote=F,row.names = FALSE)

seq <- runQTLseqAnalysis(SNPset =table,windowSize = opt$ws,popStruc =paste(opt$pop),bulkSize =opt$bs,replications = 10000,intervals = c(opt$pta, opt$ptb))
write.table(seq,file = "qtlseq.analysis.result",sep = "\t",quote=F,row.names = FALSE)


times<-Sys.time()
library(parallel)
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
## three node
c1 <- makeCluster(as.numeric(3))
A<-seq_G$AD_REF.LOW/(seq_G$AD_REF.LOW+seq_G$AD_REF.HIGH)
B<-seq_G$AD_REF.HIGH/(seq_G$AD_REF.LOW+seq_G$AD_REF.HIGH)
C<-seq_G$AD_ALT.LOW/(seq_G$AD_ALT.HIGH+seq_G$AD_ALT.LOW)
D<-seq_G$AD_ALT.HIGH/(seq_G$AD_ALT.HIGH+seq_G$AD_ALT.LOW)
ED_AD<-sqrt((B-A)^2+(D-C)^2)
ED_DP<-sqrt((seq_G$DP.HIGH -seq_G$DP.LOW)^2)
mat<-as.matrix(cbind(seq_G$AD_REF.LOW,seq_G$AD_REF.HIGH,
                     seq_G$AD_ALT.LOW,seq_G$AD_ALT.HIGH))
n<-dim(mat)[1];m<-dim(mat)[2];ED_ADG<-NULL
for (i in 1:n) {
  ED_ADG[i]<-sqrt((mat[i,2]-mat[i,1]) ^ 2+(mat[i,4]-mat[i,3])^2)
}
ED_INDEX<-sqrt((seq_G$SNPindex.LOW-seq_G$SNPindex.HIGH)^2)
chrs <- unique(seq_G$CHROM)
#chrs <- chrs[grep(chrname, chrs)]
power<-4
m<-"aicc"
#ED<-ED_ADG
ED<-ED_INDEX
#ED<-ED_AD
#ED<-ED_DP
times<-Sys.time()
for(chr in chrs) {
  e<-ED[seq_G$CHROM==chr]^power
  p<-seq_G$POS[seq_G$CHROM==chr]
  spanresults.span <- NULL
  spanresults.aicc <- NULL
  if(as.integer(length(p)) < 100){
    cat(chr, " has less than 100 snps (n=", length(p), ")\n", sep="")
    next
  }
  ###change AIC for more running times
  spanresults.span <- seq(round(50/length(p),digits=3),1, .01)#2*round(50/length(p),digits=3)
  spanresults.aicc <- parSapply(c1, spanresults.span, loessopt, e=e, p=p, m=m)
  usespan <- spanresults.span[spanresults.aicc==min(spanresults.aicc, na.rm=T)]
  lo <- loess(e ~ p, span = usespan[1], degree=1, family='symmetric', surface='direct') #designate usespan index incase of tie.
  seq_G$fitted[seq_G$CHROM==chr] <- lo$fitted
  #seq_G$unfitted[seq_G$CHROM==chr] <- lo$y
  #cat(chr, length(p), round(usespan[1], digits=3), "\n", sep="\t")
}
losefittime<-Sys.time()-times;print(losefittime);
seq_G$fitted[seq_G$fitted<0]<-0
cutoff <- 3*(sd(seq_G$fitted)+median(seq_G$fitted))
cutoff<-c(matrix(cutoff,length(seq_G$fitted)))
#seq_G$cutoff<-cutoff
ED4<-ED_AD^4

escaptime=Sys.time()-times;

result<-cbind(seq_G[,c(1:4)],seq_G[,c(10,16,18,19,20)],seq[,23],seq_G[,c(22,23,25)],ED4,seq_G[,26],cutoff)
#result<-cbind(seq_G[,c(1:4)],seq_G[,c(18,24,26,27,28)],seq[,c(31)],seq_G[,c(30,31,33)],ED4,seq_G[,c(34)],cutoff)
colnames(result)<-c("Chr","Pos","Ref","Alt","index.low",
	"index.high","delta","nSNP","WindowDeltaSNP",paste("I",95,sep=""),"Gprime","Gpvalue",
	"Gqvalue","ED4","ED4_fit","ED4_fit_cutoff")
write.table(result,file="three.methods.result",sep="\t",quote=F,row.names=F,col.names=T)
##plot
caltime<-Sys.time()-times


##result the result

SNPset<-result
##ED
dat1<-data.frame(SNPset$Chr,SNPset$Pos,SNPset$ED4_fit,SNPset$ED4_fit_cutoff)
colnames(dat1)<-c("Chr","Pos","ED4_fit","cutoff")
write.table(dat1,file="ED.result",sep="\t",quote=F,row.names=F,col.names=T)
##Gprime
getFDRThreshold <- function(pvalues, alpha = 0.01){
    sortedPvals <- sort(pvalues, decreasing = FALSE)
    pAdj <- p.adjust(sortedPvals, method = "BH")
    if (!any(pAdj < alpha)) {
        fdrThreshold <- NA
    } else {
    fdrThreshold <- sortedPvals[max(which(pAdj < alpha))]
    }
    return(fdrThreshold)
}

fdrT <- getFDRThreshold(SNPset$Gpvalue, alpha = 0.01)
GprimeT <- SNPset[which(SNPset$Gpvalue == fdrT), "Gprime"]
#SNPset$Gqvalue[SNPset$Gqvalue==0]<-0.0001
Gp<-matrix(GprimeT,length(SNPset$Gprime))
#colnames(Gp)<-"Gp"
#Gq<--log10(SNPset$Gqvalue)
#colnames(Gq)<-"Gq"
dat2<-data.frame(SNPset$Chr,SNPset$Pos,SNPset$Gprime,Gp)
colnames(dat2)<-c("Chr","Pos","Gprime","Gthreshold")
write.table(dat2,file="Gprime.result",sep="\t",quote=F,row.names=F,col.names=T)

##DeltaSNP
SNPset[,9]<-abs(SNPset[,9])
SNPset[,10]<-abs(SNPset[,10])

dat3<-data.frame(SNPset$Chr,SNPset$Pos,SNPset$WindowDeltaSNP,SNPset[,10])
colnames(dat3)<-c("Chr","Pos","WindowDeltaSNP",paste("I",95,sep="")) #opt$pta
write.table(dat3,file="Index.result",sep="\t",quote=F,row.names=F,col.names=T)


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


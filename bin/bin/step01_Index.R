#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'input','i',1,'character',
	'hb','p',1,'character',
	'lb','l',1,'character',
	'pop','t',1,'character',
	'wid','s',2,'character',
	'mid','y',2,'character',
	'output','o',1,'character',
	'opt','u','2','character',
	'ws','w','2','character',
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
	--input *vcf.table (the vcf to TABLE from GATK table with absoluted)
	--hb	HighBulk sampleID
	--lb	LowBulk sampleID
	--pop   all kinds of population
	--wid wild parent name
	--mid mutant parent name
	--output	output dir
  --ws window size default was 1M
	--opt save the file name of the p-value less 0.01
	--help		usage
\n")
	q(status=1);
}

if (is.null(opt$input)) { print_usage(spec)}
if (is.null(opt$output)){ print_usage(spec) }
if (is.null(opt$hb)){ print_usage(spec) }
if (is.null(opt$lb)){ print_usage(spec) }
if (!is.null(opt$pop)){ print(opt$pop) }else{opt$pop="F2";print(opt$pop)}
if (!is.null(opt$ws)){ print(opt$ws) }else{opt$ws=1e6;print(opt$ws)}
times<-Sys.time()
setwd(opt$output)
library(magrittr)
library(dplyr)
library(tidyr)

##read.vcf.table
table<-read.table(opt$input,header = T)
ltable<-table[c(1:17)];
popt<-opt$pop

if (!is.null(opt$wid)){
	WID<-opt$wid
	WT<-paste(WID,".GT",sep = "")
	ltable<-ltable%>%tidyr::separate(col = WT,
                                 into = c(paste(WT,"1",sep = ""),paste(WT,"2",sep = "")),
                                 sep = "/",convert = TRUE)#extra = "drop"
	}else{next}

if (!is.null(opt$mid)){
	MID<-opt$mid
	MT<-paste(MID,".GT",sep = "")
	ltable<-ltable%>%tidyr::separate(col = MT,
                                 into = c(paste(MT,"1",sep = ""),paste(MT,"2",sep = "")),
                                 sep = "/",convert = TRUE)#extra = "drop"
	}else{next}

if (!is.null(opt$hb)){
	HighBulk <- opt$hb
	HD<-paste(HighBulk,".AD",sep = "")
	HT<-paste(HighBulk,".GT",sep = "")
	ltable<-ltable%>%tidyr::separate(col = HD,
                                 into = c(paste(HD,"1",sep = ""),paste(HD,"2",sep = "")),
                                 sep = ",",convert = TRUE)#extra = "drop"
	ltable<-ltable%>%tidyr::separate(col = HT,
                                 into = c(paste(HT,"1",sep = ""),paste(HT,"2",sep = "")),
                                 sep = "/",convert = TRUE)#extra = "drop"
	}else{next}

if (!is.null(opt$lb)){
	LowBulk <- opt$lb
	LD<-paste(LowBulk ,".AD",sep = "")
	LT<-paste(LowBulk ,".GT",sep = "")
	ltable<-ltable%>%tidyr::separate(col = LD,
                                 into = c(paste(LD,"1",sep = ""),paste(LD,"2",sep = "")),
                                 sep = ",",convert = TRUE)#extra = "drop"
	ltable<-ltable%>%tidyr::separate(col = LT,
                                 into = c(paste(LT,"1",sep = ""),paste(LT,"2",sep = "")),
                                 sep = "/",convert = TRUE)#extra = "drop"
	}else{next}

Index1<-NULL;Index2<-NULL;detale<-NULL
for(i in 1:length(table[,1])){
  if (popt=="F2") {
      if (ltable[i,paste(HT,"1",sep = "")]==ltable[i,paste(MT,"1",sep = "")]){
        Index1[i]<-ltable[i,paste(HD,"2",sep = "")]/(ltable[i,paste(HD,"1",sep = "")]+ltable[i,paste(HD,"2",sep = "")])
      }else{
        Index1[i]<-ltable[i,paste(HD,"2",sep = "")]/(ltable[i,paste(HD,"1",sep = "")]+ltable[i,paste(HD,"2",sep = "")])
      }
    if (ltable[i,paste(LT,"1",sep = "")]==ltable[i,paste(WT,"1",sep = "")]) {
      Index2[i]<-ltable[i,paste(LD,"2",sep = "")]/(ltable[i,paste(LD,"1",sep = "")]+ltable[i,paste(LD,"2",sep = "")])
    }else{
      Index2[i]<-ltable[i,paste(LD,"2",sep = "")]/(ltable[i,paste(LD,"1",sep = "")]+ltable[i,paste(LD,"2",sep = "")])
    }

  }else if(is.null(PID)){
    if (ltable[i,paste(HT,"1",sep = "")]==ltable[i,paste(WT,"1",sep = "")]){
      Index1[i]<-ltable[i,paste(HD,"1",sep = "")]/(ltable[i,paste(HD,"1",sep = "")]+ltable[i,paste(HD,"2",sep = "")])
    }else{
      Index1[i]<-ltable[i,paste(HD,"2",sep = "")]/(ltable[i,paste(HD,"1",sep = "")]+ltable[i,paste(HD,"2",sep = "")])
    }
    if (ltable[i,paste(LT,"1",sep = "")]==ltable[i,paste(WT,"1",sep = "")]) {
      Index2[i]<-ltable[i,paste(LD,"1",sep = "")]/(ltable[i,paste(LD,"1",sep = "")]+ltable[i,paste(LD,"2",sep = "")])
    }else{
      Index2[i]<-ltable[i,paste(LD,"2",sep = "")]/(ltable[i,paste(LD,"1",sep = "")]+ltable[i,paste(LD,"2",sep = "")])
    }
  }else if (is.null(WID)){
    if (ltable[i,paste(HT,"1",sep = "")]==ltable[i,paste(MT,"1",sep = "")]){
      Index1[i]<-ltable[i,paste(HD,"1",sep = "")]/(ltable[i,paste(HD,"1",sep = "")]+ltable[i,paste(HD,"2",sep = "")])
    }else{
      Index1[i]<-ltable[i,paste(HD,"2",sep = "")]/(ltable[i,paste(HD,"1",sep = "")]+ltable[i,paste(HD,"2",sep = "")])
    }
    if (ltable[i,paste(LT,"1",sep = "")]==ltable[i,paste(WT,"1",sep = "")]) {
      Index2[i]<-ltable[i,paste(LD,"1",sep = "")]/(ltable[i,paste(LD,"1",sep = "")]+ltable[i,paste(LD,"2",sep = "")])
    }else{
      Index2[i]<-ltable[i,paste(LD,"2",sep = "")]/(ltable[i,paste(LD,"1",sep = "")]+ltable[i,paste(LD,"2",sep = "")])
    }
  }else if (is.null(WID)&&is.null(PID)){
    if (ltable[i,paste(HT,"1",sep = "")]==ltable[i,Ref]){
      Index1[i]<-ltable[i,paste(HD,"2",sep = "")]/(ltable[i,paste(HD,"1",sep = "")]+ltable[i,paste(HD,"2",sep = "")])
    }else{
      Index1[i]<-ltable[i,paste(HD,"1",sep = "")]/(ltable[i,paste(HD,"1",sep = "")]+ltable[i,paste(HD,"2",sep = "")])
    }
    if (ltable[i,paste(LT,"1",sep = "")]==ltable[i,Ref]) {
      Index2[i]<-ltable[i,paste(LD,"2",sep = "")]/(ltable[i,paste(LD,"1",sep = "")]+ltable[i,paste(LD,"2",sep = "")])
    }else{
      Index2[i]<-ltable[i,paste(LD,"1",sep = "")]/(ltable[i,paste(LD,"1",sep = "")]+ltable[i,paste(LD,"2",sep = "")])
    }
  }

}

if (popt=="F2" || is.null(WID) || is.null(MID)){
  deltaSNP<-Index1-Index2
}else{
  deltaSNP<-abs(Index1-Index2)
}

REF_FRQ = (ltable[,paste(HD,"1",sep = "")] + ltable[,paste(LD,"1",sep = "")]) / (ltable[,paste(HighBulk,".DP",sep="")] + ltable[,paste(LowBulk,".DP",sep="")])

ntable<-data.frame(cbind(table[,c(1:5)],ltable[,paste(HD,"1",sep ="")],
                   ltable[,paste(HD,"2",sep ="")],ltable[,paste(HighBulk,".DP",sep ="")],
                   ltable[,paste(LD,"1",sep ="")],ltable[,paste(LD,"2",sep ="")],
                   ltable[,paste(LowBulk,".DP",sep ="")],Index1,Index2,deltaSNP,REF_FRQ),stringsAsFactors=F)
colnames(ntable)<-c("CHROM","POS","REF","ALT","Vtype","AD_REF.HIGH","AD_ALT.HIGH","DP.HIGH",
                    "AD_REF.LOW","AD_ALT.LOW","DP.LOW","SNPindex.HIGH","SNPindex.LOW",
                    "deltaSNP","REF_FRQ")

n<-dim(ntable)[1];ED<-NULL

for (i in 1:n) {
  ED[i]<-sqrt(((ltable[i,paste(HD,"1",sep = "")]/ltable[i,paste(HighBulk,".DP",sep="")])-(ltable[i,paste(LD,"1",sep = "")]/ltable[i,paste(LowBulk,".DP",sep="")])) ^ 2+
                ((ltable[i,paste(HD,"2",sep = "")]/ltable[i,paste(HighBulk,".DP",sep="")])-(ltable[i,paste(LD,"2",sep = "")]/ltable[i,paste(LowBulk,".DP",sep="")]))^2)
}

ntable$ED<-ED
#setwd(file.path(opt$output)
save.image("step01.RData")
nws<-as.numeric(opt$ws)
##the number of SNPs
nnSNP<-list()

chr<-unique(ntable$CHROM)

for (k in 1:length(chr)){
	nSNPs<-NULL
	int<-round((max(ntable[which(ntable[,1] %in% chr[k] ),2])-min(ntable[which(ntable[,1] %in% chr[k] ),2]))/nws)+1
	pos2<-NULL
	pos2[1]<-ntable[which(ntable[,1] %in% chr[k] )[1],2]+nws
	for (i in 2:int){
		pos2[i]<-pos2[i-1]+nws
	}
	pos1<-matrix(0,int,1)
	pos1[1]=ntable[which(ntable[,1] %in% chr[k] )[1],2]
	pos1[c(2:int)]=pos2[c(1:int-1)]
	for (j in 1:int){
		nsnp<-length(which(ntable[which(ntable[,1] %in% chr[k] ),2]>=pos1[j] & ntable[which(ntable[,1] %in% chr[k] ),2]<pos2[j]))
		nSNPs[which(ntable[which(ntable[,1] %in% chr[k] ),2]>=pos1[j] & ntable[which(ntable[,1] %in% chr[k] ),2]<pos2[j])]<-nsnp
	}

	nnSNP[[k]]<-as.data.frame(nSNPs)
}
nSNPs<-do.call(rbind,unname(nnSNP))$nSNPs
ntable$nSNPs<-do.call(rbind,unname(nnSNP))$nSNPs
##stats
chr<-unique(ntable$CHROM)
SNP<-NULL;INDEL<-NULL;
for (i in 1:length(chr)){

	SNP[i]<-length(which(ntable[which(ntable[,1] %in% chr[i]),5] %in% "SNP" ))
	INDEL[i]<-length(which(ntable[which(ntable[,1] %in% chr[i]),5] %in% "INDEL" ))

}

stat<-cbind(SNP,INDEL)
colnames(stat)<-c("SNP_Number","INDEL_Number")
rownames(stat)<-chr
write.csv(stat,file="stat.csv",quote=F,row.names=T)

#ED_INDEX<-sqrt((seq_G$SNPindex.LOW-seq_G$SNPindex.HIGH)^2)
EDpower<-ED^4
#setwd(paste(opt$output,"RData",sep=""))
#save.image(paste(opt$output,"RData/step01.RData",sep=""))
setwd(file.path(opt$output))
save.image("step01.RData")

escaptime<-Sys.time()-times;
print("Done!");
print(escaptime)



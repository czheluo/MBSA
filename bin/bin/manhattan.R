library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'input','i',0,'character',
	'output','o',0,'character',
	'thre','t',1,'character',
	'help','m',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--input	 the input  file
	--output	the out file 
	--thre    threshold default 0.9995
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$input)) { print_usage(spec)}
if ( is.null(opt$output)){ print_usage(spec) }
if ( is.null(opt$thre)){ thre = 0.9995 }else{thre<-as.numeric(opt$thre)}
times<-Sys.time()

bsa<-read.table(opt$input,head=TRUE)
library(qqman)
chr<-bsa$CHROM
pos<-bsa$POS
pos<-pos/1000000
chrlab=unique(chr)
thres<-quantile(bsa$Gprime,probs=thre,na.rm=TRUE)
df<-data.frame(chr=as.numeric(chr),pos=as.numeric(pos),Gprime_Log10Pval=as.numeric(bsa$Gprime),threhold=rep(thres,length(pos)))
df2<-data.frame(chr=as.numeric(chr),pos=as.numeric(bsa$POS),Gprime_Log10Pval=as.numeric(bsa$Gprime),threhold=rep(thres,length(pos)))
write.table(file=paste(opt$output,".select",sep=""),df2,sep="\t",row.names=FALSE)
write.table(file=paste(opt$output,".threshold.select",sep=""),subset(df2,df2$Gprime_Log10Pval > df2$threhold),row.names=FALSE)
pdf(paste(opt$output,"pdf",sep="."),height=9,width=16)
manhattan(df,chr="chr",bp="pos",p="Gprime_Log10Pval",col=rainbow(4),logp=FALSE,chrlabs=chrlab,ylab="Log10Pval",suggestiveline=thres,genomewideline = F,cex=2)
dev.off()
png(paste(opt$output,"png",sep="."),height=900,width=1600)
manhattan(df,chr="chr",bp="pos",p="Gprime_Log10Pval",col=rainbow(4),logp=FALSE,chrlabs=chrlab,ylab="Log10Pval",suggestiveline=thres,genomewideline = F,cex=2)
dev.off()


escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
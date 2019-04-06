#!/usr/bin/env Rscript
times<-Sys.time()
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
  'input','i',0,'character',
  'outdir','o',0,'character',
  'help','h',0,'logical'
), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
      Usage:
      --input	the input file name
      --outdir	the the out director
      --help		usage
      \n")
  q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$input)) { print_usage(spec)}
if ( is.null(opt$outdir)){ print_usage(spec) }


library(ggplot2)
library(magrittr)
library(dplyr)
library(tidyr)

format_genomic <- function(...) {
      function(x) {
            limits <- c(1e0,   1e3, 1e6)
            #prefix <- c("","Kb","Mb")
            # Vector with array indices according to position in intervals
            i <- findInterval(abs(x), limits)
            # Set prefix to " " for very small values < 1e-24
            i <- ifelse(i==0, which(limits == 1e0), i)
            paste(format(round(x/limits[i], 1),
                         trim=TRUE, scientific=FALSE, ...)
                #  ,prefix[i]
            )
      }
}

SNPset<-read.table("result/three.methods.result",header=T)

CHR<-unique(SNPset$Chr)

SNPset$Chr<-factor(SNPset$Chr,levels =paste(unique(SNPset$Chr)))


var<-"tricubeDeltaSNP"

library(RColorBrewer)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
rcolor<-color[sample(1:length(color),length(color))]


bar<-ggplot(data=SNPset,aes(x=Pos,y=ED4_fit,colour=factor(Chr)))+
	scale_x_continuous(breaks = seq(from = 0,to = max(SNPset$Pos),
                                  by = 10^(floor(log10(max(SNPset$Pos))))))+
				   theme(plot.margin =margin(b = 0,l = 20,r = 20,unit = "pt"),
		         axis.title.x=element_blank(),
		         axis.text.x=element_blank(),
		         legend.position="none",#panel.background=element_blank(),
		         panel.border=element_blank(),
		         axis.ticks=element_blank(),
		          strip.text = element_text(colour = "green3"),
			  strip.text.x = element_text(size = 8,face = "bold"),
			  strip.background = element_rect(fill="dodgerblue3"))#,#panel.grid.major=element_blank(),
		         #panel.grid.minor=element_blank())#,plot.background=element_blank())


bn <- bar + facet_grid(~ Chr, scales = "free_x", space = "free_x")

times<-Sys.time()
ba<-bn+geom_bar(stat="identity")+scale_color_manual(breaks = paste(unique(SNPset$Chr)),
	values=rcolor[c(1:length(unique(SNPset$Chr)),length(unique(SNPset$Chr)))])##add one was to able plot the I95

ba
batime<-Sys.time()-times


var="WindowDeltaSNP"

pa <- ggplot(data = SNPset,aes(colour=factor(Chr))) +
  scale_x_continuous(breaks = seq(from = 0,to = max(SNPset$Pos),
                                  by = 10^(floor(log10(max(SNPset$Pos)))))) +
		   theme(plot.margin =margin(b = 0,l = 10,r = 5,unit = "pt"),
		         axis.title.x=element_blank(),
		         axis.text.x=element_blank(),
		         legend.position="none",#panel.background=element_blank(),
		         panel.border=element_blank(),
		         axis.ticks=element_blank(),
		          #strip.text = element_text(colour = "black"),
			   #strip.text.x = element_text(size = 3,face = "bold"),
			  strip.background = element_rect(fill="dodgerblue3"))#,#panel.grid.major=element_blank(),
		         #panel.grid.minor=element_blank())#,plot.background=element_blank())

pn <- pa + facet_grid(~ Chr, scales = "free_x", space = "free_x")

#scale_fill_manual(labels=unique(SNPset$Chr), values=str(rcolor[c(1:length(unique(SNPset$Chr)))]))

p1<-pn + geom_line(aes_string(x = "Pos", y = var),size=0.8 )+

scale_color_manual(breaks = factor(paste(unique(SNPset$Chr)),levels =paste(unique(SNPset$Chr))),
	values=rcolor[c(1:length(unique(SNPset$Chr)),length(unique(SNPset$Chr))+1,length(unique(SNPset$Chr))+2)])##add one was to able plot the I95

pwd<-p1+geom_line(data = SNPset,aes(x = Pos, y = I95 ,color="red"))#+geom_line(data = SNPset,aes(x = Pos, y =- I95 ,color="red"))


#PT <-select(SNPset, Chr, Pos,matches("I95")) %>% gather(key = "Interval", value = "value",-Chr,-Pos)

#pwd <- p1 + geom_line(data = PT,aes(x = Pos, y = value, color = Interval)) +
	#geom_line(data = PT,aes(x = Pos,y = -value,color = Interval))

library(grid)
g <- ggplot_gtable(ggplot_build(pwd))

stripr <- which(grepl('strip-t', g$layout$name))

fills <- rcolor[c(1:length(unique(SNPset$Chr)))]
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid.draw(g)



stripr <- which(grepl('axis-t', g$layout$name))

fills <- rcolor[c(1:length(unique(SNPset$Chr)))]
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}


grid.draw(g)






SNPset$color<-rc

var="ED4_fit"
pa <- ggplot(data = SNPset) +
  scale_x_continuous(breaks = seq(from = 0,to = max(SNPset$Pos),
                                  by = 10^(floor(log10(max(SNPset$Pos)))))) +
		   theme(plot.margin =margin(b = 0,l = 20,r = 20,unit = "pt"),
		         axis.title.x=element_blank(),
		         axis.text.x=element_blank(),
		         legend.position="none",#panel.background=element_blank(),
		         panel.border=element_blank(),
		         axis.ticks=element_blank())
		          #strip.text = element_text(colour = rcolor[c(1:length(CHR))]))#,#panel.grid.major=element_blank(),
		         #panel.grid.minor=element_blank())#,plot.background=element_blank())

pn <- pa + facet_grid(~ Chr, scales = "free_x", space = "free_x")

#  scale_fill_manual(labels=unique(SNPset$Chr), values=str(rcolor[c(1:length(unique(SNPset$Chr)))]))

p1<-pn + geom_line(aes_string(x = "Pos", y = var) )



cutoff <- 3*(sd(SNPset$ed_index)+median(SNPset$ed_index))

##Threshold
p1 <-p1+geom_hline(aes_string(yintercept = "threshold"),
    color = "red",
    size = 1,
    alpha = 0.4)


var="WindowDeltaSNP"

pb <- ggplot(data = SNPset,aes( colour = factor(Chr))) +scale_x_continuous(breaks = seq(from = 0,to = max(SNPset$Pos),by = 10^(floor(log10(max(SNPset$Pos))))),
                                                                          labels = format_genomic(), name = "Genomic Position (Mb)") +
  theme(plot.margin =margin(b = 2,l = 20,r = 20,unit = "pt"),#legend.position="none",
        strip.background = element_blank(),
        strip.text.x = element_blank())

pm <- pb + facet_grid(~ Chr, scales = "free_x", space = "free_x")

p2<-pm + geom_line(aes_string(x = "Pos", y = var) )


library(gridExtra)

grid.arrange(p1, p2, nrow = 2)

library(reshape2)
library(ggplot2)
library(grid)

png(paste("BSA",".png",sep=""),width=1600, height=900)

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
gA$widths <- gB$widths
grid.newpage()
grid.draw(arrangeGrob(gA,gB, heights = c(1/2,1/2)) )

dev.off()

p + geom_point(aes_string(x = "Pos", y = var))

var<-"Gprime"

p+geom_bar(x = "POS", y = var)


var<-"ED"
SNPset<-seq

var<-"ED4_fit"


p <- ggplot(data = SNPset) +scale_x_continuous(breaks = seq(from = 0,to = max(SNPset$POS),by = 10^(floor(log10(max(SNPset$POS))))),
                   labels = format_genomic(), name = "Genomic Position (Mb)") +
		   theme(plot.margin =margin(b = 10,l = 20,r = 20,unit = "pt"))



p <- p + facet_grid(~ CHROM, scales = "free_x", space = "free_x")


p + geom_line(aes_string(x = "POS", y = var) )


p + geom_point(aes_string(x = "POS", y = var))

p+geom_bar(x = "POS", y = var)


var<-"ED"
SNPset<-seq




escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
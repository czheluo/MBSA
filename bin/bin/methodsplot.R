
library(grid)
library(ggplot2)
library(magrittr)
library(dplyr)
library(tidyr)
library(gridExtra)

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
#SNPset<-result

#SNPset$Pos<-SNPset$Pos*1000000
#CHR<-unique(SNPset$Chr)

#SNPset$Chr<-factor(SNPset$Chr,levels =paste(unique(SNPset$Chr)))

CN<-rep(c(1:4),length(SNPset$Chr))
rcolor<-NULL
for (i in 1:length(unique(SNPset$Chr))) {

  rcolor[i]<-c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3")[CN[i]]

}


##barplot

bar<-ggplot(data=SNPset,aes(x=Pos,y=ED4_fit,colour=factor(Chr)))+
  scale_x_continuous(breaks = seq(from = 0,to = max(SNPset$Pos),
                                  by = 2*10^(floor(log10(max(SNPset$Pos))))))+
  theme(plot.margin =margin(b = 0,l = 20,r = 20,unit = "pt"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        legend.position="none",#panel.background=element_blank(),
        panel.border=element_blank(),
        axis.ticks=element_blank(),
        strip.text = element_text(colour = "black"),
        strip.text.x = element_text(size = 8,face = "bold"))
#strip.background = element_rect(fill="dodgerblue3"))#,#panel.grid.major=element_blank(),
#panel.grid.minor=element_blank())#,plot.background=element_blank())


bn <- bar + facet_grid(~ Chr, scales = "free_x", space = "free_x")

ba<-bn+geom_bar(stat="identity")+scale_color_manual(breaks = paste(unique(SNPset$Chr)),
   values=rcolor[c(1:length(unique(SNPset$Chr)),length(unique(SNPset$Chr)))])##add one was to able plot the I95
g <- ggplot_gtable(ggplot_build(ba))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- rcolor[c(1:length(unique(SNPset$Chr)))]
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
#grid.draw(g)


##line
##

pa <- ggplot(data = SNPset,aes(x=Pos,y=WindowDeltaSNP,colour=factor(Chr)))+#,aes(colour=factor(Chr))) +
  scale_x_continuous(breaks = seq(from = 0,to = max(SNPset$Pos),
                                  by = 2*10^(floor(log10(max(SNPset$Pos)))))) +
  theme(plot.margin =margin(b = 0,l = 10,r = 5,unit = "pt"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        legend.position="none",#panel.background=element_blank(),
        panel.border=element_blank(),
        axis.ticks=element_blank(),
        strip.text = element_text(colour = "black"),
        strip.text.x = element_text(size = 10,face = "bold"),
      # strip.background = element_rect(fill="dodgerblue3"))
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())#,plot.background=element_blank())
pn <- pa + facet_grid(~ Chr, scales = "free_x", space = "free_x")
#scale_fill_manual(labels=unique(SNPset$Chr), values=str(rcolor[c(1:length(unique(SNPset$Chr)))]))
p1<-pn + geom_line(size=1.4)+#aes_string(x = "Pos", y = var),
  scale_colour_manual(breaks = paste(unique(SNPset$Chr)),
                      values=rcolor[c(1:length(unique(SNPset$Chr)),
                      length(unique(SNPset$Chr))+1,length(unique(SNPset$Chr))+2)])##add one was to able plot the I95
da<-data.frame(SNPset$Chr,SNPset$Pos,SNPset$I95)
colnames(da)<-c("Chr","Pos","I95")
pwd<-p1+geom_line(data = da,aes(x = Pos, y = I95 ,color=Chr))+
  geom_line(data = da,aes(x = Pos, y =- I95 ,color=Chr))+
  labs(y="SNP-index")

g1 <- ggplot_gtable(ggplot_build(pwd))

stripr <- which(grepl('strip-t', g1$layout$name))

fills <- rcolor[c(1:length(unique(SNPset$Chr)))]
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
  g1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
#grid.draw(g1)


ed <- ggplot(data = SNPset,aes(x=Pos,y=ED4_fit,colour=factor(Chr)))+#,aes(colour=factor(Chr))) +
  scale_x_continuous(breaks = seq(from = 0,to = max(SNPset$Pos),
                                  by = 10^(floor(log10(max(SNPset$Pos)))))) +
  theme(plot.margin =margin(b = 0,l = 10,r = 5,unit = "pt"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        legend.position="none",#panel.background=element_blank(),
        panel.border=element_blank(),
        axis.ticks=element_blank(),
        strip.text = element_text(),
        #strip.text.x = element_text(),
        #strip.background = element_rect(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())#,plot.background=element_blank())
ped1 <- ed + facet_grid(~ Chr, scales = "free_x", space = "free_x")
#scale_fill_manual(labels=unique(SNPset$Chr), values=str(rcolor[c(1:length(unique(SNPset$Chr)))]))
ped2<-ped1 + geom_line(size=1.4)+#aes_string(x = "Pos", y = var),
  scale_colour_manual(breaks = factor(paste(unique(SNPset$Chr)),levels =paste(unique(SNPset$Chr))),
                      values=rcolor[c(1:length(unique(SNPset$Chr)),length(unique(SNPset$Chr))+1,length(unique(SNPset$Chr))+2)])##add one was to able plot the I95

da<-data.frame(SNPset$Chr,SNPset$Pos,SNPset$ED4_fit_cutoff)
colnames(da)<-c("Chr","Pos","cutoff")

ped3<-ped2+geom_line(data = da,aes(x = Pos, y = cutoff ,color=Chr),size=1)+
  labs(y="ED4 losse fit")


g2 <- ggplot_gtable(ggplot_build(ped3))

stripr1 <- which(grepl('strip-t', g2$layout$name))

fills <- rcolor[c(1:length(unique(SNPset$Chr)))]
k <- 1
for (i in stripr1) {
  j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
#grid.draw(g2)


SNPset$gp<-c(matrix(min(SNPset$Gprime[SNPset$Gqvalue<0.01]),length(SNPset$Chr)))

Gp <- ggplot(data = SNPset,aes(x=Pos,y=Gprime,colour=factor(Chr)))+#,aes(colour=factor(Chr))) +
  scale_x_continuous(breaks = seq(from = 0,to = max(SNPset$Pos),
                                  by = 10^(floor(log10(max(SNPset$Pos)))))) +
  theme(plot.margin =margin(b = 0,l = 10,r = 5,unit = "pt"),
        #axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        legend.position="none",#panel.background=element_blank(),
        panel.border=element_blank(),
       # axis.ticks=element_blank(),
        #strip.text = element_text(),
        #strip.text.x = element_text(),
        #strip.background = element_rect(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),panel.grid.major=element_blank(),
panel.grid.minor=element_blank())#,plot.background=element_blank())
Gp1 <- Gp + facet_grid(~ Chr, scales = "free_x", space = "free_x")
#scale_fill_manual(labels=unique(SNPset$Chr), values=str(rcolor[c(1:length(unique(SNPset$Chr)))]))
Gp2<-Gp1 + geom_line(size=1.4)+#aes_string(x = "Pos", y = var),
  scale_colour_manual(breaks = paste(unique(SNPset$Chr)),
  values=rcolor[c(1:length(unique(SNPset$Chr)),
  length(unique(SNPset$Chr))+1,length(unique(SNPset$Chr))+2)])##add one was to able plot the I95

da<-data.frame(SNPset$Chr,SNPset$Pos,SNPset$gp)
colnames(da)<-c("Chr","Pos","gp")
#da<-dplyr::select(SNPset, Chr, Pos, dplyr::matches("gp")) %>%
  #tidyr::gather(key = "Interval", value = "value",-Chr,-Pos)
Gp3<-Gp2+geom_line(data = da,aes(x = Pos, y = gp ,color=Chr),
                   size=1.2)+labs(x="Genomic Position (Mb)",y="G' Value")

g3 <- ggplot_gtable(ggplot_build(Gp3))

stripr <- which(grepl('strip-t', g3$layout$name))

fills <- rcolor[c(1:length(unique(SNPset$Chr)))]
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g3$grobs[[i]]$grobs[[1]]$childrenOrder))
  g3$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
#grid.draw(g3)

pdf(paste("three.methods","pdf",sep="."),height=9,width=16)

grid.arrange(g1, g2, g3,nrow = 3)

dev.off()

png(paste("three.methods","png",sep="."),height=900,width=1600)

grid.arrange(g1, g2, g3,nrow = 3)

dev.off()















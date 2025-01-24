library(ggplot2)
library(plyr)
library(reshape2)
library(zoo)
library(pheatmap)
#library(HiCcompare)
library(Matrix)
library(data.table)
################################################################################################  20221116

sca_raw<-fread('hichew/gw/normedAPA_reformat.txt')
row.names(sca_raw)<-10*(as.integer(row.names(sca_raw))-11)

pdf('apa_snHiChew.pdf',height=6,width=6)
pheatmap(color=colorRampPalette(c("white",'red'))(50),sca_raw,border_color=NA,cluster_rows=F,cluster_cols=F,
         show_rownames=T,show_colnames=F,breaks=seq(0.5,2.5,0.06),main = "snHiChew, 10k bin, normedAPA\n950 loops called from ngsHiC, \nPeak to Lower Left, 6.152") 
dev.off()


#library(gridExtra)
#grid.arrange(p1,p2,nrow=1,ncol=2)

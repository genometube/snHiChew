library(ggplot2)
library(plyr)
library(reshape2)
library(zoo)
library(pheatmap)
#library(HiCcompare)
library(Matrix)
library(data.table)
################################################################################################  20221116

sca_raw<-fread('bulkHiChew20241129/gw/normedAPA_.txt')
row.names(sca_raw)<-10*(as.integer(row.names(sca_raw))-11)

pdf('apa_hichew.pdf',height=6,width=6)
pheatmap(color=colorRampPalette(c("white",'red'))(50),sca_raw,border_color=NA,cluster_rows=F,cluster_cols=F,
         show_rownames=T,show_colnames=F,breaks=seq(0.5,2.5,0.06),main = "bulkHiChew, 10k bin, normedAPA\n950 loops called from ngsHiC, \nPeak to Lower Left, 5.018") 
dev.off()


ngs_raw<-fread('hek_ngsHiC_gw/normedAPA_.txt')
row.names(ngs_raw)<-10*(as.integer(row.names(ngs_raw))-11)

pdf('apa_ngs.pdf',height=6,width=6)
pheatmap(color=colorRampPalette(c("white",'red'))(50),ngs_raw,border_color=NA,cluster_rows=F,cluster_cols=F,
         show_rownames=F,show_colnames=F,breaks=seq(0.5,2.5,0.06),main = "ngsHiC, 10k bin, normedAPA\n950 loops called, \nPeak to Lower Left, 4.413") 
dev.off()

# Load libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library("dendextend") 
library(ggsci)
library(eulerr)
# Load libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library("dendextend") 
library(data.table)
library(IRanges)

#sample_group<-c('62500','125000','250000','500000')
#chr7    145813453       148118088       CNTNAP2
plot_list <- list()
#145,813,893-148,118,090
#chr7:3086310-3703105
#chr7    71244476        71877360        CALN1
#chr7    151163098       151217010       RHEB
start=145000000
end=155000000
#end-start
sample_group<-c('500000')
sample='hichew_293t'
bin=25000
clust_num=2
clust_method="ward.D2"
chrom="chr7"
validpair="500000"
unit=1000000


############################################################################################ 20240427
ins_input<-read.table(paste0(sample,'_',validpair,'_ins_cell_sum_25000.xls.gz'),sep = '\t')
colnames(ins_input)<-c('pos','ins','cell')

# add repli-seq annotation
cell_early_id<-read.table(paste0('snHiChew_HEK293wt_label_info.txt'),sep = '\t')
colnames(cell_early_id)<-c('cell','ID','early_ratio')
head(cell_early_id)

ins_input_subset<-ins_input[ins_input$pos>=floor(start/bin)&ins_input$pos<=floor(end/bin),]

ins_merge<-merge(ins_input_subset,cell_early_id,by='cell')
ins_raw <- ins_merge[order(ins_merge$early_ratio), ]
thres=1
ins_raw$ins[ins_raw$ins< -thres] <- -thres
ins_raw$ins[ins_raw$ins>thres] <- thres
#ggplot(ins_raw,aes(x=ins))+geom_density()

dat_mtx<-reshape2::dcast(ins_raw,cell~pos,value.var = "ins",fun.aggregate = sum)
dat_mtx_reorder <- dat_mtx[match(unique(ins_raw$cell), dat_mtx$cell), ]

dat_mtx_<-dat_mtx[,-1]
rownames(dat_mtx_)<-dat_mtx_reorder$cell
ins_raw_subset <- ins_raw[match(row.names(dat_mtx_), ins_raw$cell), ]
annot_row <- data.frame(row.names = ins_raw_subset$cell, early_ratio = ins_raw_subset$early_ratio)

ann_colors = list(early_ratio = c("yellow", "firebrick"))
dat_mtx__ins<-as.matrix(dat_mtx_)

cell_clust<-as.data.frame(cutree(hclust(dist(scale(dat_mtx__ins),method='euclidean'),method  = clust_method), k = clust_num))
cell_clust$cell<-rownames(cell_clust)
colnames(cell_clust)<-c('cell_cluster','cell')
cell_clust$group<-cell_clust$cell_cluster
cell_clust_1<-cell_clust[cell_clust$cell_cluster==1,2]
cell_clust_2<-cell_clust[cell_clust$cell_cluster==2,2]


############# 20240911 ############# 
cellID_dist_contact_raw<-read.table('cellID_dist_pct.xls.gz')
#head(cellID_dist_contact)
colnames(cellID_dist_contact_raw)<-c('cell','dist','contact')

#ggplot(cellID_dist_contact, aes(x=contact)) + geom_density()+xlim(0,0.0005)
cellID_dist_contact<-cellID_dist_contact_raw

dat_mtx<-reshape2::dcast(cellID_dist_contact,dist~cell,value.var = "contact",fun.aggregate = sum)
dat_mtx_df<-as.data.frame(dat_mtx)
dat_mtx_df_<-dat_mtx_df[c(-1:-20),-1]

dat_mtx_<-as.matrix(dat_mtx_df_)

length(dat_mtx[c(-1:-20),1])
length(10^(seq(log10(1000), log10(max(80000000)), by = 0.05))[-1])
10^(seq(log10(1000), log10(max(80000000)), by = 0.05))[-1:-21]

rownames(dat_mtx_)<-dat_mtx[c(-1:-20),1]
rownames(dat_mtx_)<-round(10^(seq(log10(1000), log10(max(80000000)), by = 0.05))[-1:-21],0)

dat_mtx_[dat_mtx_>2]<-2
dat_mtx_[dat_mtx_<0.5]<-0.5

cellID_barcode_repliscore<-read.table('cellID_barcode_repliscore.xls.gz')

annot_col <- data.frame(early_ratio = cellID_barcode_repliscore$V3,cell_cluster=cell_clust$group)

ann_colors = list(early_ratio = c("yellow", "firebrick"))
ann_colors = list(early_ratio = c("yellow", "firebrick"),cell_cluster=c('1'='royalblue','2'='#E18B41'))
#ann_colors = list(early_ratio = c("yellow", "firebrick"),cell_cluster=c('1'='olivedrab','2'='#B571DC'))

ann_colors = list(early_ratio = c('brown4',"#E18B41",'white', "royalblue"),cell_cluster=c('1'='#009681','2'='#CC476B'))
pdf('repli_score_contact_distance_hm.pdf',height=4,width=15)
library(RColorBrewer)
pheatmap(color=colorRampPalette(brewer.pal(11, "Spectral"))(50),
         dat_mtx_[nrow(dat_mtx_):1, ], annotation_col   = annot_col, cluster_rows=F,cluster_cols=F,show_rownames=T,show_colnames=F,
         annotation_colors = ann_colors,
         annotation_names_col=F,border_color=NA,
         main = paste0('250k validpair per cell ',length(rownames(annot_col)),' cells\n',
                       'y: genomic distance, x: repliscore ranked cell, hm value: %contact'))
dev.off()

png('repli_score_contact_distance_hm.png',height=400,width=1500)
pheatmap(color=colorRampPalette(brewer.pal(11, "Spectral"))(50),
         dat_mtx_[nrow(dat_mtx_):1, ], annotation_col   = annot_col, cluster_rows=F,cluster_cols=F,show_rownames=T,show_colnames=F,
         annotation_colors = ann_colors,
         annotation_names_col=F,border_color=NA,
         main = paste0('250k validpair per cell ',length(rownames(annot_col)),' cells\n',
                       'y: genomic distance, x: repliscore ranked cell, hm value: %contact'))
dev.off()


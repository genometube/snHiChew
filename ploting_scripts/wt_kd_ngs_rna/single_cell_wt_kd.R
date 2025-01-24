# Load libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library("dendextend") 
library(ggsci)
library(eulerr)
mypal <- pal_npg("nrc", alpha = 0.7)(9)

library(data.table)
wt_kd_sc_stat<-fread('final.filter.proportion.xls')

# Create a sample data frame
data <- data.frame( category = c("WT_CNTNAP2_expressed", "WT_CNTNAP2_not_expressed"),
                    count = c(1236,31328))

# Calculate the percentage for each category
data$percentage <- data$count / sum(data$count) * 100

p0<-ggplot(data, aes(x = "", y = percentage, fill = category)) +
  geom_bar(width = 1, stat = "identity",alpha=0.7) +
  coord_polar("y", start = 0) +
  theme_void() + # Remove background, grid, and numeric labels
  labs(title = paste0('  WT_CNTNAP2_expressed:',1236,' WT_CNTNAP2_not_expressed:',31328)) +
  theme(legend.title = element_blank())+scale_fill_manual(values=c(mypal[1],mypal[2]))
p0

# Create a sample data frame
data <- data.frame( category = c("CTCFkd_CNTNAP2_expressed", "CTCFkd_CNTNAP2_not_expressed"),
                    count = c(3626,40273))

# Calculate the percentage for each category
data$percentage <- data$count / sum(data$count) * 100

p1<-ggplot(data, aes(x = "", y = percentage, fill = category)) +
  geom_bar(width = 1, stat = "identity",alpha=0.7) +
  coord_polar("y", start = 0) +
  theme_void() + # Remove background, grid, and numeric labels
  labs(title = paste0('  WT_CNTNAP2_expressed:',3626,' WT_CNTNAP2_not_expressed:',40273)) +
  theme(legend.title = element_blank())+scale_fill_manual(values=c(mypal[1],mypal[2]))
p1

library(gridExtra)
pdf('scRNA_wt_CTCFkd.pdf',width = 9,height=9)
grid.arrange(p0,p1,nrow=2)
dev.off()


################################################################################
rna<-read.table('hg19_rna_seq.bed.gz',sep = '\t',header = T)
long_gene<-rna[(rna$end-rna$start)>300000,]
colnames(long_gene)<-c('chr','start','end','transcript','type','signal','gene','info')

melt_sum<-read.table('tx_melt_score.xls')
colnames(melt_sum)<-c('transcript','ins_mean_clust_1','ins_mean_clust_2','d','pval','cellNum_clust_1','cellNum_clust_2')

melt_sum$clust_2_clust_1_ratio<-melt_sum$cellNum_clust_1/(melt_sum$cellNum_clust_1+melt_sum$cellNum_clust_2)
#ggplot(melt_sum, aes(x = log2(melt_sum$clust_2_clust_1_ratio))) +  geom_density()
#ggplot(melt_sum, aes(x = (melt_sum$clust_2_clust_1_ratio))) +  geom_density()

melt_sum_exp<-merge(melt_sum,long_gene,by='transcript')
melt_sum_exp$exp_norm_by_length<-melt_sum_exp$signal/(melt_sum_exp$end-melt_sum_exp$start)
melt_sum_exp<-melt_sum_exp[melt_sum_exp$d<1,]

ggplot(melt_sum_exp, aes(x = melt_sum_exp$signal)) +  geom_density()

#filtering a/b
library(data.table)
eigen<-fread('snHiChew_HEK293Twt.mcool.res_500000.eig.cis.vecs.tsv')
eigen_b<-eigen[eigen$E1<0,]
colnames(eigen_b)[1]<-'chr'
eigen_b_df<-as.data.frame(eigen_b)

is_within_range <- function(chr, start, end, table2) {
  chr_table2 <- table2 %>% filter(chr == !!chr)
  any(apply(chr_table2, 1, function(range_row) {
    (end >= as.numeric(range_row['start']) && end <= as.numeric(range_row['end'])) ||
      (start >= as.numeric(range_row['start']) && start <= as.numeric(range_row['end']))  }))}

# Apply the filter
filtered_melt_sum_exp <- melt_sum_exp %>%  rowwise() %>%  filter(is_within_range(chr, start, end, eigen_b))
b_comp_high<-length(rownames(filtered_melt_sum_exp[filtered_melt_sum_exp$signal>9,]))
b_comp_low<-length(rownames(filtered_melt_sum_exp[filtered_melt_sum_exp$signal<=9,]))

melt_sum_exp_a<-melt_sum_exp[!melt_sum_exp$transcript %in% filtered_melt_sum_exp$transcript,]
a_comp_high<-length(rownames(melt_sum_exp_a[melt_sum_exp_a$signal>9,]))
a_comp_low<-length(rownames(melt_sum_exp_a[melt_sum_exp_a$signal<=9,]))

# Create a sample data frame
data <- data.frame( category = c("b_comp_high", "b_comp_low"),
  count = c(b_comp_high,b_comp_low))

# Calculate the percentage for each category
data$percentage <- data$count / sum(data$count) * 100

# Create the pie chart
#p0<-ggplot(data, aes(x = "", y = percentage, fill = category)) +
#  geom_bar(width = 1, stat = "identity",alpha=0.7) +
#  coord_polar("y", start = 0) +
#  theme_void() + # Remove background, grid, and numeric labels
#  labs(title = paste0('  b_comp_high:',b_comp_high,' b_comp_low:',b_comp_low,' \n  a_comp_high:',a_comp_high,' a_comp_low:',a_comp_low)) +
#  theme(legend.title = element_blank())+scale_fill_manual(values=c(mypal[2],mypal[4],mypal[1],mypal[3]))

p0<-ggplot(data, aes(x = "", y = percentage, fill = category)) +
  geom_bar(width = 1, stat = "identity",alpha=0.7) +
  coord_polar("y", start = 0) +
  theme_void() + # Remove background, grid, and numeric labels
  labs(title = paste0('  b_comp_high:',b_comp_high,' b_comp_low:',b_comp_low)) +
  theme(legend.title = element_blank())+scale_fill_manual(values=c(mypal[1],mypal[3]))

p0
filtered_melt_sum_exp$exp_level<-'transcript_low_exp'
filtered_melt_sum_exp[filtered_melt_sum_exp$signal>9,length(colnames(filtered_melt_sum_exp))] <-'transcript_high_exp'

melt_sum_exp_high<-filtered_melt_sum_exp[filtered_melt_sum_exp$signal>9,]

write.table(melt_sum_exp_high, file = 'melt_sum_exp_high.xls',sep = "\t",col.names=T,row.names=FALSE,quote=FALSE)


#high
pcc<-round(cor(log2(melt_sum_exp_high$exp_norm_by_length),melt_sum_exp_high$d),3)
geneNum<-length(melt_sum_exp_high$exp_norm_by_length)
p1<-ggplot(melt_sum_exp_high, aes(x = melt_sum_exp_high$d,y=log2(melt_sum_exp_high$exp_norm_by_length))) +  
  geom_point(alpha=0.5,size=3,stroke=0,color='red3')+theme_classic()+geom_smooth(method='lm',color='grey20',linetype='dashed')+
  labs(title=paste0('PCC:',pcc,', transcript Num:',geneNum),subtitle = 'high exp long transcript (normalized RPM > 9)',
       x='KS test maximum distance',y='log2 length scaled transcript RPM')
p1

p2<-ggplot(filtered_melt_sum_exp, aes(x = d,group=exp_level,color=exp_level)) +  
  geom_density(lwd=1)+theme_classic()+scale_color_manual(values=c(mypal[1],mypal[3]))+
  labs(title=paste0(subtitle = 'high exp long transcript (normalized RPM > 9)'),x='KS test maximum distance',y='density')
p2

#low
melt_sum_exp_low<-filtered_melt_sum_exp[filtered_melt_sum_exp$signal<9,]
pcc<-round(cor(log2(melt_sum_exp_low$exp_norm_by_length+0.00001),melt_sum_exp_low$d),3)
geneNum<-length(melt_sum_exp_low$exp_norm_by_length)
p3<-ggplot(melt_sum_exp_low, aes(x = melt_sum_exp_low$d,y=log2(melt_sum_exp_low$exp_norm_by_length))) +  
  geom_point(alpha=0.5,size=3,stroke=0,color=mypal[3])+theme_classic()+geom_smooth(method='lm',color='grey20',linetype='dashed')+
  labs(title=paste0('PCC:',pcc,', transcript Num:',geneNum),subtitle = 'low exp long transcript (normalized RPM < 9)',
       x='KS test maximum distance',y='log2 length scaled transcript RPM')
p3




library(gridExtra)

svg('melting_b.svg',width = 9,height=9)
grid.arrange(p0,p2,p1,p3,nrow=2)
dev.off()

pdf('melting_b.pdf',width = 9,height=9)
grid.arrange(p0,p2,p1,p3,nrow=2)
dev.off()



##################################################

bin=25000
unit=1000000
validpair='500000'
clust_num=2
clust_method="ward.D2"
ins<-read.table('gene_1_ins.xls',sep = '\t')
colnames(ins)<-c('chr','start','end','ins','cell','gene')
chrom=ins[1,1]
ins$bin_id<-ins$start/bin
dat_mtx<-dcast(ins,cell~bin_id,value.var = "ins",fun.aggregate = sum)

#dat_mtx_reorder <- dat_mtx[match(unique(ins$cell), dat_mtx$cell), ]

dat_mtx_<-dat_mtx[,-1]
rownames(dat_mtx_)<-dat_mtx$cell
dat_mtx__ins<-as.matrix(dat_mtx_)

cell_clust<-as.data.frame(cutree(hclust(dist(scale(dat_mtx__ins)),method  = clust_method), k = clust_num))
cell_clust$cell<-rownames(cell_clust)
colnames(cell_clust)<-c('cell_cluster','cell')

dend = as.dendrogram(hclust(dist(scale(dat_mtx__ins),method='euclidean'),method  = clust_method))
#dend = as.dendrogram(hclust(dist(dat_mtx__ins),method  = clust_method))
clust_num=2
dend = color_branches(dend, k = clust_num)

base_mean = rowMeans(dat_mtx_)
p<-grid.grabExpr(draw(Heatmap(dat_mtx__ins, name = "ins", cluster_rows = dend, row_split = clust_num,
                              cluster_columns = FALSE,
                              show_row_names = FALSE,show_column_names = FALSE,
                              column_title=paste0('cluster_1: ',length(cell_clust[cell_clust$group==1,3]),' cluster_2: ',length(cell_clust[cell_clust$group==2,3]),'\n',
                                                  chrom,': ',start,' - ',end,'\nvalidpair downsample: ',
                                                  validpair, '\nclust method: ',clust_method,'\ncell_num: ',
                                                  length(rownames(dat_mtx_)))),heatmap_legend_side = "top"))  

Heatmap(dat_mtx__ins,cluster_rows = dend,row_split = clust_num,cluster_columns = FALSE,show_row_names = FALSE,show_column_names = FALSE,
        col=colorRampPalette(c("darkgreen",'white','darkred'))(50))

cell_clust_1<-cell_clust[cell_clust$cell_cluster==1,2]
cell_clust_2<-cell_clust[cell_clust$cell_cluster==2,2]

head(cell_clust)
if (mean(ins[ins$cell%in%cell_clust_1,4]) < mean(ins[ins$cell%in%cell_clust_2,4])){
  cell_clust$cell_cluster <- ifelse(cell_clust$cell_cluster == 1, 2, ifelse(cell_clust$cell_cluster == 2, 1, cell_clust$cell_cluster)) 
  cell_clust_1<-cell_clust[cell_clust$cell_cluster==1,2]
  cell_clust_2<-cell_clust[cell_clust$cell_cluster==2,2]
  }

ins_clust_1<-ins[ins$cell%in%cell_clust_1,]
ins_clust_1$clust<-'cluster_1'
ins_clust_2<-ins[ins$cell%in%cell_clust_2,]
ins_clust_2$clust<-'cluster_2'
ins_clust<-rbind.data.frame(ins_clust_1,ins_clust_2)

ins_mean_clust_1<-mean(ins[ins$cell%in%cell_clust_1,4])
ins_mean_clust_2<-mean(ins[ins$cell%in%cell_clust_2,4])

ggplot(ins_clust, aes(x = ins,group=clust,color=clust)) +  stat_ecdf(geom = "step") +
  labs(title = "Cumulative Distribution Function (CDF)", x = "ins", y = "ECDF") +  theme_classic()
ks_test<-ks.test(ins_clust_2$ins,ins_clust_1$ins,alternative = c("greater"))
ks.test(ins_clust_2$ins,ins_clust_1$ins,alternative = c("greater"))
d<-ks_test[1]
pval<-ks_test[2]

ins_metrics<-data.frame(gene=ins[1,6],ins_mean_clust_1=ins_mean_clust_1,ins_mean_clust_2=ins_mean_clust_2,d=d,pval=pval)
ins_metrics
long_gene
merge(long_gene,ins_metrics,by='gene')
write.table(HiC_merge, file = out_path_to_file,sep = "\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

#########################

chrom="chr7"
#start_list<-c(1,14,69,75,145,155,145)
#end_list<-c(5,20,73,80,155,159,159)
start_list<-c(1,65,145)
end_list<-c(15,75,155)
start_list<-start_list*unit
end_list<-end_list*unit
sample='hichew_293t'
plot_list <- list()
#start=145000000
#end=155000000
metrics<-data.frame()
filtered_metrics<-read.table('C:/Users/xieyeming1/Insync/yemingxie@gmail.com/Google\ Drive/work/active_project_2024/hichew_paper_2024/hichew_NM_resubmit_figs/single_cell_compare/snHiChew_HEK293Twt/sorted_snHiChew_HEK293Twt_cut.xls',header=T)

#i=1
for(i in 1:length(start_list)){
  start=start_list[i]
  end=end_list[i]
  # ins
  ############################################################################################ 20240427
  eigen_input<-read.table(paste0(sample,'_',validpair,'_',chrom,'_ins_cell_sum_25000.xls.gz'),sep = '\t')
  colnames(eigen_input)<-c('pos','eigen','cell')
  
  # add repli-seq annotation
  early_ratio_raw<-read.table('early_ratio.xls.raw.txt',sep = '\t')
  colnames(early_ratio_raw)<-c('cell','early_ratio')
  cellID_mapping<-read.table('hichew_293t_wt_600k_cellID_mapping.xls')
  cell_early_id<-merge(early_ratio_raw,cellID_mapping,by.x='cell',by.y='V1')
  eigen_input
  head(cell_early_id)
  
  # eigen_input_subset<-eigen_input[eigen_input$pos>=floor(start/bin)&eigen_input$pos<=floor(end/bin),]
  eigen_input_subset_tmp<-eigen_input[eigen_input$pos>=floor(start/bin)&eigen_input$pos<=floor(end/bin),]
  eigen_input_subset<-eigen_input_subset_tmp[]
  typeof(eigen_input_subset_tmp)
  
  eigen_input_subset<-eigen_input_subset_tmp[eigen_input_subset_tmp$cell %in% filtered_metrics$sample,]
  
  eigen_merge<-merge(eigen_input_subset,cell_early_id,by='cell')
  eigen_raw <- eigen_merge[order(eigen_merge$early_ratio), ]
  thres=1
  eigen_raw$eigen[eigen_raw$eigen< -thres] <- -thres
  eigen_raw$eigen[eigen_raw$eigen>thres] <- thres
  #ggplot(eigen_raw,aes(x=eigen))+geom_density()
  
  dat_mtx<-dcast(eigen_raw,cell~pos,value.var = "eigen",fun.aggregate = sum)
  dat_mtx_reorder <- dat_mtx[match(unique(eigen_raw$cell), dat_mtx$cell), ]
  
  dat_mtx_<-dat_mtx[,-1]
  rownames(dat_mtx_)<-dat_mtx_reorder$cell
  eigen_raw_subset <- eigen_raw[match(row.names(dat_mtx_), eigen_raw$cell), ]
  annot_row <- data.frame(row.names = eigen_raw_subset$cell, early_ratio = eigen_raw_subset$early_ratio)
  
  ann_colors = list(early_ratio = c("yellow", "firebrick"))
  dat_mtx__ins<-as.matrix(dat_mtx_)
  
  cell_clust<-as.data.frame(cutree(hclust(dist(scale(dat_mtx__ins)),method  = clust_method), k = clust_num))
  cell_clust$cell<-rownames(cell_clust)
  colnames(cell_clust)<-c('cell_cluster','cell')

  cell_clust$group<-'1'
  if (length(rownames(cell_clust[cell_clust$cell_cluster=="1",])) < length(rownames(cell_clust[cell_clust$cell_cluster=="2",]))){
    cell_clust[cell_clust$cell_cluster==2,3]<-'2'  }
  if (length(rownames(cell_clust[cell_clust$cell_cluster=="1",])) > length(rownames(cell_clust[cell_clust$cell_cluster=="2",]))){
    cell_clust[cell_clust$cell_cluster==1,3]<-'2'  }
  cell_clust$ROI<-i
  metrics<-rbind.data.frame(metrics,cell_clust)
  metrics[metrics$group==1,2]
  cell_locus_ins<-melt(dat_mtx__ins)
  head(cell_locus_ins)
  cell_locus_ins_group<-merge(cell_locus_ins,cell_clust,by.x="Var1",by.y="cell")
  cell_locus_ins_group$Var2<-cell_locus_ins_group$Var2*bin
  #ggplot(cell_locus_ins_group,aes(y=value,x=Var2,color=group))+stat_smooth(se=TRUE,span=0.01)+
  #  theme_classic()+ geom_hline(yintercept=c(0),linetype="dotted")+
  #  scale_color_manual(values=c('red3','green4'))+theme(legend.position = "top")

  cell_locus_ins_group_1<-cell_locus_ins_group[cell_locus_ins_group$group==1,]
  summary_stats_1 <- as.data.frame(aggregate(value ~ Var2, data = cell_locus_ins_group_1, FUN = function(x) c(mean = mean(x), sd = sd(x))))
  summary_stats_df_1<-cbind.data.frame(locus=summary_stats_1$Var2,as.data.frame(summary_stats_1$value))
  summary_stats_df_1$group<-'cluster_1'
  
  cell_locus_ins_group_2<-cell_locus_ins_group[cell_locus_ins_group$group==2,]
  summary_stats_2 <- as.data.frame(aggregate(value ~ Var2, data = cell_locus_ins_group_2, FUN = function(x) c(mean = mean(x), sd = sd(x))))
  summary_stats_df_2<-cbind.data.frame(locus=summary_stats_2$Var2,as.data.frame(summary_stats_2$value))
  summary_stats_df_2$group<-'cluster_2'
  
  summary_stats_df<-rbind.data.frame(summary_stats_df_1,summary_stats_df_2)
  p1<-ggplot(summary_stats_df, aes(x=locus, y = mean, ymin = mean - sd, ymax = mean + sd,color=group)) +
    geom_line(alpha=0.8,lwd=1) + geom_errorbar(width = 0.1,alpha=0.1) + labs(x = "Locus", y = "ins") + theme_classic()+
    scale_color_manual(values=c('red3','green4'))+theme(legend.position = "bottom")
  
  dend = as.dendrogram(hclust(dist(scale(dat_mtx__ins),method='euclidean'),method  = clust_method))
  #dend = as.dendrogram(hclust(dist(dat_mtx__ins),method  = clust_method))
  clust_num=2
  dend = color_branches(dend, k = clust_num)

  #col_mean_df<-as.data.frame(colMeans(dat_mtx_))
  #col_mean_df$locus<-as.numeric(rownames(col_mean_df))*bin
  #colnames(col_mean_df)<-c('col_mean','locus')
  #p1<-ggplot(col_mean_df,aes(y=col_mean_df$col_mean,x=locus))+geom_line()+theme_classic()+ geom_hline(yintercept=c(0),linetype="dotted")
  #p1
  
  base_mean = rowMeans(dat_mtx_)
  p<-grid.grabExpr(draw(Heatmap(dat_mtx__ins, name = "ins", cluster_rows = dend, row_split = clust_num,
                                cluster_columns = FALSE,col=colorRampPalette(c("darkgreen",'white','darkred'))(50),
                                show_row_names = FALSE,show_column_names = FALSE,
                                column_title=paste0('cluster_1: ',length(cell_clust[cell_clust$group==1,3]),' cluster_2: ',length(cell_clust[cell_clust$group==2,3]),'\n',
                                  chrom,': ',start,' - ',end,'\nvalidpair downsample: ',
                                                    validpair, '\nclust method: ',clust_method,'\ncell_num: ',
                                                    length(rownames(dat_mtx_)))),heatmap_legend_side = "top"))  
  # eigen
  ############################################################################################ 20240427
  eigen_input<-read.table(paste0(sample,'_',validpair,'_eigen_cell_sum_1000000.xls.gz'),sep = '\t')
  colnames(eigen_input)<-c('pos','eigen','cell')
  eigen_input_subset<-eigen_input[eigen_input$pos>=floor(start/unit)&eigen_input$pos<=floor(end/unit),]
  eigen_input_subset$pos<-eigen_input_subset$pos*unit
  
  cell_locus_ins_group<-merge(eigen_input_subset,cell_clust,by="cell")
  
  cell_locus_ins_group_1<-cell_locus_ins_group[cell_locus_ins_group$group==1,]
  summary_stats_1 <- as.data.frame(aggregate(eigen ~ pos, data = cell_locus_ins_group_1, FUN = function(x) c(mean = mean(x), sd = sd(x))))
  summary_stats_df_1<-cbind.data.frame(locus=summary_stats_1$pos,as.data.frame(summary_stats_1$eigen))
  summary_stats_df_1$group<-'cluster_1'
  length(rownames(cell_locus_ins_group_2))
  cell_locus_ins_group_2<-cell_locus_ins_group[cell_locus_ins_group$group==2,]
  summary_stats_2 <- as.data.frame(aggregate(eigen ~ pos, data = cell_locus_ins_group_2, FUN = function(x) c(mean = mean(x), sd = sd(x))))
  summary_stats_df_2<-cbind.data.frame(locus=summary_stats_2$pos,as.data.frame(summary_stats_2$eigen))
  summary_stats_df_2$group<-'cluster_2'
  
  summary_stats_df<-rbind.data.frame(summary_stats_df_1,summary_stats_df_2)
  p3<-ggplot(summary_stats_df, aes(x=locus, y = mean, ymin = min(mean), ymax = max(mean),color=group)) +
    geom_smooth(alpha=0.8,lwd=1,se=FALSE,span=0.4) + geom_errorbar(width = 0.1,alpha=0.0) + labs(x = "Locus", y = "mean_eigen") + theme_classic()+
    scale_color_manual(values=c('red3','green4'))+theme(legend.position = "bottom")+geom_hline(yintercept=mean(summary_stats_df$mean),linetype="dotted")
  
  
  p2<-cowplot::plot_grid(p,p1,p3, ncol = 1,rel_heights = c(2,1,1), align = "v",axis = "l")
  
  plot_list[[i]]<-p2
}

pdf(paste0(sample,'_',chrom,'_3_ROI_validpair_downsample_sc_ins_141cell_clust2.pdf'),height = 10,width = 18)
cowplot::plot_grid(plotlist = plot_list, nrow = 1)
dev.off()

png(paste0(sample,'_',chrom,'_3_ROI_validpair_downsample_sc_ins_141cell_clust2.png'),height = 1000,width = 1800)
cowplot::plot_grid(plotlist = plot_list, nrow = 1)
dev.off()


group_1_cell_freq<-as.data.frame(table(metrics[metrics$group=="1",2]))
group_1_cell_freq$group<-'cluster_1'
group_2_cell_freq<-as.data.frame(table(metrics[metrics$group=="2",2]))
group_2_cell_freq$group<-'cluster_2'
cell_freq_stack<-rbind.data.frame(group_1_cell_freq,group_2_cell_freq)

ggplot(cell_freq_stack, aes(y = Freq,x = group,fill=group)) + geom_violin(position="dodge", alpha=0.5,trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="cell duplicate Freq in 3 ROI")+ geom_boxplot(width=0.1,outlier.shape = NA)+
  theme_classic()+scale_fill_manual(values=c(mypal[1],mypal[2]))

metrics$count<-1
metrics_cluster_1<-metrics[metrics$group=='1',]
metrics_cluster_1_df<-dcast(metrics_cluster_1,cell~ROI,value.var = "count",fun.aggregate = sum)
colnames(metrics_cluster_1_df)<-c('cell','ROI_1','ROI_2','ROI_3')
metrics_cluster_1_mtx<-as.matrix(metrics_cluster_1_df[,-1])
metrics_cluster_1_mtx_logical<-metrics_cluster_1_mtx ==TRUE

dev.off()

fit1<-euler(metrics_cluster_1_mtx_logical)
pdf(paste0(sample,'_',chrom,'_ROI_3_overlap_ins_509cell_clust1.pdf'),height = 6,width = 6)
plot(fit1,fills = c('lightblue','indianred2','lightgrey'),quantities = TRUE,lty = 1:3,shape = "ellipse",
     main=list(label="cluster_1_overlap", cex = 1))
dev.off()
#,labels = list(font = 20,cex=3)


metrics_cluster_2<-metrics[metrics$group=='2',]
metrics_cluster_2_df<-dcast(metrics_cluster_2,cell~ROI,value.var = "count",fun.aggregate = sum)
colnames(metrics_cluster_2_df)<-c('cell','ROI_1','ROI_2','ROI_3')
metrics_cluster_2_mtx<-as.matrix(metrics_cluster_2_df[,-1])
metrics_cluster_2_mtx_logical<-metrics_cluster_2_mtx ==TRUE

fit2<-euler(metrics_cluster_2_mtx_logical)
pdf(paste0(sample,'_',chrom,'_ROI_3_overlap_ins_509cell_clust2.pdf'),height = 6,width = 6)
plot(fit2,fills = c('lightblue','indianred2','lightgrey'),quantities = TRUE,lty = 1:3,shape = "ellipse",
     main=list(label="cluster_2_overlap", cex = 1))
dev.off()


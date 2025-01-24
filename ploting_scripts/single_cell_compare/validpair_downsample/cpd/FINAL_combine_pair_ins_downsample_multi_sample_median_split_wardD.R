# Load libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library("dendextend") 
library(data.table)
library(IRanges)

filter_overlapping_ranges <- function(df1, df2) {
  ranges1 <- IRanges(start = df1$start, end = df1$end)
  ranges2 <- IRanges(start = df2$start, end = df2$end)
  overlaps <- findOverlaps(ranges1, ranges2)
  overlapping_indices <- queryHits(overlaps)
  df1_filtered <- df1[overlapping_indices, ]
  return(df1_filtered)
}

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
clust_method="ward.D"
chrom="chr7"

unit=1000000

eigen<-fread('snHiChew_HEK293Twt.mcool.res_500000.eig.cis.vecs.tsv',header = T)
eigen_df<-as.data.frame(eigen)
eigen_df$E1[is.na(eigen_df$E1)]<-0
eigen_b<-eigen_df[eigen_df$E1<0&eigen_df$chrom==chrom,]
#metrics<-data.frame()


for(i in 1:length(sample_group)){
  validpair=sample_group[i]
  #start=start_list[i]
  #end=end_list[i]
  # ins
  ############################################################################################ 20240427
  ins_input<-read.table(paste0(sample,'_',validpair,'_ins_cell_sum_25000.xls.gz'),sep = '\t')
  colnames(ins_input)<-c('pos','ins','cell')
  
  # add repli-seq annotation
  #early_ratio_raw<-read.table('early_ratio.xls.raw.txt',sep = '\t')
  #colnames(early_ratio_raw)<-c('cell','early_ratio')
  #cellID_mapping<-read.table('hichew_293t_wt_600k_cellID_mapping.xls')
  #cell_early_id<-merge(early_ratio_raw,cellID_mapping,by.x='cell',by.y='V1')
  
  cell_early_id<-read.table(paste0(sample,'_','label_info.txt'),sep = '\t')
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
  
  
  #ins_input_subset$start<-ins_input_subset$pos*bin
  #ins_input_subset$end<-(ins_input_subset$pos+1)*bin
  #ins <- filter_overlapping_ranges(ins_input_subset, eigen_b)
  ins <- ins_input_subset

  #mean(ins[ins$cell%in%cell_clust_1,2])
  #mean(ins[ins$cell%in%cell_clust_2,2])
  cell_clust$ROI<-i
  #metrics<-rbind.data.frame(metrics,cell_clust)
  #head(metrics)
  cell_locus_ins<-reshape2::melt(dat_mtx__ins)

  cell_locus_ins_group<-merge(cell_locus_ins,cell_clust,by.x="Var1",by.y="cell")
  cell_locus_ins_group$Var2<-cell_locus_ins_group$Var2*bin
  #ggplot(cell_locus_ins_group,aes(y=value,x=Var2,color=group))+stat_smooth(se=TRUE,span=0.01)+
  #  theme_classic()+ geom_hline(yintercept=c(0),linetype="dotted")+
  #  scale_color_manual(values=c('red3','green4'))+theme(legend.position = "top")
  # if (mean(ins[ins$cell%in%cell_clust_1,3]) > mean(ins[ins$cell%in%cell_clust_2,3])){
  if (mean(cell_locus_ins_group[cell_locus_ins_group$group=="1",3]) > mean(cell_locus_ins_group[cell_locus_ins_group$group=="2",3])){
    cell_clust$cell_cluster <- ifelse(cell_clust$group == 1, 2, ifelse(cell_clust$group == 2, 1, cell_clust$cell_cluster)) 
    cell_clust_1<-cell_clust[cell_clust$cell_cluster==1,2]
    cell_clust_2<-cell_clust[cell_clust$cell_cluster==2,2]
    cell_locus_ins_group$group <- ifelse(cell_locus_ins_group$group == 1, 2, ifelse(cell_locus_ins_group$group == 2, 1, cell_locus_ins_group$group)) }
  
  cell_locus_ins_group_1<-cell_locus_ins_group[cell_locus_ins_group$group==1,]
  summary_stats_1 <- as.data.frame(aggregate(value ~ Var2, data = cell_locus_ins_group_1, FUN = function(x) c(mean = mean(x), sd = sd(x))))
  summary_stats_df_1<-cbind.data.frame(locus=summary_stats_1$Var2,as.data.frame(summary_stats_1$value))
  summary_stats_df_1$group<-'cluster_1'
  
  cell_locus_ins_group_2<-cell_locus_ins_group[cell_locus_ins_group$group==2,]
  summary_stats_2 <- as.data.frame(aggregate(value ~ Var2, data = cell_locus_ins_group_2, FUN = function(x) c(mean = mean(x), sd = sd(x))))
  summary_stats_df_2<-cbind.data.frame(locus=summary_stats_2$Var2,as.data.frame(summary_stats_2$value))
  summary_stats_df_2$group<-'cluster_2'
  
  summary_stats_df<-rbind.data.frame(summary_stats_df_1,summary_stats_df_2)
  factor(summary_stats_df$group)
  p1<-ggplot(summary_stats_df, aes(x=locus, y = mean, ymin = mean - sd, ymax = mean + sd,color=group)) +
    geom_line(alpha=0.8,lwd=1) + geom_errorbar(width = 0.1,alpha=0.1) + labs(x = "Locus", y = "ins") + theme_classic()+
    scale_color_manual(values=c('darkred','darkgreen'))+theme(legend.position = "top")+ylim(-0.3,1)

  #p1
  dend = as.dendrogram(hclust(dist(scale(dat_mtx__ins),method='euclidean'),method  = clust_method))
  #dend = as.dendrogram(hclust(dist(dat_mtx__ins,method='euclidean'),method  = clust_method))
  
  clust_num=2
  dend = color_branches(dend, k = clust_num)
    #Heatmap(dat_mtx__ins, name = "ins", cluster_rows = dend,
  #        cluster_columns = FALSE,col=colorRampPalette(c("darkgreen",'white','darkred'))(50),
  #        show_row_names = FALSE,show_column_names = FALSE,
  #        column_title=paste0(chrom,': ',start,' - ',end,'\nvalidpair downsample: ',
  #                            validpair, '\nclust method: ',clust_method,'\ncell_num: ',length(rownames(dat_mtx_))))
  
  
  #col_mean_df<-as.data.frame(colMeans(dat_mtx_))
  #col_mean_df$locus<-as.numeric(rownames(col_mean_df))*bin
  #colnames(col_mean_df)<-c('col_mean','locus')
  #p1<-ggplot(col_mean_df,aes(y=col_mean_df$col_mean,x=locus))+geom_line()+theme_classic()+ geom_hline(yintercept=c(0),linetype="dotted")
  #p1
  
  base_mean = rowMeans(dat_mtx_)
  p<-grid.grabExpr(draw(Heatmap(dat_mtx__ins, name = "ins", cluster_rows = dend, row_split = clust_num,
                                cluster_columns = FALSE,col=colorRampPalette(c("darkgreen",'white','darkred'))(50),
                                show_row_names = FALSE,show_column_names = FALSE,
                                column_title=paste0('cluster_1: ',length(cell_clust[cell_clust$group==1,3]),' cluster_2: ',
                                                    length(cell_clust[cell_clust$group==2,3]),'\n',
                                                    chrom,': ',start,' - ',end,'\nvalidpair downsample: ',
                                                    validpair, '\nclust method: ',clust_method,'\ncell_num: ',
                                                    length(rownames(dat_mtx_)))),heatmap_legend_side = "top"))  

  ############### ins diff
  colnames(cell_locus_ins_group)<-c('cell','pos','ins',"clust", "group","ROI")
  ins_clust_1<-cell_locus_ins_group[cell_locus_ins_group$group==1,]
  ins_clust_1$clust<-'cluster_1'
  ins_clust_2<-cell_locus_ins_group[cell_locus_ins_group$group==2,]
  ins_clust_2$clust<-'cluster_2'  
  
  ins_clust<-rbind.data.frame(ins_clust_1,ins_clust_2)
  ins_clust_mean<-ins_clust %>% group_by(clust,pos) %>% summarize(ins_median=median(ins), .groups = 'keep')
  ins_mean_clust1<-as.data.frame(ins_clust_mean[ins_clust_mean$clust=='cluster_1',])
  ins_mean_clust2<-as.data.frame(ins_clust_mean[ins_clust_mean$clust=='cluster_2',])
  
  sorted_ins_clust_1 <- ins_mean_clust1[order(ins_mean_clust1$ins_median), ]
  sorted_ins_clust_2 <- ins_mean_clust2[order(ins_mean_clust2$ins_median), ]
  
  sorted_ins_clust_1$cumulative_prob <- (1:nrow(sorted_ins_clust_1)) / nrow(sorted_ins_clust_1)
  sorted_ins_clust_2$cumulative_prob <- (1:nrow(sorted_ins_clust_2)) / nrow(sorted_ins_clust_2)
  sorted_ins_clust_b<-rbind.data.frame(sorted_ins_clust_1,sorted_ins_clust_2)
  
  ks_test<-ks.test(sorted_ins_clust_1$ins_median,sorted_ins_clust_2$ins_median)
  ks.test(sorted_ins_clust_1$ins_median,sorted_ins_clust_2$ins_median)
  
  d<-round(as.numeric(ks_test[1]),3)
  pval<- round(-log10(as.numeric(ks_test[2])),3)
  
  factor(sorted_ins_clust_b$clust)
  p3<-ggplot(sorted_ins_clust_b, aes(x = ins_median, y = cumulative_prob,group=clust,color=clust)) +
    geom_line() + labs(title = paste0(sample,"\nCPD Curve\n ins"),subtitle = paste0('d:',d,'  -log10(pval):',pval)) +
    theme_classic()+scale_color_manual(values=c('darkred','darkgreen'))+theme(legend.position = "top")
  p3
  p2<-cowplot::plot_grid(p,p1,p3, ncol = 1,rel_heights = c(2,1,1), align = "v",axis = "l")
  #p3
  plot_list[[i]]<-p2
}

pdf(paste0('split_validpair_downsample_sc_ins_cell_',start,'_',end,'_',clust_method,'.pdf'),height = 18,width = 6)
cowplot::plot_grid(plotlist = plot_list, nrow = 1)
dev.off()

#svg('validpair_downsample_sc_ins_cell.svg',height = 10,width = 18)
#cowplot::plot_grid(plotlist = plot_list, nrow = 1)
#dev.off()

png(paste0('split_validpair_downsample_sc_ins_cell_',start,'_',end,'_',clust_method,'.png'),height = 1800,width = 600)
cowplot::plot_grid(plotlist = plot_list, nrow = 1)
dev.off()

# cell ratio


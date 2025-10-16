library(data.table)
library(ggplot2)
library(purrr)
library(dplyr)
library(ggsci)
library(httpgd)
hgd()
setwd("/research/xieyeming1/proj_2023/hichew_paper_20230710/snHiChew/gb_revision/files/reviewer_2")

cellID_phase_barcode_annot_col<-fread('/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/EDA/time_series/dist_curve_sc_phase/cellID_phase_barcode_annot_col.xls')
hic_pro_metrics_sum<-fread('/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/benchmark/hic_pro_metrics_sum/filtered_snHiChew_HEK293Twt_hic_pro_metrics.xls')
head(hic_pro_metrics_sum)
head(cellID_phase_barcode_annot_col)

hic_pro_metrics_sum_clust<-merge(hic_pro_metrics_sum,cellID_phase_barcode_annot_col,by.x='sample',by.y='V2')
head(hic_pro_metrics_sum_clust)

hic_pro_metrics_sum_clust1<-hic_pro_metrics_sum_clust[hic_pro_metrics_sum_clust$cell_cluster=='1',]
hic_pro_metrics_sum_clust2<-hic_pro_metrics_sum_clust[hic_pro_metrics_sum_clust$cell_cluster=='2',]
write.table(hic_pro_metrics_sum_clust1, file = 'hic_pro_metrics_sum_clust1.xls',sep = "\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(hic_pro_metrics_sum_clust2, file = 'hic_pro_metrics_sum_clust2.xls',sep = "\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


mypal <- pal_npg("nrc", alpha = 0.7)(10) 
mypal <- c(mypal,'darkred','darkgreen','darkblue')

######################### 20240830
batch_f<- list(snHiChew_hek_clust1="hic_pro_metrics_sum_clust1",
               snHiChew_hek_clust2="hic_pro_metrics_sum_clust2")
#                deepSnHichew_1024_0_7Xpure="deepSnHichew_20241024_0_7Xpure",
# zhichao_293tdipC_20241025="empty_filtered_zhichao_293tdipC_20241025_hic_pro_metrics",
#                snHiChew_gm12878_0918="empty_filtered_gm12878_0918_hic_pro_metrics",
batch_data <- map2(batch_f, names(batch_f), function(file, cell_lib) {
  dt <- fread(paste0(file, '.xls'))  # Read the file into a data.table
  dt[, cell_lib := cell_lib]  # Add the label column
  return(dt)})

# Row bind the tables
group_df <- rbindlist(batch_data)
#group_df<-group_df[group_df$Total_pairs_processed<200000&group_df>Total_pairs_processed<3000,]
#group_df[group_df$cell_lib=='deepSnHichewGm12878_20241109',2]
length(which(group_df$cell_lib=='snHiChew_gm12878_0918'))


Total_pairs_processed<-cbind.data.frame(V1=group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p1<-ggplot(data = Total_pairs_processed, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="Total_pairs_processed")+
  theme_classic()+scale_fill_manual(values=mypal)+
  theme(axis.text.x = element_text(angle = 20, hjust = 1))+ theme(legend.position = "none")
p1
#dev.off()

valid_report_ratio<-cbind.data.frame(V1=group_df$valid_pairs/group_df$Reported_pairs,cell_lib=group_df$cell_lib)
p2<-ggplot(data = valid_report_ratio, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="valid / report ratio")+
  theme_classic()+scale_fill_manual(values=mypal)+
  theme(axis.text.x = element_text(angle = 20, hjust = 1))+ theme(legend.position = "none")
p2
#dev.off()

valid_pair_dup<-cbind.data.frame(V1=(group_df$valid_pairs-group_df$valid_interaction_rmdup)/group_df$valid_pairs,cell_lib=group_df$cell_lib)
p3<-ggplot(data = valid_pair_dup, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="valid_pair_dup")+
  theme_classic()+scale_fill_manual(values=mypal)+
  theme(axis.text.x = element_text(angle = 20, hjust = 1))+ theme(legend.position = "none")
p3
#dev.off()

rmdup_validpair_per_read<-cbind.data.frame(V1=group_df$valid_interaction_rmdup/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p4<-ggplot(data = rmdup_validpair_per_read, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5,trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="rmdup valid_pair per read")+ geom_boxplot(width=0.1,outlier.shape = NA)+
  theme_classic()+scale_fill_manual(values=mypal)+
  theme(axis.text.x = element_text(angle = 20, hjust = 1))+ theme(legend.position = "none")
p4
#dev.off()
median(rmdup_validpair_per_read[rmdup_validpair_per_read$cell_lib=='hichew',1])*100000

cis_trans_ratio<-cbind.data.frame(V1=log2((group_df$cis_shortRange+group_df$cis_longRange)/group_df$trans_interaction),cell_lib=group_df$cell_lib)
p5<-ggplot(data = cis_trans_ratio, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="cis_trans_ratio")+
  theme_classic()+scale_fill_manual(values=mypal)+
  theme(axis.text.x = element_text(angle = 20, hjust = 1))+ theme(legend.position = "none")
p5
#dev.off()
# ?cbind.data.frame
report_pair_ratio<-cbind.data.frame(v1=group_df$Reported_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p10<-ggplot(data = report_pair_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="report_pair_ratio")+ 
  theme_classic()+scale_fill_manual(values=mypal)+
  theme(axis.text.x = element_text(angle = 20, hjust = 1))+ theme(legend.position = "none")
p10

# Unmapped_pairs ratio
Unmapped_pair_ratio<-cbind.data.frame(v1=group_df$Unmapped_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p11<-ggplot(data = Unmapped_pair_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="Unmapped pair ratio")+ 
  theme_classic()+scale_fill_manual(values=mypal)+
  theme(axis.text.x = element_text(angle = 20, hjust = 1))+ theme(legend.position = "none")
p11

# low qual pair ratio
Low_qual_pair_ratio<-cbind.data.frame(v1=group_df$Low_qual_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p12<-ggplot(data = Low_qual_pair_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="Low_qual pair ratio")+
  theme_classic()+scale_fill_manual(values=mypal)+
  theme(axis.text.x = element_text(angle = 20, hjust = 1))+ theme(legend.position = "none")
p12

# low qual pair ratio
One_frag_aligned_pair_ratio<-cbind.data.frame(v1=group_df$Pairs_with_singleton/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p13<-ggplot(data = One_frag_aligned_pair_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="One_frag_aligned pair ratio")+
  theme_classic()+scale_fill_manual(values=mypal)+
  theme(axis.text.x = element_text(angle = 20, hjust = 1))+ theme(legend.position = "none")
p13
# ?ggplot
# dev.off()
pdf(paste0('scHiC_compare.pdf'),height=3,width=15)
plotlist<-list(p1,p10,p2,p3,p4,p5)
cowplot::plot_grid(plotlist = plotlist, ncol = 6)
dev.off()

png(paste0('scHiC_compare.png'),height=300,width=1500)
plotlist<-list(p1,p10,p2,p3,p4,p5)
cowplot::plot_grid(plotlist = plotlist, ncol = 6)
dev.off()


pdf(paste0('scHiC_validpair_yield.pdf'),height=10,width=10)
ggplot(data = group_df, aes(y = group_df$valid_interaction_rmdup,x = group_df$Total_pairs_processed,color=cell_lib))+
  geom_point(alpha=0.3)+geom_smooth(alpha=0.3)+xlim(0,2000000)+ylim(0,700000)+
  theme_classic()+scale_color_manual(values=mypal)
dev.off()




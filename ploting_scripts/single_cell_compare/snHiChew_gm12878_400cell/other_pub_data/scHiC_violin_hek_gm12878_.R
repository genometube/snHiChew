library(data.table)
library(ggplot2)
library(purrr)
library(dplyr)
library(ggsci)
mypal <- pal_npg("nrc", alpha = 0.7)(9)
pairs_thres_lower<-10000
pairs_thres_high<-10000
################## knee cut off ##################
batch_f<- list(embryo_dipC_2021="dip_C",
               gage-seq_gm12878='K562_GM12878',
               gm12878_biotin_2020="scihic_gm12878",
               gm12878_LiMCA="LiMCA",
               embryo_nagano_2017="nagano_2017",
               gm12878_SnHichew_1="empty_filtered_gm12878_0918")

batch_data <- map2(batch_f, names(batch_f), function(file, cell_lib) {
  dt <- fread(paste0(file, '_hic_pro_metrics.xls'))  # Read the file into a data.table
  dt[, cell_lib := cell_lib]  # Add the label column
  return(dt)})
group_df <- as.data.frame(rbindlist(batch_data))
group_df<-group_df[group_df$Total_pairs_processed>pairs_thres_lower,]

sample_group<-as.vector(as.data.frame(table(group_df$cell_lib))$Var1)

bw<-0.3
metrics<-data.frame()
plot_list <- list()
for(i in 1:length(sample_group)){
  label=sample_group[i]
  cell_rank_df<-cbind.data.frame(rmd_valid_pair=sort(group_df[group_df$cell_lib==sample_group[i],12],decreasing = TRUE),
                                 cell_rank=c(1:length(group_df[group_df$cell_lib==sample_group[i],12])))
  p<-ggplot(cell_rank_df, aes(x=cell_rank,y=log10(rmd_valid_pair)))+ 
    geom_line()+labs(title=label)+ theme_classic() 
  plot_list[[i]]<-p
}
library(gridExtra)
pdf('rmdValidPair_cellRank_scHiC_compare.pdf',width=8,height=8)
grid.arrange(grobs=plot_list,nrow=2)
dev.off()
svg('rmdValidPair_cellRank_scHiC_compare.svg',width=8,height=8)
grid.arrange(grobs=plot_list,nrow=2)
dev.off()
#cowplot::plot_grid(plotlist = plot_list, nrow = 2)


################## hic pro metrics ##################
color_palette<-mypal
Total_pairs_processed<-cbind.data.frame(V1=group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)

p1<-ggplot(data = Total_pairs_processed, aes(y = V1,x = cell_lib,fill=cell_lib)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="Total_pairs_processed")+theme_classic()+
  theme(axis.text.x = element_text(angle = 15, hjust = 1))+theme(legend.position = "left") +
  scale_fill_manual(values=color_palette)+guides(fill = guide_legend(ncol = 1))+ylim(0,200000)
p1
dev.off()


valid_report_ratio<-cbind.data.frame(V1=group_df$valid_pairs/group_df$Reported_pairs,cell_lib=group_df$cell_lib)
p2<-ggplot(data = valid_report_ratio, aes(y = V1,x = cell_lib,fill=cell_lib)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="valid / report ratio")+
  theme_classic()+theme(axis.text.x = element_text(angle = 15, hjust = 1))+theme(legend.position = "none") +
  scale_fill_manual(values=color_palette)+guides(fill = guide_legend(nrow = 2, byrow = TRUE))
p2
dev.off()

valid_pair_dup<-cbind.data.frame(V1=(group_df$valid_pairs-group_df$valid_interaction_rmdup)/group_df$valid_pairs,cell_lib=group_df$cell_lib)
p3<-ggplot(data = valid_pair_dup, aes(y = V1,x = cell_lib,fill=cell_lib)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="valid_pair_dup")+
  theme_classic()+theme(axis.text.x = element_text(angle = 15, hjust = 1))+theme(legend.position = "none") +
  scale_fill_manual(values=color_palette)+guides(fill = guide_legend(nrow = 2, byrow = TRUE))
p3
dev.off()

rmdup_validpair_per_read<-cbind.data.frame(V1=group_df$valid_interaction_rmdup/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p4<-ggplot(data = rmdup_validpair_per_read, aes(y = V1,x = cell_lib,fill=cell_lib)) + 
  geom_violin(position="dodge", alpha=0.5,trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="rmdup valid_pair per read")+ 
  theme_classic()+theme(axis.text.x = element_text(angle = 15, hjust = 1))+theme(legend.position = "none") +
  scale_fill_manual(values=color_palette)+guides(fill = guide_legend(nrow = 2, byrow = TRUE))
p4
dev.off()
median(rmdup_validpair_per_read[rmdup_validpair_per_read$cell_lib=='hichew',1])*100000

cis_trans_ratio<-cbind.data.frame(V1=log2((group_df$cis_shortRange+group_df$cis_longRange)/group_df$trans_interaction),cell_lib=group_df$cell_lib)
p5<-ggplot(data = cis_trans_ratio, aes(y = V1,x = cell_lib,fill=cell_lib)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="cis_trans_ratio")+
  theme_classic()+theme(axis.text.x = element_text(angle = 15, hjust = 1))+theme(legend.position = "none") +
  scale_fill_manual(values=color_palette)+guides(fill = guide_legend(nrow = 2, byrow = TRUE))
p5
dev.off()

report_pair_ratio<-cbind.data.frame(v1=group_df$Reported_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p10<-ggplot(data = report_pair_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="reported pair ratio")+
  theme_classic()+theme(axis.text.x = element_text(angle = 15, hjust = 1))+theme(legend.position = "none") +
  scale_fill_manual(values=color_palette)+guides(fill = guide_legend(nrow = 2, byrow = TRUE))
p10

Unmapped_pair_ratio<-cbind.data.frame(v1=group_df$Unmapped_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p11<-ggplot(data = Unmapped_pair_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="Unmapped_pair_ratio")+
  theme_classic()+theme(axis.text.x = element_text(angle = 15, hjust = 1))+theme(legend.position = "none") +
  scale_fill_manual(values=color_palette)+guides(fill = guide_legend(nrow = 2, byrow = TRUE))
p11

Low_qual_pair_ratio<-cbind.data.frame(v1=group_df$Low_qual_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p12<-ggplot(data = Low_qual_pair_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="Low_qual_pair_ratio")+
  theme_classic()+theme(axis.text.x = element_text(angle = 15, hjust = 1))+theme(legend.position = "none") +
  scale_fill_manual(values=color_palette)+guides(fill = guide_legend(nrow = 2, byrow = TRUE))

p12

# low qual pair ratio
One_frag_aligned_pair_ratio<-cbind.data.frame(v1=group_df$Pairs_with_singleton/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p13<-ggplot(data = One_frag_aligned_pair_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="One_frag_aligned_pair_ratio")+
  theme_classic()+theme(axis.text.x = element_text(angle = 15, hjust = 1))+theme(legend.position = "none") +
  scale_fill_manual(values=color_palette)+guides(fill = guide_legend(nrow = 2, byrow = TRUE))
p13

dev.off()


pdf(paste0('scHiC_hic_pro_compare.pdf'),height=8,width=18)
#plotlist<-list(p1,p10,p2,p3,p4,p5,p11,p12,p13)
plotlist<-list(p1,p10,p2,p3,p4,p5)
#plotlist<-list(p10,p2,p5)
cowplot::plot_grid(plotlist = plotlist, nrow = 2)
dev.off()


svg(paste0('scHiC_hic_pro_compare.svg'),height=8,width=18)
plotlist<-list(p1,p10,p2,p3,p4,p5)
cowplot::plot_grid(plotlist = plotlist, ncol = 2)
dev.off()







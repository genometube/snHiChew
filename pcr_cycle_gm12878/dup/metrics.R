library(ggplot2)
library(gridExtra)
library(ggsci)
mypal <- pal_npg("nrc", alpha = 0.7)(9)


s1<-read.table('empty_filtered_pcrCycle18_20241012_hic_pro_metrics.xls',header=TRUE)
s1_lab<-'cycle_18'
s1<-s1[s1$Total_pairs_processed>200000&s1$Total_pairs_processed<500000,]

s2<-read.table('empty_filtered_pcrCycle20_20241012_hic_pro_metrics.xls',header=TRUE)
s2_lab<-'cycle_20'
s2<-s2[s2$Total_pairs_processed>200000&s2$Total_pairs_processed<500000,]

# s3<-read.table('cycle_21_hic_pro_metrics.xls',header=TRUE)
# s3_lab<-'cycle_21'
# s3<-s3[s3$Total_pairs_processed>100000,]

s4<-read.table('empty_filtered_gm12878_0918_hic_pro_metrics.xls',header=TRUE)
s4_lab<-'cycle_15'
s4<-s4[s4$Total_pairs_processed>200000&s4$Total_pairs_processed<500000,]

s1$cell_lib<-s1_lab
s2$cell_lib<-s2_lab
# s3$cell_lib<-s3_lab
s4$cell_lib<-s4_lab

bw<-0.3
#group_df<-rbind.data.frame(s1,s2,s3,s4)
#sample_group<-c(s1_lab,s2_lab,s3_lab,s4_lab)

group_df<-rbind.data.frame(s1,s2,s4)
sample_group<-c(s1_lab,s2_lab,s4_lab)

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

pdf('rmdValidPair_cellRank.pdf',width=10,height=10)
grid.arrange(grobs=plot_list,nrow=2)
dev.off()

svg('rmdValidPair_cellRank.svg',width=10,height=10)
grid.arrange(grobs=plot_list,nrow=2)
dev.off()

#cowplot::plot_grid(plotlist = plot_list, nrow = 2)
################## hic pro metrics ##################

color_palette_raw<-c(mypal[1],mypal[2],mypal[4],mypal[9],mypal[5],mypal[3],mypal[6])
color_palette<-color_palette_raw[1:length(sample_group)]
color_palette
group_df$Total_pairs_processed<-as.integer(group_df$Total_pairs_processed)
Total_pairs_processed<-cbind.data.frame(V1=group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
#svg('valid_report_ratio.svg')
p1<-ggplot(data = Total_pairs_processed, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="Total_pairs_processed")+ylim(0,1000000)+
  theme_classic()+scale_fill_manual(values=color_palette)+theme(axis.text.x = element_text(angle = 15, hjust = 1))
p1
dev.off()

#pdf('valid_report_ratio.pdf')
valid_report_ratio<-cbind.data.frame(V1=group_df$valid_pairs/group_df$Reported_pairs,cell_lib=group_df$cell_lib)
#svg('valid_report_ratio.svg')
p2<-ggplot(data = valid_report_ratio, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="valid / report ratio")+
  theme_classic()+scale_fill_manual(values=color_palette)+theme(axis.text.x = element_text(angle = 15, hjust = 1))
p2
dev.off()

#pdf('05.valid_pair_dup_100k.pdf')
valid_pair_dup<-cbind.data.frame(V1=(group_df$valid_pairs-group_df$valid_interaction_rmdup)/group_df$valid_pairs,cell_lib=group_df$cell_lib)
p3<-ggplot(data = valid_pair_dup, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="valid_pair_dup")+ylim(0,1)+
  theme_classic()+scale_fill_manual(values=color_palette)+theme(axis.text.x = element_text(angle = 15, hjust = 1))
p3
dev.off()

#pdf('05.rmdup_validpair_per_read_100k.pdf')
rmdup_validpair_per_read<-cbind.data.frame(V1=group_df$valid_interaction_rmdup/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p4<-ggplot(data = rmdup_validpair_per_read, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5,trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="rmdup valid_pair per read")+ geom_boxplot(width=0.1,outlier.shape = NA)+
  theme_classic()+scale_fill_manual(values=color_palette)+theme(axis.text.x = element_text(angle = 15, hjust = 1))
p4
dev.off()
median(rmdup_validpair_per_read[rmdup_validpair_per_read$cell_lib=='hichew',1])*100000

#pdf('cis_trans_ratio.pdf')
cis_trans_ratio<-cbind.data.frame(V1=log2((group_df$cis_shortRange+group_df$cis_longRange)/group_df$trans_interaction),cell_lib=group_df$cell_lib)
p5<-ggplot(data = cis_trans_ratio, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="cis_trans_ratio")+
  theme_classic()+scale_fill_manual(values=color_palette)+theme(axis.text.x = element_text(angle = 15, hjust = 1))
p5
dev.off()

report_pair_ratio<-cbind.data.frame(v1=group_df$Reported_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p10<-ggplot(data = report_pair_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="reported pair ratio")+theme_classic()+scale_fill_manual(values=color_palette)+
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
p10
dev.off()

Unmapped_pairs_ratio<-cbind.data.frame(v1=group_df$Unmapped_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p6<-ggplot(data = Unmapped_pairs_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="Unmapped pairs ratio")+theme_classic()+scale_fill_manual(values=color_palette)+
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
p6
dev.off()

Low_qual_pairs_ratio<-cbind.data.frame(v1=group_df$Low_qual_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p7<-ggplot(data = Low_qual_pairs_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="Low qual pairs ratio")+theme_classic()+scale_fill_manual(values=color_palette)+
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
p7
dev.off()

Intra_fragment_pairs_ratio<-cbind.data.frame(v1=group_df$Intra_fragment_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p8<-ggplot(data = Intra_fragment_pairs_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="Intra fragment pairs ratio")+theme_classic()+scale_fill_manual(values=color_palette)+
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
p8
dev.off()


pdf(paste0('scHiC_compare_juicer.pdf'),height=8,width=16)
plotlist<-list(p1,p10,p2,p3,p4,p5,p6,p7)
cowplot::plot_grid(plotlist = plotlist, ncol = 3)
dev.off()

svg(paste0('scHiC_compare.svg'),height=8,width=16)
plotlist<-list(p1,p10,p2,p3,p4,p5,p6,p7)
cowplot::plot_grid(plotlist = plotlist, ncol = 3)
dev.off()








################################################################################
################################################################################

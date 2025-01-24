library(ggplot2)
library(ggsci)
mypal <- pal_npg("nrc", alpha = 0.7)(9)
pairs_thres<-5000
s1<-read.table('hic_pro_metrics.xls',header=TRUE)
left<-read.table('cluster_left.txt')
right<-read.table('cluster_right.txt')
color_palette<-mypal

left_hic_pro_metrics<-merge(s1,left,by.x='sample',by.y='V6')
left_hic_pro_metrics$cell_lib<-'left'
right_hic_pro_metrics<-merge(s1,right,by.x='sample',by.y='V6')
right_hic_pro_metrics$cell_lib<-'right'

group_df<-rbind.data.frame(left_hic_pro_metrics,right_hic_pro_metrics)

cis_trans_ratio<-cbind.data.frame(V1=log2((group_df$cis_shortRange+group_df$cis_longRange)/group_df$trans_interaction),cell_lib=group_df$cell_lib)
p5<-ggplot(data = cis_trans_ratio, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="cis_trans_ratio")+
  theme_classic()+scale_fill_manual(values=color_palette)+theme(axis.text.x = element_text(angle = 15, hjust = 1))
p5

log2_trans_cis<-cbind.data.frame(V1=log2(group_df$trans_interaction/(group_df$cis_shortRange+group_df$cis_longRange)),cell_lib=group_df$cell_lib)
p5<-ggplot(data = log2_trans_cis, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="log2_trans_cis")+
  theme_classic()+scale_fill_manual(values=color_palette)+theme(axis.text.x = element_text(angle = 15, hjust = 1))

library(gridExtra)
pdf('trans_cis_snHiChew_testis.pdf',width=4,height=4)
p5
dev.off()
################## knee cut off ##################
s1<-read.table('20240203_testis_L01_hic_pro_metrics.xls',header=TRUE)
s1<-s1[s1$Total_pairs_processed>=pairs_thres,]
s1$read_len<-200
lab1<-'snHiChew_testis_replicate_1'
s1$cell_lib<-lab1

s2<-read.table('testis_replicate_2_hic_pro_metrics.xls',header=TRUE)
s2<-s2[s2$Total_pairs_processed>=pairs_thres,]
s2$read_len<-200
lab2<-'snHiChew_testis_replicate_2'
s2$cell_lib<-lab2

s3<-read.table('tanchen_scHiC10x_20231025_hic_pro_metrics.filter.xls',header=TRUE)
s3<-s3[s3$Total_pairs_processed>=pairs_thres,]
s3$read_len<-200
lab3<-'tanchen_hek293t_dipC_scHiC10x'
s3$cell_lib<-lab3

s4<-read.table('gm12878_0918_hic_pro_metrics.xls',header=TRUE)
s4<-s4[s4$Total_pairs_processed>=pairs_thres,]
s4$read_len<-200
lab4<-'snHiChew_gm12878_0918'
s4$cell_lib<-lab4

s5<-read.table('K562_GM12878_hic_pro_metrics.xls',header=TRUE)
s5<-s5[s5$Total_pairs_processed>=pairs_thres,]
s5$read_len<-200
lab5<-'K562_GM12878'
s5$cell_lib<-lab5

s6<-read.table('scihic_gm12878_hic_pro_metrics.xls',header=TRUE)
s6<-s6[s6$Total_pairs_processed>=pairs_thres,]
s6$read_len<-200
lab6<-'scihic_gm12878'
s6$cell_lib<-lab6

s7<-read.table('LiMCA_hic_pro_metrics.xls',header=TRUE)
s7<-s7[s7$Total_pairs_processed>=pairs_thres,]
s7$read_len<-200
lab7<-'LiMCA_gm12878'
s7$cell_lib<-lab7

group_df<-rbind.data.frame(s1,s2,s3,s4,s5,s6,s7)
sample_group<-c(lab1,lab2,lab3,lab4,lab5,lab6,lab7)

group_df<-rbind.data.frame(s1,s2)
sample_group<-c(lab1,lab2)

bw<-0.3
#group_df<-rbind.data.frame(s1,s3)
#sample_group<-c(lab1,lab3)
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
pdf('rmdValidPair_cellRank_snHiChew_testis.pdf',width=8,height=4)
grid.arrange(grobs=plot_list,nrow=1)
dev.off()
svg('rmdValidPair_cellRank_snHiChew_testis.svg',width=8,height=4)
grid.arrange(grobs=plot_list,nrow=1)
dev.off()
#cowplot::plot_grid(plotlist = plot_list, nrow = 2)


label_testis<-read.table('testis_label_info.txt')
ggplot(label_testis, aes(x=log10(label_testis$V4)))+ 
  geom_density()+labs(title='')+ theme_classic() 
################## hic pro metrics ##################
color_palette<-mypal
Total_pairs_processed<-cbind.data.frame(V1=group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)

p1<-ggplot(data = Total_pairs_processed, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="Total_pairs_processed")+
  theme_classic()+scale_fill_manual(values=color_palette)+theme(axis.text.x = element_text(angle = 15, hjust = 1))
p1
dev.off()


valid_report_ratio<-cbind.data.frame(V1=group_df$valid_pairs/group_df$Reported_pairs,cell_lib=group_df$cell_lib)
p2<-ggplot(data = valid_report_ratio, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="valid / report ratio")+
  theme_classic()+scale_fill_manual(values=color_palette)+theme(axis.text.x = element_text(angle = 15, hjust = 1))
p2
dev.off()

valid_pair_dup<-cbind.data.frame(V1=(group_df$valid_pairs-group_df$valid_interaction_rmdup)/group_df$valid_pairs,cell_lib=group_df$cell_lib)
p3<-ggplot(data = valid_pair_dup, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="valid_pair_dup")+
  theme_classic()+scale_fill_manual(values=color_palette)+theme(axis.text.x = element_text(angle = 15, hjust = 1))
p3
dev.off()

rmdup_validpair_per_read<-cbind.data.frame(V1=group_df$valid_interaction_rmdup/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p4<-ggplot(data = rmdup_validpair_per_read, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5,trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="rmdup valid_pair per read")+ geom_boxplot(width=0.1,outlier.shape = NA)+
  theme_classic()+scale_fill_manual(values=color_palette)+theme(axis.text.x = element_text(angle = 15, hjust = 1))
p4
dev.off()
median(rmdup_validpair_per_read[rmdup_validpair_per_read$cell_lib=='hichew',1])*100000

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

Unmapped_pair_ratio<-cbind.data.frame(v1=group_df$Unmapped_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p11<-ggplot(data = Unmapped_pair_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="Unmapped_pair_ratio")+theme_classic()+scale_fill_manual(values=color_palette)+
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
p11

Low_qual_pair_ratio<-cbind.data.frame(v1=group_df$Low_qual_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p12<-ggplot(data = Low_qual_pair_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="Low_qual_pair_ratio")+theme_classic()+scale_fill_manual(values=color_palette)+
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

p12

# low qual pair ratio
One_frag_aligned_pair_ratio<-cbind.data.frame(v1=group_df$Pairs_with_singleton/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p13<-ggplot(data = One_frag_aligned_pair_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="One_frag_aligned_pair_ratio")+theme_classic()+scale_fill_manual(values=color_palette)+
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
p13

dev.off()


pdf(paste0('scHiC_hic_pro_compare.pdf'),height=10,width=10)
plotlist<-list(p1,p10,p2,p3,p4,p5,p11,p12,p13)
plotlist<-list(p1,p10,p2,p3,p4,p5)
plotlist<-list(p10,p2,p5)
cowplot::plot_grid(plotlist = plotlist, ncol = 2)
dev.off()


svg(paste0('scHiC_hic_pro_compare.svg'),height=10,width=10)
plotlist<-list(p10,p2,p5)
plotlist<-list(p1,p10,p2,p3,p4,p5)
cowplot::plot_grid(plotlist = plotlist, ncol = 2)
dev.off()






pdf(paste0('scHiC_read_validpair_test.pdf'),height=10,width=20)
ggplot(data = group_df, aes(y = group_df$valid_interaction_rmdup,x = group_df$Total_pairs_processed,color=group_df$cell_lib)) + 
  geom_point(alpha=0.3)+geom_smooth(method='loess',span=1)+theme_classic()+xlim(0,20000000)+ylim(0,1000000)+scale_color_manual(values=color_palette)
dev.off()


pdf(paste0('scHiC_validpair_dup_test.pdf'),height=10,width=20)
ggplot(data = group_df, aes(y = (group_df$valid_pairs - group_df$valid_interaction_rmdup)/group_df$valid_pairs,x = group_df$Total_pairs_processed,color=group_df$cell_lib)) + 
  geom_point(alpha=0.3)+geom_smooth(method='loess',span=1)+theme_classic()+xlim(0,5000000)
dev.off()






################################################################################
################################################################################
################################################################################
# valid_pairs per gb
contact_gb_df<-cbind.data.frame(contact=group_df$valid_pairs,gb=group_df$Total_pairs_processed*group_df$read_len/1000000000,cell_lib=group_df$cell_lib)
contact_per_gb<-cbind.data.frame(V1=contact_gb_df$contact/contact_gb_df$gb/100000,cell_lib=group_df$cell_lib)
ggplot(data = contact_per_gb, aes(y = V1,x = cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="valid_pair per Gb (1 * 100000)")


rmdup_contact_gb_df<-cbind.data.frame(contact=group_df$valid_interaction_rmdup,gb=group_df$Total_pairs_processed*group_df$read_len/1000000000,cell_lib=group_df$cell_lib)
rmdup_contact_per_gb<-cbind.data.frame(V1=rmdup_contact_gb_df$contact/rmdup_contact_gb_df$gb/100000,cell_lib=group_df$cell_lib)
p1<-ggplot(data = rmdup_contact_per_gb, aes(y = V1,x = cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="rmdup valid_pair per Gb (1 * 100000)")
p1
# valid_pair_per_read
validpair_per_read<-cbind.data.frame(V1=group_df$valid_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p15<-ggplot(data = validpair_per_read, aes(y = V1,x = cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
   geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="valid_pair per read")
p15
# valid_pair frag size
validpair_per_read<-cbind.data.frame(V1=group_df$mean_frag_size,cell_lib=group_df$cell_lib)
p2<-ggplot(data = validpair_per_read, aes(y = V1,x = cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="valid_pair frag size")
p2

rmdup_validpair_per_read<-cbind.data.frame(V1=group_df$valid_interaction_rmdup/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p4<-ggplot(data = rmdup_validpair_per_read, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5,trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="rmdup valid_pair per read")+ geom_boxplot(width=0.1,outlier.shape = NA)+
  theme_classic()+scale_fill_manual(values=c(mypal[1],mypal[2],mypal[4]))
p4

# log2(cis-to-trans ratio)
cis_trans_ratio<-cbind.data.frame(V1=log2((group_df$cis_shortRange+group_df$cis_longRange)/group_df$trans_interaction),cell_lib=group_df$cell_lib)
p3<-ggplot(data = cis_trans_ratio, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5,trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="log2 cis_trans_ratio")+ geom_boxplot(width=0.1,outlier.shape = NA)+
  theme_classic()+scale_fill_manual(values=c(mypal[1],mypal[2],mypal[4]))
p3


# rmdup valid_pair
rmdup_validpair_per_cell<-cbind.data.frame(v1=log10(group_df$valid_interaction_rmdup),cell_lib=group_df$cell_lib)
p6<-ggplot(data = rmdup_validpair_per_cell, aes(y=v1,x = cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="rmd_valid_pair log10")
p6

# coverage per cell
cov<-cbind.data.frame(v1=group_df$validpair_genome_frag_cov,cell_lib=group_df$cell_lib)
p5<-ggplot(data = cov, aes(y=v1,x = cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="cov")
p5
cov_per_read<-cbind.data.frame(v1=group_df$validpair_genome_frag_cov/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p6<-ggplot(data = cov_per_read, aes(y=v1,x = cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="cov_per_read")
p6

# validpair dup
validpair_dup<-cbind.data.frame(v1=((group_df$valid_pairs-group_df$valid_interaction_rmdup)/group_df$valid_pairs),cell_lib=group_df$cell_lib)
p7<-ggplot(data = validpair_dup, aes(y=v1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5,trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="valid_pair dup")+ geom_boxplot(width=0.1,outlier.shape = NA)+
  theme_classic()+scale_fill_manual(values=c(mypal[1],mypal[2],mypal[4]))
p7
# totalread
total<-cbind.data.frame(v1=log10(group_df$Total_pairs_processed),cell_lib=group_df$cell_lib)
p8<-ggplot(data = total, aes(y=v1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.02) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="total read log10")+ geom_boxplot(width=0.1,outlier.shape = NA)+
  theme_classic()+scale_fill_manual(values=c(mypal[1],mypal[2],mypal[4]))
p8

# DE gatc ratio
DE_gatc_ratio<-cbind.data.frame(v1=group_df$Dangling_end_pairs_gatc_ratio,cell_lib=group_df$cell_lib)
p9<-ggplot(data = DE_gatc_ratio, aes(y=v1,x = cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="dangling end pair gatc ratio")
p9

# reported pair ratio
report_pair_ratio<-cbind.data.frame(v1=group_df$Reported_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p10<-ggplot(data = report_pair_ratio, aes(y=v1,x = cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="reported pair ratio")
p10

# Unmapped_pairs ratio
Unmapped_pair_ratio<-cbind.data.frame(v1=group_df$Unmapped_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p11<-ggplot(data = Unmapped_pair_ratio, aes(y=v1,x = cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="Unmapped pair ratio")
p11

# low qual pair ratio
Low_qual_pair_ratio<-cbind.data.frame(v1=group_df$Low_qual_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p12<-ggplot(data = Low_qual_pair_ratio, aes(y=v1,x = cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="Low_qual pair ratio")
p12

# low qual pair ratio
One_frag_aligned_pair_ratio<-cbind.data.frame(v1=group_df$Pairs_with_singleton/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p13<-ggplot(data = Low_qual_pair_ratio, aes(y=v1,x = cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="One_frag_aligned pair ratio")
p13

# valid pair / report pair ratio
valid_report_ratio<-cbind.data.frame(v1=group_df$valid_pairs/group_df$Reported_pairs,cell_lib=group_df$cell_lib)
p14<-ggplot(data = valid_report_ratio, aes(y=v1,x = cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="valid_report ratio")
p14

dev.off()
pdf('scHiC_metrics_',label_hichew_test,'.pdf',height=18,width=18)
plotlist<-list(p8,p14,p10,p15,p11,p12,p13,p6,p4,p1,p2,p7,p3,p5,p9)
cowplot::plot_grid(plotlist = plotlist, nrow = 5)
dev.off()

pdf('scHiC_compare_',label_hichew_test,'.pdf',height=10,width=15)
plotlist<-list(p10,p11,p12,p13,p14)
cowplot::plot_grid(plotlist = plotlist, nrow = 2)
dev.off()

pdf('scHiC_metrics_hichew_wengzhe_20240112.pdf',height=10,width=10)
plotlist<-list(p8,p4,p7,p3)
cowplot::plot_grid(plotlist = plotlist, nrow = 2)
dev.off()



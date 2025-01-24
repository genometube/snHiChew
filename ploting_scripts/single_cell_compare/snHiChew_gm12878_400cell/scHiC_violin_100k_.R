library(ggplot2)

library(ggsci)
mypal <- pal_npg("nrc", alpha = 0.7)(9)

################## knee cut off ##################
nagano_2017_raw<-read.table('nagano_2017_hic_pro_metrics_100k.xls',header=TRUE)
nagano_2017<-nagano_2017_raw[nagano_2017_raw$Total_pairs_processed>=100000,]

dip_C_raw<-read.table('dip_C_hic_pro_metrics_100k.xls',header=TRUE)
dip_C<-dip_C_raw[dip_C_raw$Total_pairs_processed>=100000,]

hichew_raw<-read.table('filtered_snHiChew_HEK293Twt_hic_pro_metrics.xls',header = TRUE)
hichew<-hichew_raw[hichew_raw$Total_pairs_processed>=100000,]
#scCARE<-read.table('../hichew_paper_2023/scCAREseq/single_cell_metrics/downsample/hic_pro_metrics.xls',header=TRUE)
#hc_1013_chr7<-read.table('../../../zhichao_porecHiC_20231011/cell_split/sc_metrics_chr7/hic_pro_metrics.xls',header=TRUE)
library(data.table)
#hichew_test<-read.table('single_cut_0729_hic_pro_metrics.xls',header = TRUE)

hichew_test<-read.table('empty_filtered_zhichao_400cellgm12878b1_20240729_hic_pro_metrics.xls',header = TRUE)
# hichew_test_raw<-as.data.frame(hichew_test_dt[-which(is.na(hichew_test_dt$Total_pairs_processed)),])
# hichew_test<-as.data.frame(hichew_test_raw)
label_hichew_test<-'zhichao_400cellgm12878b1'
mean(hichew_test$Total_pairs_processed)

hichew_test2<-read.table('lane2_empty_filtered_zhichao_200cell293Tb1_20240729_hic_pro_metrics.xls',header = TRUE)
#hichew_test_raw2<-as.data.frame(hichew_test_dt2[-which(is.na(hichew_test_dt2$Total_pairs_processed)),])
#hichew_test2<-as.data.frame(hichew_test_raw2)
label_hichew_test2<-'snHiChew_200cell_lane2'

deepSnHichew_0428<-read.table('sorted_deepSnHichew_0428_cut.xls',header = TRUE)
noIP_10x<-read.table('tanchen_scHiC10x_20231025_hic_pro_metrics.filter.xls',header = TRUE)


nagano_2017$read_len<-300
hichew$read_len<-200
dip_C$read_len<-300
hichew_test$read_len<-200
hichew_test2$read_len<-200

noIP_10x$read_len<-200
deepSnHichew_0428$read_len<-200

nagano_2017$cell_lib<-'nagano_2017'
hichew$cell_lib<-'hichew'
dip_C$cell_lib<-'dip_C'
hichew_test$cell_lib<-label_hichew_test
hichew_test2$cell_lib<-label_hichew_test2

noIP_10x$cell_lib<-'noIP_10x'
deepSnHichew_0428$cell_lib<-'deepSnHichew_0428'

bw<-0.3
group_df<-rbind.data.frame(nagano_2017,dip_C,hichew,hichew_test,hichew_test2,noIP_10x,deepSnHichew_0428)
sample_group<-c('nagano_2017','dip_C','hichew',label_hichew_test,label_hichew_test2,'noIP_10x','deepSnHichew_0428')

sort(nagano_2017$valid_interaction_rmdup,decreasing = TRUE)
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
pdf('rmdValidPair_cellRank.pdf',width=10,height=10)
grid.arrange(grobs=plot_list,nrow=2)
dev.off()
#cowplot::plot_grid(plotlist = plot_list, nrow = 2)





################## hic pro metrics ##################

#group_df<-rbind.data.frame(nagano_2017,hichew,dip_C,sorted_cut_hichew_test,sorted_cut_hichew_test2,noIP_10x,deepSnHichew_0428)
sample_group<-c('nagano_2017','hichew','dip_C',label_hichew_test,label_hichew_test2,'noIP_10x','deepSnHichew_0428')


#mean(group_df$valid_interaction_rmdup)
# valid_report_ratio
#median(valid_report_ratio[valid_report_ratio$cell_lib=="noIP_10x",1])
#sample_group<-c('nagano_2017','dip_C','hichew',label_hichew_test,'noIP_10x')
#pdf('Total_pairs_processed.pdf')
color_palette<-c(mypal[1],mypal[2],mypal[4],mypal[9],mypal[5],mypal[3],mypal[6])
Total_pairs_processed<-cbind.data.frame(V1=group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
#svg('valid_report_ratio.svg')
p1<-ggplot(data = Total_pairs_processed, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="Total_pairs_processed")+
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
  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="valid_pair_dup")+
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

Unmapped_pair_ratio<-cbind.data.frame(v1=group_df$Unmapped_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p11<-ggplot(data = Unmapped_pair_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="Unmapped_pair_ratio")+theme_classic()+scale_fill_manual(values=color_palette)+
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
p11

# low qual pair ratio
Low_qual_pair_ratio<-cbind.data.frame(v1=group_df$Low_qual_pairs/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p12<-ggplot(data = Low_qual_pair_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="Low_qual_pair_ratio")+theme_classic()+scale_fill_manual(values=color_palette)+
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

p12

# low qual pair ratio
One_frag_aligned_pair_ratio<-cbind.data.frame(v1=group_df$Pairs_with_singleton/group_df$Total_pairs_processed,cell_lib=group_df$cell_lib)
p13<-ggplot(data = Low_qual_pair_ratio, aes(y=v1,x = cell_lib,fill=cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE,bw=0.005) +
  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="One_frag_aligned_pair_ratio")+theme_classic()+scale_fill_manual(values=color_palette)+
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
p11
p13


dev.off()
pdf(paste0('scHiC_compare_',label_hichew_test,'.pdf'),height=8,width=16)
plotlist<-list(p1,p10,p2,p3,p4,p5,p11,p12,p13)
cowplot::plot_grid(plotlist = plotlist, ncol = 3)
dev.off()



ggplot(data = group_df, aes(y = group_df$valid_interaction_rmdup,x = group_df$Total_pairs_processed,color=group_df$cell_lib)) + 
  geom_point(alpha=0.3)+geom_smooth(method='loess',span=1)+theme_classic()+xlim(0,5000000)


ggplot(data = group_df, aes(y = (group_df$valid_pairs - group_df$valid_interaction_rmdup)/group_df$valid_pairs,x = group_df$Total_pairs_processed,color=group_df$cell_lib)) + 
  geom_point(alpha=0.3)+geom_smooth(method='loess',span=1)+theme_classic()+xlim(0,5000000)

?geom_smooth




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

# valid_pair
#validpair_per_cell<-cbind.data.frame(v1=log10(group_df$valid_pairs),cell_lib=group_df$cell_lib)
#p15<-ggplot(data = validpair_per_cell, aes(y=v1,x = cell_lib)) + geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",trim=FALSE) +
#  geom_boxplot(width=0.1,outlier.shape = NA) +labs(title="valid_pair log10")
#p15
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



library(plotly)
plot_ly(
  hichew_H_M, x = ~hg19, y = ~mm10, z = ~hg19_mm10, size=0.05,
  colors = c('#BF382A', '#0C4B8E','yellow')
) %>%
  add_markers() %>%
  layout(
    scene = list(xaxis = list(title = 'hg19'),
                 yaxis = list(title = 'mm10'),
                 zaxis = list(title = 'hg19_mm10'))
  )


################## knee cut off ##################
library(kneedle)
library(ggsci)
raw_metrics<-hichew_test
cell_rank_df<-cbind.data.frame(rmd_valid_pair=sort(raw_metrics$valid_interaction_rmdup,decreasing = TRUE),
                               cell_rank=c(1:length(raw_metrics$valid_interaction_rmdup)))
rank_df_y_axis<-sort(raw_metrics$valid_interaction_rmdup,decreasing = TRUE)
rank_df_y_axis[10:length(rank_df_y_axis)]
knee1<-kneedle(c(1:length(raw_metrics$valid_interaction_rmdup)),rank_df_y_axis, 
               concave = FALSE, decreasing = TRUE,sensitivity = 20)
knee1
start2<-knee1[1]+100
knee2<-kneedle(c(start2:length(raw_metrics$valid_interaction_rmdup)),rank_df_y_axis[start2:length(rank_df_y_axis)] ,
               concave = FALSE, decreasing = TRUE,sensitivity = 15)
knee2
pdf(paste0(label_hichew_test,'_cell_rank.pdf'),width=4,height=4)
ggplot(cell_rank_df, aes(x=cell_rank,y=log10(rmd_valid_pair)))+ 
  geom_line()+labs(title=knee2[1])+ theme_classic()+geom_vline(xintercept= knee2[1], linetype='dotted')
dev.off()

sorted_metrics <- raw_metrics[order(raw_metrics$valid_interaction_rmdup,decreasing=TRUE), ]
rownames(sorted_metrics)<-c(1:length(rownames(sorted_metrics)))
sorted_cut_hichew_test <- sorted_metrics[c(1:knee2[1]),]
write.table(sorted_cut_hichew_test,file = paste0('sorted_',label_hichew_test,'_cut.xls'),quote = FALSE,sep="\t",row.names = FALSE)


raw_metrics<-hichew_test2
cell_rank_df<-cbind.data.frame(rmd_valid_pair=sort(raw_metrics$valid_interaction_rmdup,decreasing = TRUE),
                               cell_rank=c(1:length(raw_metrics$valid_interaction_rmdup)))
rank_df_y_axis<-sort(raw_metrics$valid_interaction_rmdup,decreasing = TRUE)
rank_df_y_axis[10:length(rank_df_y_axis)]
knee1<-kneedle(c(1:length(raw_metrics$valid_interaction_rmdup)),rank_df_y_axis, 
               concave = FALSE, decreasing = TRUE,sensitivity = 20)
knee1
start2<-knee1[1]+100
knee2<-kneedle(c(start2:length(raw_metrics$valid_interaction_rmdup)),rank_df_y_axis[start2:length(rank_df_y_axis)] ,
               concave = FALSE, decreasing = TRUE,sensitivity = 15)
knee2
pdf(paste0(label_hichew_test2,'_cell_rank.pdf'),width=4,height=4)
ggplot(cell_rank_df, aes(x=cell_rank,y=log10(rmd_valid_pair)))+ 
  geom_line()+labs(title=knee2[1])+ theme_classic()+geom_vline(xintercept= knee2[1], linetype='dotted')
dev.off()

sorted_metrics <- raw_metrics[order(raw_metrics$valid_interaction_rmdup,decreasing=TRUE), ]
rownames(sorted_metrics)<-c(1:length(rownames(sorted_metrics)))
sorted_cut_hichew_test2 <- sorted_metrics[c(1:knee2[1]),]
write.table(sorted_cut_hichew_test2,file = paste0('sorted_',label_hichew_test2,'_cut.xls'),quote = FALSE,sep="\t",row.names = FALSE)


################## knee cut off ##################
library(kneedle)
library(ggsci)
raw_metrics<-hichew_test
cell_rank_df<-cbind.data.frame(rmd_valid_pair=sort(raw_metrics$valid_interaction_rmdup,decreasing = TRUE),
                               cell_rank=c(1:length(raw_metrics$valid_interaction_rmdup)))
rank_df_y_axis<-sort(raw_metrics$valid_interaction_rmdup,decreasing = TRUE)
rank_df_y_axis[10:length(rank_df_y_axis)]
knee1<-kneedle(c(1:length(raw_metrics$valid_interaction_rmdup)),rank_df_y_axis, 
               concave = FALSE, decreasing = TRUE,sensitivity = 20)
knee1
start2<-knee1[1]+100
knee2<-kneedle(c(start2:length(raw_metrics$valid_interaction_rmdup)),rank_df_y_axis[start2:length(rank_df_y_axis)] ,
               concave = FALSE, decreasing = TRUE,sensitivity = 15)
knee2
pdf(paste0(label_hichew_test,'_cell_rank.pdf'),width=4,height=4)
ggplot(cell_rank_df, aes(x=cell_rank,y=log10(rmd_valid_pair)))+ 
  geom_line()+labs(title=knee2[1])+ theme_classic()+geom_vline(xintercept= knee2[1], linetype='dotted')
dev.off()

sorted_metrics <- raw_metrics[order(raw_metrics$valid_interaction_rmdup,decreasing=TRUE), ]
rownames(sorted_metrics)<-c(1:length(rownames(sorted_metrics)))
sorted_cut_hichew_test <- sorted_metrics[c(1:knee2[1]),]
write.table(sorted_cut_hichew_test,file = paste0('sorted_',label_hichew_test,'_cut.xls'),quote = FALSE,sep="\t",row.names = FALSE)

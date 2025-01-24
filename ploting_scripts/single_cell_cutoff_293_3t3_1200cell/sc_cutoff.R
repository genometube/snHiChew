library(patchwork)
library(gridExtra)
library(gtable)
library(ggplot2)
library(tibble)
library(kneedle)
library(ggsci)
mypal <- pal_npg("nrc", alpha = 0.7)(9)

################## knee cut off human mouse ##################
hc_293_3t3<-read.table('hichew_293_3t3_hic_pro_metrics.xls',header=TRUE)

cell_rank_df<-cbind.data.frame(rmd_valid_pair=sort(hc_293_3t3$valid_interaction_rmdup,decreasing = TRUE),
                               cell_rank=c(1:length(hc_293_3t3$valid_interaction_rmdup)))
rank_df_y_axis<-sort(hc_293_3t3$valid_interaction_rmdup,decreasing = TRUE)
rank_df_y_axis[10:length(rank_df_y_axis)]
knee1<-kneedle(c(1:length(hc_293_3t3$valid_interaction_rmdup)),log10(rank_df_y_axis[1:length(hc_293_3t3$valid_interaction_rmdup)]), 
               concave = FALSE, decreasing = TRUE,sensitivity = 6)
knee1
start2<-knee1[1]+200
knee2<-kneedle(c(start2:length(hc_293_3t3$valid_interaction_rmdup)),log10(rank_df_y_axis[start2:length(rank_df_y_axis)]) ,
               concave = FALSE, decreasing = TRUE,sensitivity = 20)
knee2
pdf('293_3t3_cell_rank.pdf',width=4,height=4)
ggplot(cell_rank_df, aes(x=cell_rank,y=log10(rmd_valid_pair)))+xlim(0,1735)+
  geom_line()+labs(title=paste0('N = 1735, \ncut off rank 1: ',knee1[1],'\ncut off rank 2: ',knee2[1],'\ncut off 2 validpair:',10^knee2[2]))+ 
  theme_classic()+geom_vline(xintercept= c(knee1[1],knee2[1]), linetype='dotted')

p1<-ggplot(cell_rank_df, aes(x=cell_rank,y=(rmd_valid_pair)))+xlim(0,1735)+ xlim(0,1900)+
  geom_line()+labs(title=paste0('N = 1735, \ncut off rank 1: ',knee1[1],'\ncut off rank 2: ',knee2[1],'\ncut off 2 validpair:',10^knee2[2]))+ 
  theme_classic()+geom_vline(xintercept= c(knee1[1],knee2[1]), linetype='dotted')
dev.off()
svg('293_3t3_cell_rank.svg',width=4,height=4)
p1
dev.off()
################## doublet_density ##################
sorted_hc_293_3t3 <- hc_293_3t3[order(hc_293_3t3$valid_interaction_rmdup,decreasing=TRUE), ]
rownames(sorted_hc_293_3t3)<-c(1:length(rownames(sorted_hc_293_3t3)))
#clean_hc_293_3t3 <- sorted_hc_293_3t3
#subset_hc_293_3t3 <- sorted_hc_293_3t3
#sorted_hc_293_3t3 <- clean_hc_293_3t3

sorted_hc_293_3t3$group<-'empty'
knee1

knee2
knee2_thres<-10^(knee2[2])-4000

knee2_thres
sorted_hc_293_3t3[sorted_hc_293_3t3$valid_interaction_rmdup>knee2_thres&sorted_hc_293_3t3$hg19/sorted_hc_293_3t3$mm10 > 20,length(colnames(sorted_hc_293_3t3))] <-'hg19'
sorted_hc_293_3t3[sorted_hc_293_3t3$valid_interaction_rmdup>knee2_thres&sorted_hc_293_3t3$hg19/sorted_hc_293_3t3$mm10 < 0.5,length(colnames(sorted_hc_293_3t3))] <-'mm10'
sorted_hc_293_3t3[sorted_hc_293_3t3$valid_interaction_rmdup>knee2_thres+4000&sorted_hc_293_3t3$hg19/sorted_hc_293_3t3$mm10 < 20 & sorted_hc_293_3t3$hg19/sorted_hc_293_3t3$mm10 > 0.5,length(colnames(sorted_hc_293_3t3))] <-'doublet'

hg19_num<-length(which(sorted_hc_293_3t3$group=="hg19"))
mm10_num<-length(which(sorted_hc_293_3t3$group=="mm10"))
doublet_num<-length(which(sorted_hc_293_3t3$group=="doublet"))
empty_num<-length(which(sorted_hc_293_3t3$group=="empty"))
double_rate<-doublet_num/(doublet_num+mm10_num+hg19_num)
doublet_num+mm10_num+hg19_num
double_rate
mm10_num+hg19_num

doublet_den<-cbind.data.frame(cell_rank=as.integer(rownames(sorted_hc_293_3t3)),cell_type=sorted_hc_293_3t3$group,doublet_count=0)
doublet_den[doublet_den$cell_type=="doublet",3]<-1
pdf('doublet_density.pdf',width=4,height=4)
p2<-ggplot(doublet_den, aes(x=cell_rank,y=doublet_count))+xlim(0,1735)+geom_smooth(alpha=0.7,span=0.5,method='loess',color=mypal[1])+
  labs(title=paste0('doublet density distribution   N = 1735, \ncut off rank 1: ',knee1[1],'\ncut off rank 2: ',knee2[1]))+ 
  theme_classic()+geom_vline(xintercept= c(knee1[1],knee2[1]), linetype='dotted')
dev.off()
svg('doublet_density.svg',width=4,height=4)
p2
dev.off()

################## 2 knee cut off human mouse ##################
sorted_hc_293_3t3 <- hc_293_3t3[order(hc_293_3t3$valid_interaction_rmdup,decreasing=TRUE), ]
rownames(sorted_hc_293_3t3)<-c(1:length(rownames(sorted_hc_293_3t3)))
clean_hc_293_3t3 <- sorted_hc_293_3t3[c(knee1[1]:length(rownames(sorted_hc_293_3t3))),]
subset_hc_293_3t3 <- sorted_hc_293_3t3
sorted_hc_293_3t3 <- clean_hc_293_3t3

sorted_hc_293_3t3$group<-'empty'
knee1

knee2
knee2_thres<-10^(knee2[2])-4000

knee2_thres
sorted_hc_293_3t3[sorted_hc_293_3t3$valid_interaction_rmdup>knee2_thres&sorted_hc_293_3t3$hg19/sorted_hc_293_3t3$mm10 > 20,length(colnames(sorted_hc_293_3t3))] <-'hg19'
sorted_hc_293_3t3[sorted_hc_293_3t3$valid_interaction_rmdup>knee2_thres&sorted_hc_293_3t3$hg19/sorted_hc_293_3t3$mm10 < 0.5,length(colnames(sorted_hc_293_3t3))] <-'mm10'
sorted_hc_293_3t3[sorted_hc_293_3t3$valid_interaction_rmdup>knee2_thres+4000&sorted_hc_293_3t3$hg19/sorted_hc_293_3t3$mm10 < 20 & sorted_hc_293_3t3$hg19/sorted_hc_293_3t3$mm10 > 0.5,length(colnames(sorted_hc_293_3t3))] <-'doublet'

hg19_num<-length(which(sorted_hc_293_3t3$group=="hg19"))
mm10_num<-length(which(sorted_hc_293_3t3$group=="mm10"))
doublet_num<-length(which(sorted_hc_293_3t3$group=="doublet"))
empty_num<-length(which(sorted_hc_293_3t3$group=="empty"))
double_rate<-doublet_num/(doublet_num+mm10_num+hg19_num)
doublet_num+mm10_num+hg19_num+empty_num
double_rate
mm10_num+hg19_num


ggplot(sorted_hc_293_3t3, aes(x=hg19,y=mm10,color=group)) + geom_point(alpha=0.5,stroke=NA)+
  theme_classic() +scale_color_manual(values=c(mypal[1],'grey',mypal[3],mypal[4]))+
  labs(title=paste0('loaded cell: 1200, demultiplexed barcode: top 1735
hg19: ',hg19_num,'  mm10:',mm10_num,'  empty:',empty_num,'  doublet:',doublet_num,'\ndoublet rate ',round(double_rate,3), ' doublet cut off: ', knee1[1] ))


p3<-ggplot(sorted_hc_293_3t3, aes(x=hg19,y=mm10,color=group)) + geom_point(alpha=0.5,stroke=NA,size=1.2)+
  theme_classic() +scale_color_manual(values=c(mypal[1],'grey',mypal[3],mypal[4]))+
  labs(title=paste0('loaded cell: 1200, demultiplexed barcode: top 1735
hg19: ',hg19_num,'  mm10:',mm10_num,'  empty:',empty_num,'  doublet:',doublet_num,'\ndoublet rate ',round(double_rate,3), ' doublet cut off: ', knee1[1] ))

pdf('cell_1200_293_3t3_doublet.pdf',width=5,height=5)
p3
dev.off()
svg('cell_1200_293_3t3_doublet.svg',width=5,height=5)
p3
dev.off()

write.table(sorted_hc_293_3t3,file = 'sorted_hc_293_3t3.xls',quote = FALSE,sep="\t",row.names = FALSE)

################## cis short, cis long, trans, interspecies ##################
cis_trans<-cbind.data.frame(cis=clean_hc_293_3t3$cis_longRange + clean_hc_293_3t3$cis_shortRange,
                 trans=clean_hc_293_3t3$trans_interaction-clean_hc_293_3t3$hg19_mm10,
                 interspecies=clean_hc_293_3t3$hg19_mm10,
                 rmd_validpair=clean_hc_293_3t3$valid_interaction_rmdup)
cis_trans_pct<-cis_trans/cis_trans$rmd_validpair
cis_trans_pct<-cis_trans_pct[,-4]
cis_pct<-round(mean(cis_trans_pct$cis),3)
trans_pct<-round(mean(cis_trans_pct$trans),3)
inter_pct<-round(mean(cis_trans_pct$interspecies),3)
library(tidyr)
cis_trans_pct_l<-gather(cis_trans_pct, key = "type_of_interaction", value = "fraction_of_unique_interaction")
cis_pct/trans_pct
pdf('293_3t3_cis_trans_interspecies_QC.pdf',width=5,height=4)
ggplot(cis_trans_pct_l, aes(x = type_of_interaction, y = fraction_of_unique_interaction, fill = type_of_interaction)) +
  geom_violin() + geom_boxplot(width=0.1,outlier.shape = NA)+
  theme_classic()+scale_fill_manual(values=c(mypal[1],mypal[2],mypal[4])) +
  labs(x = "Column Name", y = "Value", title = paste0('cis: ',cis_pct, '  trans:',trans_pct,'  interspecies:',inter_pct))
dev.off()
################## cis short, cis long, trans, interspecies ##################
















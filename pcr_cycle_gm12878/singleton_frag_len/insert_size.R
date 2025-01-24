library(patchwork)
library(gridExtra)
library(gtable)
library(ggplot2)
library(tibble)
library(kneedle)
library(ggsci)

library(purrr)
library(dplyr)

mypal <- pal_npg("nrc", alpha = 0.7)(9)

################## gc cov 20241016 ##################
dipC<-read.table('dipC_insert_size.xls',header=TRUE)
cycle15<-read.table('cycle_15_insert_size.xls',header=TRUE)
cycle18<-read.table('cycle_18_insert_size.xls',header=TRUE)
cycle20<-read.table('cycle_20_insert_size.xls',header=TRUE)

stack<-rbind.data.frame(dipC,cycle15,cycle18,cycle20)

pdf('insert_size.pdf',width=6,height=3)
ggplot(stack, aes(x=insert_size,color=cycle_number))+ geom_density(lwd=1)+theme_classic()+
  scale_color_manual(values=c(mypal[2],mypal[1],mypal[4],mypal[9]))+
  labs(title='singleton insert_size, N=10')
dev.off()

ggplot(stack, aes(x=log10(insert_size),color=cycle_number))+ geom_density(lwd=1)+theme_classic()+
  scale_color_manual(values=c(mypal[2],mypal[1],mypal[4],mypal[9]))+
  labs(title='singleton insert_size, N=10')

svg('insert_size.svg',width=6,height=3)
ggplot(stack, aes(x=insert_size,color=cycle_number))+ geom_density(lwd=1)+theme_classic()+
  scale_color_manual(values=c(mypal[2],mypal[1],mypal[4],mypal[9]))+
  labs(title='singleton insert_size, N=10')
dev.off()

################## gc cov 20241016 ##################



stack<-rbind.data.frame(cbind.data.frame(sample='dipC',chrM_chr_count=dipC$NM_count),
                        cbind.data.frame(sample='cycle15',chrM_chr_count=cycle15$NM_count),
                        cbind.data.frame(sample='cycle18',chrM_chr_count=cycle18$NM_count),
                        cbind.data.frame(sample='cycle20',chrM_chr_count=cycle20$NM_count))
stack<-stack[stack$chrM_chr_count>1,]
pdf('mito_noise_dist.pdf',width=5,height=5)
ggplot(stack, aes(x=chrM_chr_count,color=sample))+ geom_density(lwd=1)+ theme_classic()+scale_color_manual(values=c(mypal[2],mypal[1],mypal[4],mypal[9]))+
  labs(title='chrM vs chr1-22 X Y pair, 100000 rmd validpairs
       cycle15 N=283 , cycle18 N=95, cycle20 N=104')
dev.off()

svg('mito_noise_dist.svg',width=5,height=5)
ggplot(stack, aes(x=chrM_chr_count,color=sample))+ geom_density(lwd=1)+ theme_classic()+scale_color_manual(values=c(mypal[2],mypal[1],mypal[4],mypal[9]))+
  labs(title='gm12878 chrM vs chr1-22 X Y pair, 100000 rmd validpairs
       cycle15 N=283 , cycle18 N=95, cycle20 N=104')
dev.off()
# +scale_color_manual(values=c(mypal[1],mypal[2],mypal[4],mypal[3],mypal[7],mypal[8]))
hichew_mito
hichew_mito[nrow(hichew_mito):1,][1:8]

################## RF dup dist ##################
dupLevel_func<-function(sample_name,in_dir){
  dup<-read.table(paste(in_dir,'/',sample_name,'_uniq_contact_per_RF.xls', sep=""),header=FALSE)
  dupLevel_df<-cbind.data.frame(dupLevel=c(1:5),pct=dup[nrow(dup):1,3][1:5])
  #dupLevel_df$pct[is.na(dupLevel_df$pct)]<-0
  #dupLevel_df$sample<-sample_name
  dupLevel_list<<- append(dupLevel_list,list(dupLevel_df)) 
}

in_dir<-'hichew_uniq_contact_per_RF'
sample_list<-read.table(paste0('hichew_sample.list'))
suffix<-"_uniq_contact_per_RF.xls"
file_name<-paste(sample_list$V1,suffix , sep="")

hicf_raw<-as.list(sample_list$V1)
names(hicf_raw)<-as.list(sample_list$V1)
hicf<-hicf_raw[1:499]

label<-paste0('hichew_uniq_contact_per_RF, n=',length(hicf))
dupLevel_list <- list()
hic_data<- map2(hicf,in_dir,dupLevel_func)
dupLevel_stack <- bind_rows(dupLevel_list)
p2<-ggplot(data = dupLevel_stack, aes(x = dupLevel, y = pct)) + geom_line(alpha=1)+labs(title=label)+theme_classic()
p2
hichew_dupLevel_stack<-dupLevel_stack
hichew_dupLevel_stack$sample<-'hichew'



in_dir<-'dip_C_uniq_contact_per_RF'
sample_list<-read.table(paste0('dip_C_sample.list'))
suffix<-"_uniq_contact_per_RF.xls"
file_name<-paste(sample_list$V1,suffix , sep="")

hicf_raw<-as.list(sample_list$V1)
names(hicf_raw)<-as.list(sample_list$V1)
hicf<-hicf_raw[1:561]

label<-paste0('dip_C_uniq_contact_per_RF, n=',length(hicf))
dupLevel_list <- list()
hic_data<- map2(hicf,in_dir,dupLevel_func)
dupLevel_stack <- bind_rows(dupLevel_list)
p3<-ggplot(data = dupLevel_stack, aes(x = dupLevel, y = pct)) + geom_line(alpha=1)+labs(title=label)+theme_classic()
p3

plotlist<-list(p2,p3)

pdf('hichew_dip_C_uniq_contact_per_RF.pdf',width = 10,height = 5)
cowplot::plot_grid(plotlist = plotlist, ncol = 2)
dev.off()

dip_C_dupLevel_stack<-dupLevel_stack
dip_C_dupLevel_stack$sample<-'dip_C'
dip_C_hichew_dupLevel_stack<-rbind.data.frame(hichew_dupLevel_stack,dip_C_dupLevel_stack)

pdf('hichew_dip_C_uniq_contact_per_RF_stack.pdf',width = 5,height = 5)
ggplot(data = dip_C_hichew_dupLevel_stack, aes(x = dupLevel, y = pct, color=sample)) + geom_line(alpha=1)+labs(title='hichew_dip_C_uniq_contact_per_RF')+theme_classic()
dev.off()


################## knee cut off test ##################
cell_rank_df<-cbind.data.frame(rmd_valid_pair=sort(hichew_1203$valid_interaction_rmdup,decreasing = TRUE),
                               cell_rank=c(1:length(hichew_1203$valid_interaction_rmdup)))
rank_df_y_axis<-sort(hichew_1203$valid_interaction_rmdup,decreasing = TRUE)
rank_df_y_axis[10:length(rank_df_y_axis)]
knee1<-kneedle(c(1:length(hichew_1203$valid_interaction_rmdup)),rank_df_y_axis, 
               concave = FALSE, decreasing = TRUE,sensitivity = 20)
knee1
start2<-knee1[1]+100
knee2<-kneedle(c(start2:length(hichew_1203$valid_interaction_rmdup)),rank_df_y_axis[start2:length(rank_df_y_axis)] ,
               concave = FALSE, decreasing = TRUE,sensitivity = 20)
ggplot(cell_rank_df, aes(x=cell_rank,y=rmd_valid_pair))+ 
  geom_line()+labs(title=knee2[1])+ theme_classic()+geom_vline(xintercept= knee2[1], linetype='dotted')
################## knee cut off test ##################


################## 3D ##################
hc_293_3t3<-read.table('../../raw_table/sc_hic_pro_metrics/hichew_293_3t3_hic_pro_metrics.xls',header=TRUE)
library(plotly)
plot_ly(
  hc_293_3t3, x = ~log10(hg19), y = ~log10(mm10), z = ~log10(hg19_mm10), size=0.05
) %>%
  add_markers() %>%
  layout(
    scene = list(xaxis = list(title = 'hg19'),
                 yaxis = list(title = 'mm10'),
                 zaxis = list(title = 'hg19_mm10')))
################## 3D ##################


################## knee cut off human mouse ##################
cell_rank_df<-cbind.data.frame(rmd_valid_pair=sort(hc_293_3t3$valid_interaction_rmdup,decreasing = TRUE),
                               cell_rank=c(1:length(hc_293_3t3$valid_interaction_rmdup)))
rank_df_y_axis<-sort(hc_293_3t3$valid_interaction_rmdup,decreasing = TRUE)
rank_df_y_axis[10:length(rank_df_y_axis)]
knee1<-kneedle(c(1:length(hc_293_3t3$valid_interaction_rmdup)),rank_df_y_axis, 
               concave = FALSE, decreasing = TRUE,sensitivity = 20)
knee1
start2<-knee1[1]+100
knee2<-kneedle(c(start2:length(hc_293_3t3$valid_interaction_rmdup)),log10(rank_df_y_axis[start2:length(rank_df_y_axis)]) ,
               concave = FALSE, decreasing = TRUE,sensitivity = 20)
pdf('293_3t3_cell_rank.pdf',width=4,height=4)
ggplot(cell_rank_df, aes(x=cell_rank,y=log10(rmd_valid_pair)))+ xlim(0,1900)+
  geom_line()+labs(title=paste0('N = 1900, cut off: ',knee2[1]))+ theme_classic()+geom_vline(xintercept= knee2[1], linetype='dotted')
dev.off()
################## knee cut off human mouse ##################


################## empty doublet human mouse ##################
sorted_hc_293_3t3 <- hc_293_3t3[order(hc_293_3t3$valid_interaction_rmdup,decreasing=TRUE), ]
rownames(sorted_hc_293_3t3)<-c(1:length(rownames(sorted_hc_293_3t3)))
clean_hc_293_3t3 <- sorted_hc_293_3t3[c(1:knee2[1]),]
subset_hc_293_3t3 <- sorted_hc_293_3t3[c(1:1500),]

sorted_hc_293_3t3$group<-'empty'
knee2
knee2_thres<-10^knee2[2]
sorted_hc_293_3t3[sorted_hc_293_3t3$valid_interaction_rmdup>knee2_thres&sorted_hc_293_3t3$hg19/sorted_hc_293_3t3$mm10 > 8,length(colnames(sorted_hc_293_3t3))] <-'hg19'
sorted_hc_293_3t3[sorted_hc_293_3t3$valid_interaction_rmdup>knee2_thres&sorted_hc_293_3t3$hg19/sorted_hc_293_3t3$mm10 < 0.5,length(colnames(sorted_hc_293_3t3))] <-'mm10'
sorted_hc_293_3t3[sorted_hc_293_3t3$valid_interaction_rmdup>knee2_thres&sorted_hc_293_3t3$hg19/sorted_hc_293_3t3$mm10 < 8 & sorted_hc_293_3t3$hg19/sorted_hc_293_3t3$mm10 > 0.5,length(colnames(sorted_hc_293_3t3))] <-'doublet'

hg19_num<-length(which(sorted_hc_293_3t3$group=="hg19"))
mm10_num<-length(which(sorted_hc_293_3t3$group=="mm10"))
doublet_num<-length(which(sorted_hc_293_3t3$group=="doublet"))
empty_num<-length(which(sorted_hc_293_3t3$group=="empty"))

pdf('293_3t3_doublet.pdf',width=5,height=4)
ggplot(sorted_hc_293_3t3, aes(x=log10(hg19),y=log10(mm10),color=group)) + geom_point(alpha=0.5,stroke=NA) +xlim(1,5)+ylim(1,5)+
  theme_classic() +scale_color_manual(values=c(mypal[1],'grey',mypal[3],mypal[4]))+
  labs(title=paste0('loaded cell: 1200, demultiplexed barcode: top 1900
hg19: ',hg19_num,'  mm10:',mm10_num,'  empty:',empty_num,'  doublet:',doublet_num))
dev.off()

write.table(sorted_hc_293_3t3,file = 'sorted_hc_293_3t3.xls',quote = FALSE,sep="\t",row.names = FALSE)

################## empty doublet human mouse ##################


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

pdf('293_3t3_cis_trans_interspecies_QC.pdf',width=5,height=4)
ggplot(cis_trans_pct_l, aes(x = type_of_interaction, y = fraction_of_unique_interaction, fill = type_of_interaction)) +
  geom_violin() + geom_boxplot(width=0.1,outlier.shape = NA)+
  theme_classic()+scale_fill_manual(values=c(mypal[1],mypal[2],mypal[4])) +
  labs(x = "Column Name", y = "Value", title = paste0('cis: ',cis_pct, '  trans:',trans_pct,'  interspecies:',inter_pct))
dev.off()
################## cis short, cis long, trans, interspecies ##################
















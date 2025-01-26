library(ggplot2)
library(ggsci)
mypal <- pal_npg("nrc", alpha = 0.7)(9)
library(kneedle)
library(ggsci)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
in_metrics=args[1]
out_path_to_file=args[2]
label_hichew_test=args[3]
knee1_sensitivity=as.integer(args[4])
knee2_sensitivity=as.integer(args[5])
knee1_knee2_interval=as.integer(args[6])

# raw_metrics<-read.table('snHiChew_HEK293Twt_hic_pro_metrics.xls',header = TRUE)
# label_hichew_test<-'snHiChew_HEK293Twt'
raw_metrics<-read.table(in_metrics,header = TRUE)

################## knee cut off ##################
cell_rank_df<-cbind.data.frame(rmd_valid_pair=sort(raw_metrics$valid_interaction_rmdup,decreasing = TRUE),
                               cell_rank=c(1:length(raw_metrics$valid_interaction_rmdup)))
rank_df_y_axis<-sort(raw_metrics$valid_interaction_rmdup,decreasing = TRUE)
rank_df_y_axis[10:length(rank_df_y_axis)]
knee1<-kneedle(c(1:length(raw_metrics$valid_interaction_rmdup)),rank_df_y_axis, 
               concave = FALSE, decreasing = TRUE,sensitivity = knee1_sensitivity)
knee1
start2<-knee1[1]+knee1_knee2_interval
knee2<-kneedle(c(start2:length(raw_metrics$valid_interaction_rmdup)),rank_df_y_axis[start2:length(rank_df_y_axis)],
               concave = FALSE, decreasing = TRUE,sensitivity = knee2_sensitivity)
knee2
pdf(paste0(label_hichew_test,'_cell_rank_log10.pdf'),width=4,height=4)
ggplot(cell_rank_df, aes(x=cell_rank,y=log10(rmd_valid_pair)))+ 
  geom_line()+labs(title=paste0('cut-off 1 :',knee1[1],'\ncut-off 2 :',knee2[1]))+ 
  theme_classic()+geom_vline(xintercept= c(knee1[1],knee2[1]), linetype='dotted')
dev.off()

svg(paste0(label_hichew_test,'_cell_rank_log10.svg'),width=4,height=4)
ggplot(cell_rank_df, aes(x=cell_rank,y=log10(rmd_valid_pair)))+
  geom_line()+labs(title=paste0('cut-off 1 :',knee1[1],'\ncut-off 2 :',knee2[1]))+
  theme_classic()+geom_vline(xintercept= c(knee1[1],knee2[1]), linetype='dotted')
dev.off()

png(paste0(label_hichew_test,'_cell_rank_log10.png'))
ggplot(cell_rank_df, aes(x=cell_rank,y=log10(rmd_valid_pair)))+
  geom_line()+labs(title=paste0('cut-off 1 :',knee1[1],'\ncut-off 2 :',knee2[1]))+
  theme_classic()+geom_vline(xintercept= c(knee1[1],knee2[1]), linetype='dotted')
dev.off()

pdf(paste0(label_hichew_test,'_cell_rank.pdf'),width=4,height=4)
ggplot(cell_rank_df, aes(x=cell_rank,y=rmd_valid_pair))+ 
  geom_line()+labs(title=paste0('cut-off 1 :',knee1[1],'\ncut-off 2 :',knee2[1]))+ 
  theme_classic()+geom_vline(xintercept= c(knee1[1],knee2[1]), linetype='dotted')
dev.off()

svg(paste0(label_hichew_test,'_cell_rank.svg'),width=4,height=4)
ggplot(cell_rank_df, aes(x=cell_rank,y=rmd_valid_pair))+
  geom_line()+labs(title=paste0('cut-off 1 :',knee1[1],'\ncut-off 2 :',knee2[1]))+
  theme_classic()+geom_vline(xintercept= c(knee1[1],knee2[1]), linetype='dotted')
dev.off()

png(paste0(label_hichew_test,'_cell_rank.png'))
ggplot(cell_rank_df, aes(x=cell_rank,y=rmd_valid_pair))+
  geom_line()+labs(title=paste0('cut-off 1 :',knee1[1],'\ncut-off 2 :',knee2[1]))+
  theme_classic()+geom_vline(xintercept= c(knee1[1],knee2[1]), linetype='dotted')
dev.off()

sorted_metrics <- raw_metrics[order(raw_metrics$valid_interaction_rmdup,decreasing=TRUE), ]
rownames(sorted_metrics)<-c(1:length(rownames(sorted_metrics)))
sorted_cut_hichew_test <- sorted_metrics[c(1:knee2[1]),]

# write.table(sorted_cut_hichew_test,file = paste0('filtered_',label_hichew_test,'_hic_pro_metrics.xls'),quote = FALSE,sep="\t",row.names = FALSE)
write.table(sorted_cut_hichew_test,file = paste0(out_path_to_file),quote = FALSE,sep="\t",row.names = FALSE)

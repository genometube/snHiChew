library(ggplot2)
library(purrr)
library(dplyr)
library(ggsci)
mypal <- pal_npg("nrc", alpha = 0.7)(9)

dip_C<-read.table('dipC_50cell_20240923_downsample_rmdup_validpair.xls',header=FALSE)
# dip_C$sample<-'dip_C_50cell'
colnames(dip_C)<-c('readNum_per_cell','rmdup_valid_pair_dip_C_50cell_1','dup_dip_C_50cell_1')

hichew<-read.table('hichew_50cell_20240923_downsample_rmdup_validpair.xls',header=FALSE)
# hichew$sample<-'hichew_50cell'
colnames(hichew)<-c('readNum_per_cell','rmdup_valid_pair_hichew_50cell_1','dup_hichew_50cell_1')

dip_C_<-read.table('PCR_combine_160M_downsample_rmdup_validpair_dup.xls',header=FALSE)
#dip_C_$sample<-'dip_C_50cell'
colnames(dip_C_)<-c('readNum_per_cell','rmdup_valid_pair_dip_C_50cell_2','dup_dip_C_50cell_2')

hichew_<-read.table('hichew_combine_160M_downsample_rmdup_validpair_dup.xls',header=FALSE)
#hichew_$sample<-'hichew_50cell'
colnames(hichew_)<-c('readNum_per_cell','rmdup_valid_pair_hichew_50cell_2','dup_hichew_50cell_2')


dip_C_merge<-merge(dip_C,dip_C_,by=c('readNum_per_cell'))
dip_C_merge$mean_validpair<-(dip_C_merge$rmdup_valid_pair_dip_C_50cell_1 + dip_C_merge$rmdup_valid_pair_dip_C_50cell_2)/2
row_sd <- function(x) {apply(x, 1, sd)}
# Calculate row-wise standard deviation for columns A and B
dip_C_merge$validpair_SD <- row_sd(dip_C_merge[, c("rmdup_valid_pair_dip_C_50cell_1", "rmdup_valid_pair_dip_C_50cell_2")])
dip_C_merge$mean_dup<-(dip_C_merge$dup_dip_C_50cell_1 + dip_C_merge$dup_dip_C_50cell_2)/2
# Calculate row-wise standard deviation for columns A and B
dip_C_merge$dup_SD <- row_sd(dip_C_merge[, c("dup_dip_C_50cell_1", "dup_dip_C_50cell_2")])


hichew_merge<-merge(hichew,hichew_,by=c('readNum_per_cell'))
hichew_merge$mean_validpair<-(hichew_merge$rmdup_valid_pair_hichew_50cell_1 + hichew_merge$rmdup_valid_pair_hichew_50cell_2)/2
row_sd <- function(x) {apply(x, 1, sd)}
# Calculate row-wise standard deviation for columns A and B
hichew_merge$validpair_SD <- row_sd(hichew_merge[, c("rmdup_valid_pair_hichew_50cell_1", "rmdup_valid_pair_hichew_50cell_2")])
hichew_merge$mean_dup<-(hichew_merge$dup_hichew_50cell_1 + hichew_merge$dup_hichew_50cell_2)/2
# Calculate row-wise standard deviation for columns A and B
hichew_merge$dup_SD <- row_sd(hichew_merge[, c("dup_hichew_50cell_1", "dup_hichew_50cell_2")])

hichew_merge$sample<-'hichew_50cell'
dip_C_merge$sample<-'dip_C_50cell'

stack_df<-rbind.data.frame(hichew_merge[,c(1,6:10)],dip_C_merge[,c(1,6:10)])
p2<-ggplot(stack_df, aes(x = readNum_per_cell, y = mean_validpair,group=sample,color=sample))+ 
  geom_line(alpha=0.8,lwd=1,aes(linetype=sample))+ geom_point()+scale_linetype_manual(values=c(3,1,5))+
  geom_pointrange( aes(x=readNum_per_cell, ymin=mean_validpair-validpair_SD, ymax=mean_validpair+validpair_SD),
                   alpha=0.6, size=0.5)+
  #scale_x_continuous(breaks = c(1e5,1e6,1e7))+
  theme_classic()+scale_color_manual(values=c(mypal[1],mypal[2],mypal[4]))+
  labs(title='rmdvalidpair')
p2


p3<-ggplot(stack_df, aes(x = readNum_per_cell, y = mean_dup,group=sample,color=sample))+ 
  geom_line(alpha=0.8,lwd=1,aes(linetype=sample))+ geom_point()+scale_linetype_manual(values=c(3,1,5))+
  geom_pointrange( aes(x=readNum_per_cell, ymin=mean_dup-dup_SD, ymax=mean_dup+dup_SD),
                   alpha=0.6, size=0.5)+
  #scale_x_continuous(breaks = c(1e5,1e6,1e7))+
  theme_classic()+scale_color_manual(values=c(mypal[1],mypal[2],mypal[4]))+
  labs(title='dup')
p3

library(gridExtra)

pdf('rmd_validpair_saturation_curve.pdf',height=6,width=6)
grid.arrange(p2,p3,nrow=2,ncol=1)
dev.off()

svg('rmd_validpair_saturation_curve.svg',height=6,width=6)
grid.arrange(p2,p3,nrow=2,ncol=1)
dev.off()

write.table(stack_df,file = 'with_error_bar_metrics.xls',quote = FALSE,sep="\t",row.names = FALSE)


######################### 20240505


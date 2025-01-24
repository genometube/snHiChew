library(ggplot2)
library(purrr)
library(dplyr)
library(ggsci)
mypal <- pal_npg("nrc", alpha = 0.7)(9)

############ single run ############ 
dup_readnum<-read.table('nagano2017_downsample_dup.xls',header=FALSE)
dup_readnum$sample<-'dup_readnum'
colnames(dup_readnum)<-c('readNum','rmdup_valid_pair','sd','sample')
dup_readnum$readNum<-log10(dup_readnum$readNum)
dup_readnum$rmdup_valid_pair<-log10(dup_readnum$rmdup_valid_pair)
dup_readnum$sd<-log10(dup_readnum$sd)
ggplot(dup_readnum, aes(x = readNum, y = rmdup_valid_pair))+ 
  geom_line()+ geom_point()+scale_x_continuous(breaks = c(1e5,1e6,1e7))+
  geom_pointrange(aes(x=readNum, ymin=rmdup_valid_pair-sd, ymax=rmdup_valid_pair+sd),
                 width=0.4, colour="orange", alpha=1, size=1)+
  theme_classic()

############ stack plot ############ 

nagano2017<-read.table('nagano2017_downsample_dup.xls',header=FALSE)
nagano2017$sample<-'nagano2017'
colnames(nagano2017)<-c('readNum_per_cell','rmdvalidpair_ratio','sd','sample')
nagano2017$rmdvalidpair_ratio<- 1-nagano2017$rmdvalidpair_ratio

dip_C<-read.table('dip_C_downsample_dup.xls',header=FALSE)
dip_C$sample<-'dip_C'
colnames(dip_C)<-c('readNum_per_cell','rmdvalidpair_ratio','sd','sample')
dip_C$rmdvalidpair_ratio<-1-dip_C$rmdvalidpair_ratio
#scCARE<-read.table('scCARE_seq_downsample_dup.xls',header=FALSE)
#scCARE$sample<-'scCARE'
#colnames(scCARE)<-c('readNum_per_cell','mean_rmdup_valid_pair','sd','sample')

hichew<-read.table('hichew_20_lanes_downsample_dup.xls',header=FALSE)
hichew$sample<-'hichew'
colnames(hichew)<-c('readNum_per_cell','rmdvalidpair_ratio','sd','sample')
hichew$rmdvalidpair_ratio<- 1-hichew$rmdvalidpair_ratio

stack_df<-rbind.data.frame(nagano2017,dip_C,hichew)
#stack_df$mean_rmdup_valid_pair<-log2(stack_df$mean_rmdup_valid_pair)
#stack_df$sd<-log10(stack_df$sd)
1321651/214280

pdf('rmd_validpair_saturation_curve.pdf',height=5,width=8)
ggplot(stack_df, aes(x = readNum_per_cell, y = rmdvalidpair_ratio,group=sample,color=sample))+ 
  geom_line()+ geom_point()+
  geom_pointrange( aes(x=readNum_per_cell, ymin=rmdvalidpair_ratio-sd, ymax=rmdvalidpair_ratio+sd),
                     alpha=0.3, size=0.3)+
  #scale_x_continuous(breaks = c(1e5,1e6,1e7))+
  theme_classic()+scale_color_manual(values=c(mypal[1],mypal[2],mypal[4]))+
  labs(title='rmdvalidpair / validpair',subtitle='gradient: 25000 50000 100000 200000 400000 
       800000 1600000 3200000 6400000 12800000\nnagano2017 515 cells \ndipC 1071 cells\nhichew 1178 cells')+xlim(0,15000000)
dev.off()

'#E64B35B2'
'#4DBBD5B2'
'#3C5488B2'

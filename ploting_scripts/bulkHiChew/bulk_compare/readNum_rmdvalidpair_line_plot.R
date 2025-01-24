library(ggplot2)
library(purrr)
library(dplyr)
library(ggsci)
mypal <- pal_npg("nrc", alpha = 0.7)(9)

dip_C<-read.table('dipC_50cell_20240923_downsample_rmdup_validpair.xls',header=FALSE)
dip_C$sample<-'dip_C_50cell'
colnames(dip_C)<-c('readNum_per_cell','mean_rmdup_valid_pair','sample')

hichew<-read.table('hichew_50cell_20240923_downsample_rmdup_validpair.xls',header=FALSE)
hichew$sample<-'hichew_50cell'
colnames(hichew)<-c('readNum_per_cell','mean_rmdup_valid_pair','sample')


stack_df<-rbind.data.frame(dip_C,hichew)
#stack_df<-stack_df[stack_df$readNum_per_cell<6400000,]

ggplot(stack_df, aes(x = readNum_per_cell, y = mean_rmdup_valid_pair,group=sample,color=sample))+ 
  geom_line(alpha=0.8,lwd=1,aes(linetype=sample))+scale_linetype_manual(values=c(3,1,5))+ 
  geom_point()+
  #scale_x_continuous(breaks = c(1e5,1e6,1e7))+
  theme_classic()+scale_color_manual(values=c(mypal[1],mypal[2]))+
  #scale_y_continuous(labels = scales::scientific_format())+ theme(axis.title.y = element_blank())+
  labs(title='rmd_validpair',subtitle='gradient: 781250 1562500 3125000 6250000 \n12.5M 25M 50M 100M 200M 300M 350M\nbulk dipC 50 cells\nbulk hichew 50 cells')

#pdf('rmdvalidpair_saturation_curve_log10.pdf',height=5,width=8)

p1<-ggplot(stack_df, aes(x = readNum_per_cell, y = mean_rmdup_valid_pair,group=sample,color=sample))+ 
  geom_line(alpha=0.8,lwd=1,aes(linetype=sample))+scale_linetype_manual(values=c(3,1,5))+ 
  geom_point()+
  #scale_x_continuous(breaks = c(1e5,1e6,1e7))+
  theme_classic()+scale_color_manual(values=c(mypal[1],mypal[2]))+
  #scale_y_continuous(labels = scales::scientific_format())+ theme(axis.title.y = element_blank())+
  labs(title='rmd_validpair',subtitle='gradient: 781250 1562500 3125000 6250000 \n12.5M 25M 50M 100M 200M 300M 350M\nbulk dipC 50 cells\nbulk hichew 50 cells')

p1
#dev.off()


######################### 20240505
dip_C<-read.table('hichew_50cell_20240923_downsample_dup.xls',header=FALSE)
dip_C$sample<-'hichew_50cell'
colnames(dip_C)<-c('readNum_per_cell','rmdvalidpair_dup','sample')

hichew<-read.table('dipC_50cell_20240923_downsample_dup.xls',header=FALSE)
hichew$sample<-'dipC_50cell'
colnames(hichew)<-c('readNum_per_cell','rmdvalidpair_dup','sample')

stack_df<-rbind.data.frame(dip_C,hichew)
#stack_df$mean_rmdup_valid_pair<-log2(stack_df$mean_rmdup_valid_pair)
#stack_df$sd<-log10(stack_df$sd)
#1321651/214280

#pdf('rmd_validpair_saturation_curve.pdf',height=5,width=8)
p2<-ggplot(stack_df, aes(x = readNum_per_cell, y = rmdvalidpair_dup,group=sample,color=sample))+ 
  geom_line(alpha=0.8,lwd=1,aes(linetype=sample))+ geom_point()+scale_linetype_manual(values=c(3,1,5))+
  theme_classic()+scale_color_manual(values=c(mypal[1],mypal[2],mypal[4]))+
  labs(title='rmdvalidpair_dup',subtitle='gradient: 781250 1562500 3125000 6250000 \n12.5M 25M 50M 100M 200M 300M 350M\nbulk dipC 50 cells\nbulk hichew 50 cells')

p2
#dev.off()

library(gridExtra)

pdf('rmd_validpair_saturation_curve.pdf',height=8,width=6)
grid.arrange(p1,p2,nrow=2,ncol=1)
dev.off()

svg('rmd_validpair_saturation_curve.svg',height=8,width=6)
grid.arrange(p1,p2,nrow=2,ncol=1)
dev.off()









############ single run ############ 
rmdup_validpair_readnum<-read.table('nagano2017_downsample_rmdup_validpair_log10.xls',header=FALSE)
rmdup_validpair_readnum$sample<-'rmdup_validpair_readnum'
colnames(rmdup_validpair_readnum)<-c('readNum','rmdup_valid_pair','sd','sample')
rmdup_validpair_readnum$readNum<-log10(rmdup_validpair_readnum$readNum)
rmdup_validpair_readnum$rmdup_valid_pair<-log10(rmdup_validpair_readnum$rmdup_valid_pair)
rmdup_validpair_readnum$sd<-log10(rmdup_validpair_readnum$sd)
ggplot(rmdup_validpair_readnum, aes(x = readNum, y = rmdup_valid_pair))+ 
  geom_line()+ geom_point()+scale_x_continuous(breaks = c(1e5,1e6,1e7))+
  geom_pointrange( aes(x=readNum, ymin=rmdup_valid_pair-sd, ymax=rmdup_valid_pair+sd),
                 width=0.4, colour="orange", alpha=1, size=1)+
  theme_classic()

############ stack plot ############ 

nagano2017<-read.table('nagano2017_downsample_rmdup_validpair_log.xls',header=FALSE)
nagano2017$sample<-'nagano2017'
colnames(nagano2017)<-c('readNum_per_cell','mean_rmdup_valid_pair','sd','sample')

dip_C<-read.table('dip_C_downsample_rmdup_validpair_log.xls',header=FALSE)
dip_C$sample<-'dip_C'
colnames(dip_C)<-c('readNum_per_cell','mean_rmdup_valid_pair','sd','sample')

hichew<-read.table('hichew_20_lanes_downsample_rmdvalidpair_log.xls',header=FALSE)
hichew$sample<-'hichew'
colnames(hichew)<-c('readNum_per_cell','mean_rmdup_valid_pair','sd','sample')

stack_df<-rbind.data.frame(nagano2017,dip_C,hichew)

stack_df<-stack_df[stack_df$readNum_per_cell<6400000,]

pdf('rmdvalidpair_saturation_curve_log10.pdf',height=5,width=8)

ggplot(stack_df, aes(x = readNum_per_cell, y = mean_rmdup_valid_pair,group=sample,color=sample))+ 
  geom_line()+ geom_point()+
  geom_pointrange( aes(x=readNum_per_cell, ymin=mean_rmdup_valid_pair-sd, ymax=mean_rmdup_valid_pair+sd),
                     alpha=0.3, size=0.3)+
  #scale_x_continuous(breaks = c(1e5,1e6,1e7))+
  theme_classic()+scale_color_manual(values=c(mypal[1],mypal[2],mypal[4]))+
  labs(subtitle='gradient: 25000 50000 100000 200000 400000 
       800000 1600000 3200000 6400000 12800000\nnagano2017 515 cells \ndipC 1071 cells\nhichew 1178 cells')+xlim(0,15000000)
dev.off()

'#E64B35B2'
'#4DBBD5B2'
'#3C5488B2'

ggplot(stack_df, aes(x = readNum_per_cell, y = mean_rmdup_valid_pair,group=sample,color=sample))+ 
  geom_line()+ geom_point()+
  geom_pointrange( aes(x=readNum_per_cell, ymin=mean_rmdup_valid_pair-sd, ymax=mean_rmdup_valid_pair+sd),
                   width=0.4,  alpha=0.3, size=0.3)+
  scale_x_continuous(breaks = c(1e5,1e6,1e7))+
  theme_classic()+
  labs(title='gradient: 25000 50000 100000 200000 400000 
       800000 1600000 3200000 6400000 12800000')+xlim(20000,1600000)

ggplot(stack_df, aes(x = readNum_per_cell, y = mean_rmdup_valid_pair,group=sample,color=sample))+ 
  geom_line()+ geom_point()+scale_x_continuous(breaks = c(1e5,1e6,1e7))+
  theme_classic()+
  labs(title='gradient: 25000 50000 100000 200000 400000 
       800000 1600000 3200000 6400000 12800000')+xlim(20000,3200000)





######################### 20240505

######################### 20240830

nagano2017<-read.table('hichew_downsample_rmdup_validpair.xls',header=FALSE)
nagano2017$sample<-'hichew'
colnames(nagano2017)<-c('readNum_per_cell','mean_rmdup_valid_pair','sd','sample')

dip_C<-read.table('dipC_downsample_rmdup_validpair.xls',header=FALSE)
dip_C$sample<-'dip_C'
colnames(dip_C)<-c('readNum_per_cell','mean_rmdup_valid_pair','sd','sample')

hichew<-read.table('hichew_Rep2_downsample_rmdup_validpair.xls',header=FALSE)
hichew$sample<-'hichew_rep2'
colnames(hichew)<-c('readNum_per_cell','mean_rmdup_valid_pair','sd','sample')


stack_df<-rbind.data.frame(nagano2017,dip_C,hichew)
#stack_df<-stack_df[stack_df$readNum_per_cell<6400000,]

ggplot(stack_df, aes(x = readNum_per_cell, y = mean_rmdup_valid_pair,group=sample,color=sample))+ 
  geom_line(alpha=0.8,lwd=1,aes(linetype=sample))+scale_linetype_manual(values=c(3,1,5))+ 
  geom_point()+  geom_pointrange( aes(x=readNum_per_cell, ymin=mean_rmdup_valid_pair-sd, ymax=mean_rmdup_valid_pair+sd), alpha=0.2, size=0.5)+
  #scale_x_continuous(breaks = c(1e5,1e6,1e7))+
  theme_classic()+scale_color_manual(values=c(mypal[1],mypal[2],mypal[4]))+
  #scale_y_continuous(labels = scales::scientific_format())+ theme(axis.title.y = element_blank())+
  labs(title='rmd_validpair',subtitle='gradient: 25000 50000 100000 200000 400000 800000\nnagano2017 400 cells \ndipC 400 cells\nhichew 400 cells')

#pdf('rmdvalidpair_saturation_curve_log10.pdf',height=5,width=8)

p1<-ggplot(stack_df, aes(x = readNum_per_cell, y = mean_rmdup_valid_pair,group=sample,color=sample))+ 
  geom_line(alpha=0.8,lwd=1,aes(linetype=sample))+scale_linetype_manual(values=c(3,1,5))+ 
  geom_point()+  geom_pointrange( aes(x=readNum_per_cell, ymin=mean_rmdup_valid_pair-sd, ymax=mean_rmdup_valid_pair+sd), alpha=0.2, size=0.5)+
  #scale_x_continuous(breaks = c(1e5,1e6,1e7))+
  theme_classic()+scale_color_manual(values=c(mypal[1],mypal[2],mypal[4]))+
 scale_y_continuous(labels = scales::scientific_format())+ theme(axis.title.y = element_blank())+
  labs(title='rmd_validpair',subtitle='gradient: 25000 50000 100000 200000 400000 
       800000\nnagano2017 8 cells \ndipC 8 cells\nhichew 8 cells')

p1
#dev.off()



######################### 20240505
nagano2017<-read.table('hichew_downsample_dup.xls',header=FALSE)
nagano2017$sample<-'hichew'
colnames(nagano2017)<-c('readNum_per_cell','rmdvalidpair_dup','sd','sample')

dip_C<-read.table('hichew_Rep2_downsample_dup.xls',header=FALSE)
dip_C$sample<-'hichew_rep2'
colnames(dip_C)<-c('readNum_per_cell','rmdvalidpair_dup','sd','sample')
#scCARE<-read.table('scCARE_seq_downsample_dup.xls',header=FALSE)
#scCARE$sample<-'scCARE'
#colnames(scCARE)<-c('readNum_per_cell','mean_rmdup_valid_pair','sd','sample')

hichew<-read.table('dipC_downsample_dup.xls',header=FALSE)
hichew$sample<-'dipC'
colnames(hichew)<-c('readNum_per_cell','rmdvalidpair_dup','sd','sample')

stack_df<-rbind.data.frame(nagano2017,dip_C,hichew)
#stack_df$mean_rmdup_valid_pair<-log2(stack_df$mean_rmdup_valid_pair)
#stack_df$sd<-log10(stack_df$sd)
#1321651/214280

#pdf('rmd_validpair_saturation_curve.pdf',height=5,width=8)
p2<-ggplot(stack_df, aes(x = readNum_per_cell, y = rmdvalidpair_dup,group=sample,color=sample))+ 
  geom_line(alpha=0.8,lwd=1,aes(linetype=sample))+ geom_point()+scale_linetype_manual(values=c(3,1,5))+
  geom_pointrange( aes(x=readNum_per_cell, ymin=rmdvalidpair_dup-sd, ymax=rmdvalidpair_dup+sd),
                   alpha=0.2, size=0.5)+
  #scale_x_continuous(breaks = c(1e5,1e6,1e7))+
  theme_classic()+scale_color_manual(values=c(mypal[1],mypal[2],mypal[4]))+
  labs(title='rmdvalidpair_dup')
p2
#dev.off()

library(gridExtra)
grid.arrange(p1,p2,nrow=2,ncol=1)

dev.off()

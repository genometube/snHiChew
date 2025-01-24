library(ggplot2)
library(purrr)
library(dplyr)
library(ggsci)
mypal <- pal_npg("nrc", alpha = 0.7)(9)
############################################################################################################ 20241108

hichew_bin_cov<-read.table('snHiChew_HEK293Twt_bin_cov.xls',header=FALSE)
colnames(hichew_bin_cov)<-c('cell','bin','validpair','non_zero_cov_pct')
hichew_bin_cov$validpair<-factor(hichew_bin_cov$validpair, levels = c('validpair_62500','validpair_125000','validpair_250000','validpair_500000'))

#hichew_bin_cov$sample<-'snHiChew_gm12878'
p1<-ggplot(hichew_bin_cov, aes(x = validpair, y = non_zero_cov_pct,group=validpair,color=validpair))+ geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=0.4, alpha=0.2)+facet_wrap(~bin,nrow =1)+scale_color_manual(values=c(mypal[1],mypal[4],mypal[7],mypal[9]))+
  theme_light()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                       strip.background = element_rect(color="black", fill="white"),
                       strip.text.x = element_text(size = 12, color = "black"),
                       axis.text.x = element_text(angle = 20,  hjust=1))+ggtitle(label=c('287cell, validpair downsample, snHiChew_HEK293Twt bin cov'))
#p1


hichew_bin_cov<-read.table('snHiChew_HEK293Twt_bin_cov_PE.xls',header=FALSE)
colnames(hichew_bin_cov)<-c('cell','bin','validpair','non_zero_cov_pct')
hichew_bin_cov$validpair<-factor(hichew_bin_cov$validpair, levels = c('validpair_62500','validpair_125000','validpair_250000','validpair_500000'))

p2<-ggplot(hichew_bin_cov, aes(x = validpair, y = non_zero_cov_pct,group=validpair,color=validpair))+ geom_violin()+
  geom_jitter(size=0.4, alpha=0.2)+facet_wrap(~bin,nrow =1)+scale_color_manual(values=c(mypal[1],mypal[4],mypal[7],mypal[9]))+
  theme_light()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                       strip.background = element_rect(color="black", fill="white"),
                       strip.text.x = element_text(size = 12, color = "black"),
                       axis.text.x = element_text(angle = 20,  hjust=1))+ggtitle(label=c('287cell, validpair downsample, snHiChew_HEK293Twt bin cov promoter enhancer'))
#p2

library(gridExtra)



pdf('snHiChew_HEK293Twt_bin_cov.pdf',width = 15,height = 10)
grid.arrange(p1, p2,nrow=2)
dev.off()

svg('snHiChew_HEK293Twt_bin_cov.svg',width = 15,height = 10)
grid.arrange(p1, p2,nrow=2)
dev.off()

############################################################################################################ 20241108



hichew_bin_cov<-read.table('snHiChew_gm12878_bin_cov.xls',header=FALSE)
colnames(hichew_bin_cov)<-c('cell','bin','non_zero_cov_pct')
hichew_bin_cov$sample<-'snHiChew_gm12878'

dip_C_bin_cov<-read.table('gageseq_gm12878_bin_cov.xls',header=FALSE)
colnames(dip_C_bin_cov)<-c('cell','bin','non_zero_cov_pct')
dip_C_bin_cov$sample<-'gageseq_gm12878'

min_perturbation <- -0.1
max_perturbation <- 0.1
dip_C_bin_cov$non_zero_cov_pct <- dip_C_bin_cov$non_zero_cov_pct + runif(n = nrow(dip_C_bin_cov), min = min_perturbation, max = max_perturbation)

sciHiC_bin_cov<-read.table('sciHiC_gm12878_bin_cov.xls',header=FALSE)
colnames(sciHiC_bin_cov)<-c('cell','bin','non_zero_cov_pct')
sciHiC_bin_cov$sample<-'sciHiC_gm12878'

stack_df<-rbind.data.frame(sciHiC_bin_cov,dip_C_bin_cov,hichew_bin_cov)
stack_df$sample<-factor(stack_df$sample,levels=c('sciHiC_gm12878','gageseq_gm12878','snHiChew_gm12878'))
stack_df<-stack_df[stack_df$bin!=100000,]
stack_df$bin<-as.factor(stack_df$bin)

pdf('top200_bin_cov.pdf',width = 10,height = 5)
p1<-ggplot(stack_df, aes(x = sample, y = non_zero_cov_pct,group=sample,color=sample))+ geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=0.4, alpha=0.9)+facet_wrap(~bin,nrow =1)+scale_color_manual(values=c('#3C5488B2','#E64B35B2','#4DBBD5B2'))+
  theme_light()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                       strip.background = element_rect(color="black", fill="white"),
                       strip.text.x = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(label=c('cell of top200 reads'))
p1
dev.off()

svg('top200_bin_cov.svg',width = 10,height = 5)
p1
dev.off()



hichew_bin_cov<-read.table('snHiChew_gm12878_bin_cov_PE.xls',header=FALSE)
colnames(hichew_bin_cov)<-c('cell','bin','non_zero_cov_pct')
hichew_bin_cov$sample<-'snHiChew_gm12878'

dip_C_bin_cov<-read.table('gageseq_gm12878_bin_cov_PE.xls',header=FALSE)
colnames(dip_C_bin_cov)<-c('cell','bin','non_zero_cov_pct')
dip_C_bin_cov$sample<-'gageseq_gm12878'

min_perturbation <- -0.1
max_perturbation <- 0.1
dip_C_bin_cov$non_zero_cov_pct <- dip_C_bin_cov$non_zero_cov_pct + runif(n = nrow(dip_C_bin_cov), min = min_perturbation, max = max_perturbation)

sciHiC_bin_cov<-read.table('sciHiC_gm12878_bin_cov_PE.xls',header=FALSE)
colnames(sciHiC_bin_cov)<-c('cell','bin','non_zero_cov_pct')
sciHiC_bin_cov$sample<-'sciHiC_gm12878'

stack_df<-rbind.data.frame(sciHiC_bin_cov,dip_C_bin_cov,hichew_bin_cov)
stack_df$sample<-factor(stack_df$sample,levels=c('sciHiC_gm12878','gageseq_gm12878','snHiChew_gm12878'))
stack_df<-stack_df[stack_df$bin!=100000,]
stack_df$bin<-as.factor(stack_df$bin)

pdf('top200_bin_cov_PE.pdf',width = 10,height = 5)
p2<-ggplot(stack_df, aes(x = sample, y = non_zero_cov_pct,group=sample,color=sample))+ geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=0.4, alpha=0.9)+facet_wrap(~bin,nrow =1)+scale_color_manual(values=c('#3C5488B2','#E64B35B2','#4DBBD5B2'))+
  theme_light()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                       strip.background = element_rect(color="black", fill="white"),
                       strip.text.x = element_text(size = 12, color = "black"),
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(label=c('enhancer promoter overlapping bins, cell of top200 reads'))
p2
dev.off()

svg('top200_bin_cov_PE.svg',width = 10,height = 5)
p2
dev.off()












############ stack plot ############ 

nagano2017<-read.table('nagano2017_downsample_rmdup_validpair.xls',header=FALSE)
nagano2017$sample<-'nagano2017'
colnames(nagano2017)<-c('readNum_per_cell','mean_rmdup_valid_pair','sd','sample')

dip_C<-read.table('dip_C_downsample_rmdup_validpair.xls',header=FALSE)
dip_C$sample<-'dip_C'
colnames(dip_C)<-c('readNum_per_cell','mean_rmdup_valid_pair','sd','sample')

#scCARE<-read.table('scCARE_seq_downsample_rmdup_validpair.xls',header=FALSE)
#scCARE$sample<-'scCARE'
#colnames(scCARE)<-c('readNum_per_cell','mean_rmdup_valid_pair','sd','sample')

hichew<-read.table('hichew_20_lanes_downsample_rmdvalidpair.xls',header=FALSE)
hichew$sample<-'hichew'
colnames(hichew)<-c('readNum_per_cell','mean_rmdup_valid_pair','sd','sample')

stack_df<-rbind.data.frame(nagano2017,dip_C,hichew)
#stack_df$mean_rmdup_valid_pair<-log2(stack_df$mean_rmdup_valid_pair)
#stack_df$sd<-log10(stack_df$sd)
1321651/214280

pdf('rmdvalidpair_saturation_curve.pdf',height=5,width=8)
ggplot(stack_df, aes(x = readNum_per_cell, y = mean_rmdup_valid_pair,group=sample,color=sample))+ 
  geom_line()+ geom_point()+
  geom_pointrange( aes(x=readNum_per_cell, ymin=mean_rmdup_valid_pair-sd, ymax=mean_rmdup_valid_pair+sd),
                     alpha=0.3, size=0.3)+
  #scale_x_continuous(breaks = c(1e5,1e6,1e7))+
  theme_classic()+scale_color_manual(values=c(mypal[1],mypal[2],mypal[4]))+
  labs(subtitle='gradient: 25000 50000 100000 200000 400000 
       800000 1600000 3200000 6400000 12800000')+xlim(0,15000000)
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




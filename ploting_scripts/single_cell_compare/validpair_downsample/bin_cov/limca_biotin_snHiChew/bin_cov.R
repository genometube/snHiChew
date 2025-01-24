library(ggplot2)
library(purrr)
library(dplyr)
library(ggsci)
library(reshape2)

mypal <- pal_npg("nrc", alpha = 0.7)(9)

hichew_bin_cov<-read.table('snHiChew_gm12878_4M_read_bin_cov.xls',header=FALSE)
hichew_bin_cov<-hichew_bin_cov[,-3]
colnames(hichew_bin_cov)<-c('cell','bin','non_zero_cov_pct')
hichew_bin_cov$sample<-'snHiChew_gm12878'

#dip_C_bin_cov<-read.table('LiMCA_4M_read_bin_cov.xls',header=FALSE)
#dip_C_bin_cov<-dip_C_bin_cov[,-3]
#colnames(dip_C_bin_cov)<-c('cell','bin','non_zero_cov_pct')
#dip_C_bin_cov$sample<-'LiMCA'

dip_C_bin_cov<-read.table('zhichao_dipC_gm12878_4M_read_bin_cov.xls',header=FALSE)
dip_C_bin_cov<-dip_C_bin_cov[,-3]
colnames(dip_C_bin_cov)<-c('cell','bin','non_zero_cov_pct')
dip_C_bin_cov$sample<-'zhichao_dipC_gm12878'

#min_perturbation <- -0.1
#max_perturbation <- 0.1
#dip_C_bin_cov$non_zero_cov_pct <- dip_C_bin_cov$non_zero_cov_pct + runif(n = nrow(dip_C_bin_cov), min = min_perturbation, max = max_perturbation)

sciHiC_bin_cov<-read.table('sciHiC_gm12878_bin_cov.xls',header=FALSE)
colnames(sciHiC_bin_cov)<-c('cell','bin','non_zero_cov_pct')
sciHiC_bin_cov$sample<-'sciHiC_gm12878'

stack_df<-rbind.data.frame(sciHiC_bin_cov,dip_C_bin_cov,hichew_bin_cov)
stack_df$sample<-factor(stack_df$sample,levels=c('sciHiC_gm12878','zhichao_dipC_gm12878','snHiChew_gm12878'))
stack_df<-stack_df[stack_df$bin!=100000,]
stack_df$bin<-as.factor(stack_df$bin)

stack_df_validpair_mean<-stack_df %>% group_by(sample,bin) %>% summarize(mean_pct=mean(non_zero_cov_pct))
ggplot(stack_df_validpair_mean,aes(x=sample,y=mean_pct))+geom_point()+facet_wrap(~bin,nrow=1)

pdf('bin_cov.pdf',width = 10,height = 5)
p1<-ggplot(stack_df, aes(x = sample, y = non_zero_cov_pct,group=sample,color=sample))+ geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=0.4, alpha=0.9)+facet_wrap(~bin,nrow =1)+scale_color_manual(values=c('#3C5488B2','#E64B35B2','#4DBBD5B2'))+
  theme_light()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                       strip.background = element_rect(color="black", fill="white"),
                       strip.text.x = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(label=c('sciHiC: cell of top200 reads, \nsnHiChew: 4M reads, 33 cells \n zhichao dipC: 4M reads, 215 cells'))
p1
dev.off()

svg('bin_cov.svg',width = 10,height = 5)
p1
dev.off()



hichew_bin_cov<-read.table('snHiChew_gm12878_4M_read_bin_cov_PE.xls',header=FALSE)
hichew_bin_cov<-hichew_bin_cov[,-3]
colnames(hichew_bin_cov)<-c('cell','bin','non_zero_cov_pct')
hichew_bin_cov$sample<-'snHiChew_gm12878'

#dip_C_bin_cov<-read.table('LiMCA_4M_read_bin_cov_PE.xls',header=FALSE)
#dip_C_bin_cov<-dip_C_bin_cov[,-3]
#colnames(dip_C_bin_cov)<-c('cell','bin','non_zero_cov_pct')
#dip_C_bin_cov$sample<-'LiMCA'

dip_C_bin_cov<-read.table('zhichao_dipC_gm12878_4M_read_bin_cov_PE.xls',header=FALSE)
dip_C_bin_cov<-dip_C_bin_cov[,-3]
colnames(dip_C_bin_cov)<-c('cell','bin','non_zero_cov_pct')
dip_C_bin_cov$sample<-'zhichao_dipC_gm12878'


#min_perturbation <- -0.1
#max_perturbation <- 0.1
#dip_C_bin_cov$non_zero_cov_pct <- dip_C_bin_cov$non_zero_cov_pct + runif(n = nrow(dip_C_bin_cov), min = min_perturbation, max = max_perturbation)

sciHiC_bin_cov<-read.table('sciHiC_gm12878_bin_cov_PE.xls',header=FALSE)
colnames(sciHiC_bin_cov)<-c('cell','bin','non_zero_cov_pct')
sciHiC_bin_cov$sample<-'sciHiC_gm12878'

stack_df<-rbind.data.frame(sciHiC_bin_cov,dip_C_bin_cov,hichew_bin_cov)
stack_df$sample<-factor(stack_df$sample,levels=c('sciHiC_gm12878','zhichao_dipC_gm12878','snHiChew_gm12878'))
stack_df<-stack_df[stack_df$bin!=100000,]
stack_df$bin<-as.factor(stack_df$bin)


stack_df_validpair_mean<-stack_df %>% group_by(sample,bin) %>% summarize(mean_pct=mean(non_zero_cov_pct))
ggplot(stack_df_validpair_mean,aes(x=sample,y=mean_pct))+geom_point()+facet_wrap(~bin,nrow=1)


pdf('bin_cov_PE.pdf',width = 10,height = 5)
p2<-ggplot(stack_df, aes(x = sample, y = non_zero_cov_pct,group=sample,color=sample))+ geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=0.4, alpha=0.9)+facet_wrap(~bin,nrow =1)+scale_color_manual(values=c('#3C5488B2','#E64B35B2','#4DBBD5B2'))+
  theme_light()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                       strip.background = element_rect(color="black", fill="white"),
                       strip.text.x = element_text(size = 12, color = "black"),
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(label=c('Promoter enhancer bin\nsciHiC: cell of top200 reads, \nsnHiChew: 4M reads, 33 cells \n LiMCA: 4M reads, 215 cells'))
p2
dev.off()

svg('bin_cov_PE.svg',width = 10,height = 5)
p2
dev.off()











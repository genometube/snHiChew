library(ggplot2)
library(purrr)
library(dplyr)
library(ggsci)
mypal <- pal_npg("nrc", alpha = 0.7)(9)

hichew_bin_cov<-read.table('hichew_bin_cov_PE.xls',header=FALSE)
colnames(hichew_bin_cov)<-c('cell','bin','non_zero_cov_pct')
hichew_bin_cov$sample<-'hichew'

dip_C_bin_cov<-read.table('dip_C_bin_cov_PE.xls',header=FALSE)
colnames(dip_C_bin_cov)<-c('cell','bin','non_zero_cov_pct')
dip_C_bin_cov$sample<-'dip_C'

nagano_2017_bin_cov<-read.table('nagano2017_bin_cov_PE.xls',header=FALSE)
colnames(nagano_2017_bin_cov)<-c('cell','bin','non_zero_cov_pct')
nagano_2017_bin_cov$sample<-'nagano_2017'

stack_df<-rbind.data.frame(nagano_2017_bin_cov,dip_C_bin_cov,hichew_bin_cov)
stack_df$sample<-factor(stack_df$sample,levels=c('nagano_2017','dip_C','hichew'))
stack_df$bin<-as.factor(stack_df$bin)

pdf('top200_bin_cov_PE.pdf',width = 8,height = 5)
ggplot(stack_df, aes(x = sample, y = non_zero_cov_pct,group=sample,color=sample))+ geom_boxplot()+
  geom_jitter(size=0.4, alpha=0.9)+facet_wrap(~bin,nrow =1)+scale_color_manual(values=c('#3C5488B2','#E64B35B2','#4DBBD5B2'))+
  theme_light()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                       strip.background = element_rect(color="black", fill="white"),
                       strip.text.x = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(label=c('Promoter/Enhancer bin, cell of top200 reads'))

dev.off()



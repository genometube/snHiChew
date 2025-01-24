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

################## mito noise dist ##################
dipC<-read.table('dipC_out_mito_false.xls',header=TRUE)
cycle15<-read.table('cycle_15_out_mito_false.xls',header=TRUE)
cycle18<-read.table('cycle_18_out_mito_false.xls',header=TRUE)
cycle20<-read.table('cycle_20_out_mito_false.xls',header=TRUE)

stack<-rbind.data.frame(cbind.data.frame(sample='dipC',NM_MM_ratio=dipC$NM_MM_ratio),
                        cbind.data.frame(sample='cycle15',NM_MM_ratio=cycle15$NM_MM_ratio),
                        cbind.data.frame(sample='cycle18',NM_MM_ratio=cycle18$NM_MM_ratio),
                        cbind.data.frame(sample='cycle20',NM_MM_ratio=cycle20$NM_MM_ratio))

stack<-rbind.data.frame(cbind.data.frame(sample='dipC',NM_MM_ratio=dipC$NM_MM_ratio),
                        cbind.data.frame(sample='cycle15',NM_MM_ratio=cycle15$NM_MM_ratio))
#stack<-stack[stack$chrM_chr_count>1,]
pdf('mito_noise_dist.pdf',width=8,height=8)
ggplot(stack, aes(x=NM_MM_ratio,color=sample))+ geom_density(lwd=1)+ theme_classic()+scale_color_manual(values=c(mypal[2],mypal[1],mypal[4],mypal[9]))+
  labs(title='chrM vs chr1-22 X Y pair normalized by mito vs mito pair num,
       cycle15 N=356 , cycle18 N=186, cycle20 N=188, dipC=29')
dev.off()


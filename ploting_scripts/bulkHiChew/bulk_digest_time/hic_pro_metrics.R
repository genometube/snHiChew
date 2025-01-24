library(ggplot2)
library(gridExtra)
library(ggsci)
library(dplyr)

mypal <- pal_npg("nrc", alpha = 0.7)(9)

s1<-read.table('L01_hic_pro_metrics.xls',header=TRUE)
s1<-s1[c(-3:-5),]
s1$sample<-c('bulkHiChew_293T_2h_5000','bulkHiChew_293T_4h_5000')

s2<-read.table('L02_hic_pro_metrics.xls',header=TRUE)
s2<-s2[c(-3:-5),]
s2$sample<-c('bulkHiChew_293T_2h_5000','bulkHiChew_293T_4h_5000')

#s3<-read.table('L03_hic_pro_metrics.xls',header=TRUE)
#s3<-s3[c(-3:-5),]
#s3$sample<-c('bulkHiChew_293T_2h_5000','bulkHiChew_293T_4h_5000')

s4<-read.table('bulkHiChew_5k_500_50_hic_pro_metrics.xls',header=TRUE)
s4<-s4[c(1:2),]
s4$sample<-c('bulkHiChew_293T_overnight_5000')


color_palette_raw<-c(mypal[1],mypal[2],mypal[4],mypal[9],mypal[5],mypal[3],mypal[6])
color_palette<-color_palette_raw[1:3]

group_df<-rbind.data.frame(s1,s2,s4)
valid_report_ratio<-cbind.data.frame(V1=group_df$valid_pairs/group_df$Reported_pairs,cell_lib=group_df$sample)
summary_df <- valid_report_ratio %>%
  group_by(cell_lib) %>%
  summarise(mean_value = mean(V1),
            se = sd(V1) / sqrt(n()))
p2 <- ggplot(valid_report_ratio, aes(x = cell_lib, y = V1,color=cell_lib)) +
  geom_point(position = position_jitter(width = 0.1), size = 3) +  # Add dots with jitter
  geom_errorbar(data = summary_df, aes(y = mean_value, ymin = mean_value - se, ymax = mean_value + se), width = 0.2) +ylim(0,1)+  # Add error bars
  labs(title="valid / report ratio")+
  theme_classic()+scale_fill_manual(values=color_palette)+theme(axis.text.x = element_text(angle = 25, hjust = 1))+scale_color_manual(values=color_palette)
p2

#p2<-ggplot(data = valid_report_ratio, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_jitter(width=0.3,alpha=0.5) +
#  geom_boxplot(width=0.1,outlier.shape = NA)+labs(title="valid / report ratio")+
#  theme_classic()+scale_fill_manual(values=color_palette)+theme(legend.position = "top") +
#  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
#  theme(axis.text.x = element_text(angle = 15, hjust = 1))
pdf(paste0('digest_time_5k_cell.pdf'),height=8,width=6)
p2
dev.off()
svg(paste0('digest_time_5k_cell.svg'),height=8,width=6)
p2
dev.off()



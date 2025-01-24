library(ggplot2)
library(gridExtra)
library(ggsci)
library(dplyr)

mypal <- pal_npg("nrc", alpha = 0.7)(9)
color_palette<-c(mypal[1],mypal[2],mypal[3],mypal[4],mypal[5],mypal[6],mypal[7],mypal[8],mypal[9],'darkred','darkgreen','darkblue')

group_df<-read.table('bulkHiChew_5k_500_50_hic_pro_metrics.xls',header=TRUE)
remove_last_underscore <- function(x) {sub("_[^_]*$", "", x)}
group_df$sample_ <- sapply(group_df$sample, remove_last_underscore)

valid_report_ratio<-cbind.data.frame(V1=group_df$valid_pairs/group_df$Reported_pairs,cell_lib=group_df$sample_)
summary_df <- valid_report_ratio %>% group_by(cell_lib) %>%  summarise(mean_value = mean(V1), se = sd(V1) / sqrt(n()))

# Create the dot plot
p2 <- ggplot(valid_report_ratio, aes(x = cell_lib, y = V1,color=cell_lib)) +
  geom_point(position = position_jitter(width = 0.1), size = 3) +  # Add dots with jitter
  geom_errorbar(data = summary_df, aes(y = mean_value, ymin = mean_value - se, ymax = mean_value + se), width = 0.2) +ylim(0,1)+  # Add error bars
  labs(title="valid / report ratio")+
  theme_classic()+scale_fill_manual(values=color_palette)+theme(axis.text.x = element_text(angle = 25, hjust = 1))+scale_color_manual(values=color_palette)
p2

#p2<-ggplot(data = valid_report_ratio, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_jitter(width=0.3,alpha=0.5) +
#  geom_boxplot(width=0.05,outlier.shape = NA)+labs(title="valid / report ratio")+
#  theme_classic()+scale_fill_manual(values=color_palette)+theme(axis.text.x = element_text(angle = 25, hjust = 1))+scale_color_manual(values=color_palette)

# Calculate percentages
group_df$gatc_non <- (group_df$gatc_non / group_df$Reported_pairs) 
group_df$gatc_mid <- (group_df$gatc_mid / group_df$Reported_pairs) 
group_df$gatc_side <- (group_df$gatc_side / group_df$Reported_pairs) 
library(reshape2)
group_df_<-group_df[c(1:6),]
raw_df_melted <- melt(group_df_, id.vars = 'sample', measure.vars = c('gatc_non', 'gatc_mid', 'gatc_side'))

p1<-ggplot(raw_df_melted, aes(x = sample, y = value, fill = variable)) +
  geom_bar(stat = 'identity',alpha=0.8) +
  labs(x = 'Samples', y = 'Percentage', title = 'gatc_side: gatc at start or end 5bp, \n3 types of dangling end pairs / Reported_pairs') +
  theme_classic()+scale_fill_manual(values=color_palette)+theme(axis.text.x = element_text(angle = 25, hjust = 1))+ylim(0,0.7)+scale_color_manual(values=color_palette)
p1
dev.off()

cis_trans_ratio<-cbind.data.frame(V1=log2((group_df$cis_shortRange+group_df$cis_longRange)/group_df$trans_interaction),cell_lib=group_df$sample_)
summary_df <- cis_trans_ratio %>% group_by(cell_lib) %>%  summarise(mean_value = mean(V1), se = sd(V1) / sqrt(n()))

# Create the dot plot
p5 <- ggplot(cis_trans_ratio, aes(x = cell_lib, y = V1,color=cell_lib)) +
  geom_point(position = position_jitter(width = 0.1), size = 3) +  # Add dots with jitter
  geom_errorbar(data = summary_df, aes(y = mean_value, ymin = mean_value - se, ymax = mean_value + se), width = 0.2) +ylim(-0,4)+  # Add error bars
  labs(title="cis trans ratio")+
  theme_classic()+scale_fill_manual(values=color_palette)+theme(axis.text.x = element_text(angle = 25, hjust = 1))+scale_color_manual(values=color_palette)
p5

#p5<-ggplot(data = cis_trans_ratio, aes(y = V1,x = cell_lib,fill=cell_lib)) + geom_jitter(width=0.2,alpha=0.5) +
#  geom_boxplot(width=0.05,outlier.shape = NA)+labs(title="cis_trans_ratio")+ylim(-0,4)+
#  theme_classic()+scale_fill_manual(values=color_palette)+theme(axis.text.x = element_text(angle = 25, hjust = 1))+scale_color_manual(values=color_palette)
#p5
dev.off()

pdf(paste0('hichew_pcr_50_500_5k.pdf'),height=6,width=15)
plotlist<-list(p1,p2,p5)
cowplot::plot_grid(plotlist = plotlist, ncol = 3)
dev.off()

svg(paste0('hichew_pcr_50_500_5k.svg'),height=6,width=15)
plotlist<-list(p1,p2,p5)
cowplot::plot_grid(plotlist = plotlist, ncol = 3)
dev.off()


#############################
library(ggplot2)
library(dplyr)
# Sample data frame
df <- data.frame(
  group = rep(c("A", "B"), each = 2),
  value = c(5, 7, 3, 4)
)

# Calculate mean and standard error for each group
summary_df <- df %>%
  group_by(group) %>%
  summarise(mean_value = mean(value),
    se = sd(value) / sqrt(n()))

# Create the dot plot
p <- ggplot(df, aes(x = group, y = value)) +
  geom_point(position = position_jitter(width = 0.1), size = 3) +  # Add dots with jitter
  geom_errorbar(data = summary_df, aes(y = mean_value, ymin = mean_value - se, ymax = mean_value + se), width = 0.2) +  # Add error bars
  geom_point(data = summary_df, aes(y = mean_value), color = "red", size = 4) +  # Add mean points
  theme_minimal() +
  labs(title = "Dot Plot with Mean and Error Bars", x = "Group", y = "Value")
p

library(ggplot2)
library(plyr)
library(reshape2)
library(zoo)
library(Matrix)
library(data.table)
library(dplyr)
library(ggsci)
mypal <- pal_npg("nrc", alpha = 0.7)(9)

### RNAseq CTCF expression ###
wt_kd_RNA<-read.table('groupA-vs-groupC.gene.DEseq2.Filter.xls',header = T)
wt_kd_RNA_subset<-wt_kd_RNA[wt_kd_RNA$GeneID=='10664',c(20,21,22,36,37,38)]
df<-reshape2::melt(wt_kd_RNA_subset)

df$condition<-sub(".*_(.).*", "\\1", df$variable)
colnames(df)<-c('sample','expression','condition')
df$condition<-c('wt','wt','wt','ctcfkd','ctcfkd','ctcfkd')
df$condition<-factor(df$condition, levels = c('wt','ctcfkd'))

summary_data <- df %>%  group_by(condition) %>%  summarise(    mean_expression = mean(expression),
    se_expression = sd(expression) / sqrt(n())  )

t_test_result <- t.test(expression ~ condition, data = df)

summary_df <- df %>%
  group_by(condition) %>%
  summarise(mean_value = mean(expression),
            se = sd(expression) / sqrt(n()))
mypal <- pal_npg("nrc", alpha = 0.7)(9)

p1 <- ggplot(df, aes(x = condition, y = expression,color=condition)) +
  geom_point(position = position_jitter(width = 0.1), size = 2,alpha=0.7) +  # Add dots with jitter
  geom_errorbar(data = summary_df, aes(y = mean_value, ymin = mean_value - se, ymax = mean_value + se), width = 0.2) +  # Add error bars
  labs(title=paste('RNAseq CTCF expression\n',"p-value:", format(t_test_result$p.value, digits = 3)))+
  theme_classic()+theme(axis.text.x = element_text(angle = 0, hjust = 1))+ scale_color_manual(values=c('#999999','#0ABAB5'))
p1

### qPCR CTCF expression ###
wt_kd_subset<-cbind.data.frame(wt_1=0.973,wt_2=0.863,wt_3=1.192,ctcfkd_1=0.389,ctcfkd_2=0.574,ctcfkd_3=0.741)
df<-reshape2::melt(wt_kd_subset)

df$condition <- sub("^(.*?)_.*", "\\1", df$variable)
colnames(df)<-c('sample','expression','condition')
df$condition<-factor(df$condition, levels = c('wt','ctcfkd'))

summary_data <- df %>%  group_by(condition) %>%  summarise(mean_expression = mean(expression), se_expression = sd(expression) / sqrt(n())  )
t_test_result <- t.test(expression ~ condition, data = df)

summary_df <- df %>%
  group_by(condition) %>%
  summarise(mean_value = mean(expression),
            se = sd(expression) / sqrt(n()))

p2 <- ggplot(df, aes(x = condition, y = expression,color=condition)) +
  geom_point(position = position_jitter(width = 0.1), size = 2,alpha=0.7) +  # Add dots with jitter
  geom_errorbar(data = summary_df, aes(y = mean_value, ymin = mean_value - se, ymax = mean_value + se), width = 0.2) +  # Add error bars
  labs(title=paste('qPCR CTCF expression\n',"p-value:", format(t_test_result$p.value, digits = 3)))+
  theme_classic()+theme(axis.text.x = element_text(angle = 0, hjust = 1))+ scale_color_manual(values=c('#999999','#0ABAB5'))
p2

### sfig 5e ab efficiency ###
stack<-read.table('bulkHiChew_5k_500_50_hic_pro_metrics.xls',header=T)
stack$gatc_mid_ratio<-round(stack$gatc_mid/stack$Dangling_end_pairs,2)
df<-stack[c(1,2,7,8),c(1,21)]

df$condition <- sub("^(.*?)_.*", "\\1", df$sample)

colnames(df)<-c('sample','gatc_mid_ratio','condition')
df$condition<-factor(df$condition, levels = c('bulkHiChew','bulkPCR'))

summary_data <- df %>%  group_by(condition) %>%  summarise(mean_expression = mean(gatc_mid_ratio), se_expression = sd(gatc_mid_ratio) / sqrt(n())  )
t_test_result <- t.test(gatc_mid_ratio ~ condition, data = df)

summary_df <- df %>%  group_by(condition) %>% summarise(mean_value = mean(gatc_mid_ratio), se = sd(gatc_mid_ratio) / sqrt(n()))

p3 <- ggplot(df, aes(x = condition, y = gatc_mid_ratio,color=condition)) +
  geom_point(position = position_jitter(width = 0.1), size = 2,alpha=0.7) +  # Add dots with jitter
  geom_errorbar(data = summary_df, aes(y = mean_value, ymin = mean_value - se, ymax = mean_value + se), width = 0.2) +  # Add error bars
  labs(title=paste('sfig 5e ab efficiency\n',"p-value:", format(t_test_result$p.value, digits = 3)))+
  theme_classic()+theme(axis.text.x = element_text(angle = 0, hjust = 1))+ scale_color_manual(values=c('#999999','#0ABAB5'))
p3


### sfig 5d ab pulldown ###
wt_kd_subset<-cbind.data.frame(cell5k_oneIP_1=2.51,cell5k_oneIP_2=2.6,cell500_oneIP_1=2.82,cell500_oneIP_2=2.7,cell50_oneIP_1=3.18,cell50_oneIP_2=3.14,
                               cell5k_twoIP_1=2.77,cell5k_twoIP_2=2.85,cell500_twoIP_1=3.28,cell500_twoIP_2=3.13,cell50_twoIP_1=3.46,cell50_twoIP_2=3.43,
                               cell5k_sdsIP_1=2.53,cell5k_sdsIP_2=2.62,cell500_sdsIP_1=2.84,cell500_sdsIP_2=2.72,cell50_sdsIP_1=3.21,cell50_sdsIP_2=3.16)

df<-reshape2::melt(wt_kd_subset)
df$cond1 <- sub("^(.*?)_.*", "\\1", df$variable)
df$cond2 <- sub("^[^_]*_(.*?)_.*", "\\1", df$variable)
df$cond3 <- sub("^([^_]*_[^_]*)_.*", "\\1", df$variable)
colnames(df) <- c('sample', 'value', 'cond1', 'cond2', 'cond3')
summary_df <- df %>%  group_by(cond3) %>%  summarise(mean_value = mean(value), se = sd(value) / sqrt(n()))

p4 <- ggplot(df, aes(x = cond3, y = value,color=cond3)) +
  geom_point(position = position_jitter(width = 0.1), size = 2,alpha=0.7) +  # Add dots with jitter
  geom_errorbar(data = summary_df, aes(y = mean_value, ymin = mean_value - se, ymax = mean_value + se), width = 0.2) +  # Add error bars
  labs(title=paste('sfig 5d ab pulldown\n'))+
  theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ ylim(0,4)+
  scale_color_manual(values=c('#999999','#0ABAB5','coral3','#999999','#0ABAB5','coral3','#999999','#0ABAB5','coral3'))
p4

### frag size ###
wt_kd_subset<-cbind.data.frame(frag550_1=0.403,frag550_2=0.421,frag550_3=0.411,
                               frag450_1=0.371,frag450_2=0.377,frag450_3=0.387,
                               frag350_1=0.310,frag350_2=0.314,frag350_3=0.308)
df<-reshape2::melt(wt_kd_subset)

df$condition <- sub("^(.*?)_.*", "\\1", df$variable)
colnames(df)<-c('sample','value','condition')
#df$condition<-factor(df$condition, levels = c('wt',''))

summary_df <- df %>%  group_by(condition) %>%  summarise(mean_value = mean(value), se = sd(value) / sqrt(n())  )
#t_test_result <- t.test(expression ~ condition, data = df)

p5 <- ggplot(df, aes(x = condition, y = value,color=condition)) +
  geom_point(position = position_jitter(width = 0.1), size = 2,alpha=0.7) +  # Add dots with jitter
  geom_errorbar(data = summary_df, aes(y = mean_value, ymin = mean_value - se, ymax = mean_value + se), width = 0.2) +  # Add error bars
  labs(title=paste('frag size valid pair ratio\n'))+ ylim(0,1)+
  theme_classic()+theme(axis.text.x = element_text(angle = 0, hjust = 1))+ scale_color_manual(values=c('#999999','#0ABAB5','coral3'))
p5


### frag size ###
wt_kd_subset<-cbind.data.frame(frag550_1=0.403,frag550_2=0.421,frag550_3=0.411,
                               frag450_1=0.371,frag450_2=0.377,frag450_3=0.387,
                               frag350_1=0.310,frag350_2=0.314,frag350_3=0.308)
df<-reshape2::melt(wt_kd_subset)

df$condition <- sub("^(.*?)_.*", "\\1", df$variable)
colnames(df)<-c('sample','value','condition')
#df$condition<-factor(df$condition, levels = c('wt',''))

summary_df <- df %>%  group_by(condition) %>%  summarise(mean_value = mean(value), se = sd(value) / sqrt(n())  )
#t_test_result <- t.test(expression ~ condition, data = df)

p5 <- ggplot(df, aes(x = condition, y = value,color=condition)) +
  geom_point(position = position_jitter(width = 0.1), size = 2,alpha=0.7) +  # Add dots with jitter
  geom_errorbar(data = summary_df, aes(y = mean_value, ymin = mean_value - se, ymax = mean_value + se), width = 0.2) +  # Add error bars
  labs(title=paste('sfig 5g frag size valid pair ratio\n'))+ ylim(0,1)+
  theme_classic()+theme(axis.text.x = element_text(angle = 0, hjust = 1))+ scale_color_manual(values=c('#999999','#0ABAB5','coral3'))
p5


### frag size ###
wt_kd_subset<-cbind.data.frame(millipore_1=0.453,millipore_2=0.421,millipore_3=0.481,
                               abcam_1=0.321,abcam_2=0.357,abcam_3=0.306)
df<-reshape2::melt(wt_kd_subset)

df$condition <- sub("^(.*?)_.*", "\\1", df$variable)
colnames(df)<-c('sample','value','condition')
#df$condition<-factor(df$condition, levels = c('wt',''))

summary_df <- df %>%  group_by(condition) %>%  summarise(mean_value = mean(value), se = sd(value) / sqrt(n())  )
t_test_result <- t.test(value ~ condition, data = df)
p6 <- ggplot(df, aes(x = condition, y = value,color=condition)) +
  geom_point(position = position_jitter(width = 0.1), size = 2,alpha=0.7) +  # Add dots with jitter
  geom_errorbar(data = summary_df, aes(y = mean_value, ymin = mean_value - se, ymax = mean_value + se), width = 0.2) +  # Add error bars
  labs(title=paste('sfig 5c antibody valid pair ratio\n',"p-value:", format(t_test_result$p.value, digits = 3)))+ ylim(0,0.7)+
  theme_classic()+theme(axis.text.x = element_text(angle = 0, hjust = 1))+ scale_color_manual(values=c('#999999','#0ABAB5','coral3'))
p6

### lib compare ###
wt_kd_subset<-cbind.data.frame(mdaInvitro_1=0.32,mdaInvitro_2=0.28,mdaInvitro_3=0.35,
                               tn5Invitro_1=0.34,tn5Invitro_2=0.27,tn5Invitro_3=0.26,
                               tn5Invivo_1=0.678,tn5Invivo_2=0.68,tn5Invivo_3=0.71)
df<-reshape2::melt(wt_kd_subset)

df$condition <- sub("^(.*?)_.*", "\\1", df$variable)
colnames(df)<-c('sample','value','condition')
#df$condition<-factor(df$condition, levels = c('wt',''))

summary_df <- df %>%  group_by(condition) %>%  summarise(mean_value = mean(value), se = sd(value) / sqrt(n())  )
#t_test_result <- t.test(value ~ condition, data = df)
p7 <- ggplot(df, aes(x = condition, y = value,color=condition)) +
  geom_point(position = position_jitter(width = 0.1), size = 2,alpha=0.7) +  # Add dots with jitter
  geom_errorbar(data = summary_df, aes(y = mean_value, ymin = mean_value - se, ymax = mean_value + se), width = 0.2) +  # Add error bars
  labs(title=paste('sfig 5j lib compare valid pair ratio\n',"p-value:", format(t_test_result$p.value, digits = 3)))+ ylim(0,0.7)+
  theme_classic()+theme(axis.text.x = element_text(angle = 0, hjust = 1))+ scale_color_manual(values=c('#999999','#0ABAB5','coral3'))
p7

library(gridExtra)
pdf('suppl_dotplots_250113.pdf',width = 16,height=10)
grid.arrange(p1,p2,p3,p5,p6,p7,p4,layout_matrix=rbind(c(1,2,3,4),c(5,6,7,7)))
dev.off()

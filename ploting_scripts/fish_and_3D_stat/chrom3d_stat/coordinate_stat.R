library(ggplot2)
library(purrr)
library(dplyr)
library(ggsci)
mypal <- pal_npg("nrc", alpha = 0.7)(9)

######################################## 20241021
library(data.table)

coordinate_stat_raw<-fread('coordinate_stat_lamin.xls.gz',header=F)
coordinate_stat_lamin<-coordinate_stat_raw[,c(1:9)]
colnames(coordinate_stat_lamin)<-c('chr','start','end','x','y','z','cell','type','lamin')
coordinate_stat_lamin<-coordinate_stat_lamin[coordinate_stat_lamin$lamin!=".",]

coordinate_stat<-coordinate_stat_raw[,c(1:8)]
colnames(coordinate_stat)<-c('chr','start','end','x','y','z','cell','type')
coordinate_stat_subset<-coordinate_stat[coordinate_stat$chr=="chr7" & coordinate_stat$start>145000000 & coordinate_stat$end<148000000]

# Define a function to calculate the Euclidean distance
euclidean_distance <- function(x1, y1, z1, x2, y2, z2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)}

# Calculate the closest distance for each row in table1
coordinate_stat_subset[, closest_distance := sapply(1:.N, function(i) {
  distances <- euclidean_distance(x[i], y[i], z[i], coordinate_stat_lamin$x, coordinate_stat_lamin$y, coordinate_stat_lamin$z)
  min(distances)
})]

coordinate_stat_subset<-coordinate_stat_subset[coordinate_stat_subset$closest_distance>0,]
coordinate_stat_subset_df<-as.data.frame(coordinate_stat_subset)
wilcox.test(coordinate_stat_subset_df[coordinate_stat_subset_df$type=="melting",9],coordinate_stat_subset_df[coordinate_stat_subset_df$type=="concretion",9])

mean(coordinate_stat_subset_df[coordinate_stat_subset_df$type=="melting",9])
mean(coordinate_stat_subset_df[coordinate_stat_subset_df$type=="concretion",9])

p1<-ggplot(coordinate_stat_subset_df, aes(x = cell, y = closest_distance,color=type))+ 
  geom_point(stat = "identity",width=0.5,alpha=0.5)+ 
  theme_classic()+ theme(axis.text.x = element_text(angle = 10, hjust = 1))+
  geom_boxplot(alpha=0.3,width=0.4)+scale_color_manual(values=c(mypal[3],mypal[1]))+
  labs(title='chr7:145M-148M TAD beads distance to lamin clad',subtitle='3 melting cells and 3 concretion cells, wilcox p-value = 5.781e-07\n
       mean melting: 0.538, mean concretion: 0.274')
p1

p2<-ggplot(coordinate_stat_subset_df, aes(x = type, y = closest_distance,color=type))+ 
  geom_point(stat = "identity",width=0.5,alpha=0.5)+
  theme_classic()+ geom_boxplot(alpha=0.3,width=0.4)+scale_color_manual(values=c(mypal[3],mypal[1]))+
  labs(title='chr7:145M-148M TAD beads distance to lamin clad',subtitle='3 melting cells and 3 concretion cells, wilcox p-value = 5.781e-07\n
       mean melting: 0.538, mean concretion: 0.274')
p2
dev.off()


library(gridExtra)

pdf('tad_bead_distance_to_lamin.pdf',height=16,width=8)
grid.arrange(p1,p2,nrow=2,ncol=1)
dev.off()

svg('tad_bead_distance_to_lamin.svg',height=16,width=8)
grid.arrange(p1,p2,nrow=2,ncol=1)
dev.off()


svg('euclidean_distance_to_nucleolus.svg',height=6,width=6)
p2
dev.off()

pdf('euclidean_distance_to_nucleolus.pdf',height=6,width=6)
p2
dev.off()


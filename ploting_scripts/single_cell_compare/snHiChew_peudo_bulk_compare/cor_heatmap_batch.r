library(ggplot2)
library(plyr)
library(reshape2)
library(zoo)
library(pheatmap)
#library(HiCcompare)
library(Matrix)
library(data.table)
library(ComplexHeatmap)

library(ggplot2)
library(purrr)
library(cowplot)

################################################################################################ a/b compartment track
# Create example data
batch_f<- list(biotin_bulk_HEK293T="hek293_ngsHiC_700M",hek293t200cellb1="snHiChew_HEK293Twt",hek293t200cellb2="snHiChew_200cell293T",
            gm12878_400cellb1="snHiChew_400cellgm12878_20240918",gm12878_400cellb2="snHiChew_400cellgm12878_20240930",bulk_gm12878_4dn="4dn_bulk_gm12878")
batch_data<- map(batch_f,~fread(paste0(.x,'.mcool.res_500000.eig.cis.vecs.tsv')))
eigen_bin=500000

library(tidyverse)
#merge all data frames in list
df_merge<-batch_data %>% reduce(full_join, by=c('chrom','start'))
df_merge[is.na(df_merge)]<-0
df_E1<-df_merge[,c(6,10,14,18,22,26)]
df_E1[is.na(df_E1)]<-0

colnames(df_E1)<-c('biotin_bulk_HEK293T','200cell293Tb1','200cell293Tb2','400cellgm12878b1','400cellgm12878b2','4dn_bulk_gm12878')
corr <- round(cor(df_E1), 3)
typeof(corr)

color_scale =c("royalblue3", "#FDDBC7", "indianred3")
library(ComplexHeatmap)
pdf('snHiChew_cor_eigen_heatmap.pdf',height = 6,width=6)
Heatmap(corr,column_title = paste0('snHiChew eigen correlation'),cluster_rows = F,cluster_columns = F,col=colorRampPalette(color_scale)(150),
        layer_fun = function(j, i, x, y, width, height, fill) {
          # since grid.text can also be vectorized
          grid.text(sprintf("%.2f", pindex(corr, i, j)), x, y, 
                    gp = gpar(fontsize = 10))})
dev.off()

svg('snHiChew_cor_eigen_heatmap.svg',height = 6,width=6)
Heatmap(corr,column_title = paste0('snHiChew eigen correlation'),cluster_rows = F,cluster_columns = F,col=colorRampPalette(color_scale)(150),
        layer_fun = function(j, i, x, y, width, height, fill) {
          # since grid.text can also be vectorized
          grid.text(sprintf("%.2f", pindex(corr, i, j)), x, y, 
                    gp = gpar(fontsize = 10))})
dev.off()


df_E1$x<-as.integer(rownames(df_E1))*eigen_bin
df_E1_long <- reshape2::melt(df_E1, id.vars = "x", variable.name = "sample", value.name = "y")
df_E1_long_sub<-df_E1_long

df_E1_long_sub<-df_E1_long[df_E1_long$x>40000000&df_E1_long$x<70000000,]

colors <- c("darkred", "darkred", "darkred","darkgreen", "darkgreen", "darkgreen")

p2<-ggplot(df_E1_long_sub, aes(x = x, y = y, fill=sample,group = sample)) +
  geom_area(alpha=0.7) + facet_wrap(~ sample, scales = "free_y",ncol=1) +
  labs(title = "eigen score, chr7:1-40M, 500k bin",
       x = "Genomic position", y = "eigen") + theme_classic()+  scale_fill_manual(values=colors)

pdf('snHiChew_cor_eigen_track.pdf',height = 10,width=6)
p2
dev.off()

svg('snHiChew_cor_ins_track.svg',height = 10,width=6)
p2
dev.off()





################################################################################################ ins track
batch_f<- list(biotin_bulk_HEK293T="hek293_ngsHiC_700M",hek293t200cellb1="snHiChew_HEK293Twt",hek293t200cellb2="snHiChew_200cell293T",
               gm12878_400cellb1="snHiChew_400cellgm12878_20240918",gm12878_400cellb2="snHiChew_400cellgm12878_20240930",bulk_4dn_gm12878="4dn_bulk_gm12878")
batch_data<- map(batch_f,~fread(paste0(.x,'.mcool.res_50k.window_100k.ins')))
ins_bin=50000
library(tidyverse)
#merge all data frames in list
df_merge<-batch_data %>% reduce(full_join, by=c('chrom','start'))
df_merge[is.na(df_merge)]<-0
df_merge<-df_merge[complete.cases(df_merge), ]
df_ins<-df_merge[,c(5,10,15,20,25,30)]
colnames(df_ins)<-c('biotin_bulk_HEK293T','200cell293Tb1','200cell293Tb2','400cellgm12878b1','400cellgm12878b2','4dn_bulk_gm12878')

corr <- round(cor(df_ins), 3)

color_scale =c("royalblue3", "#FDDBC7", "indianred3")
library(ComplexHeatmap)
pdf('snHiChew_cor_ins_heatmap.pdf',height = 6,width=6)
Heatmap(corr,column_title = paste0('snHiChew ins correlation'),cluster_rows = F,cluster_columns = F,col=colorRampPalette(color_scale)(150),
        layer_fun = function(j, i, x, y, width, height, fill) {
          # since grid.text can also be vectorized
          grid.text(sprintf("%.2f", pindex(corr, i, j)), x, y, 
                    gp = gpar(fontsize = 10))})
dev.off()

svg('snHiChew_cor_ins_heatmap.svg',height = 6,width=6)
Heatmap(corr,column_title = paste0('snHiChew eigen correlation'),cluster_rows = F,cluster_columns = F,col=colorRampPalette(color_scale)(150),
        layer_fun = function(j, i, x, y, width, height, fill) {
          # since grid.text can also be vectorized
          grid.text(sprintf("%.2f", pindex(corr, i, j)), x, y, 
                    gp = gpar(fontsize = 10))})
dev.off()

df_ins$x<-as.integer(rownames(df_ins))*ins_bin
df_ins_long <- reshape2::melt(df_ins, id.vars = "x", variable.name = "sample", value.name = "y")
df_ins_long_sub<-df_ins_long[df_ins_long$x>41000000&df_ins_long$x<43000000,]

colors <- c("darkred", "darkred", "darkred","darkgreen", "darkgreen", "darkgreen")

p2<-ggplot(df_ins_long_sub, aes(x = x, y = y, fill=sample,group = sample)) +
  geom_area(alpha=0.7) + facet_wrap(~ sample, scales = "free_y",ncol=1) +
  labs(title = "insulation score, chr7:41-43M, 50k bin",
       x = "Genomic position", y = "ins") + theme_classic()+  scale_fill_manual(values=colors)

pdf('snHiChew_cor_ins_track.pdf',height = 10,width=6)
p2
dev.off()

svg('snHiChew_cor_ins_track.svg',height = 10,width=6)
p2
dev.off()




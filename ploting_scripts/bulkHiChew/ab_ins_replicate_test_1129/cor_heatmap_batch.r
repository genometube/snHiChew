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
batch_f<- list(biotin_bulk_hek_5Mcell="hek293_ngsHiC_700M",
               bulkHiChew_5k_hek_b1="bulkHiChew_5000_1",bulkHiChew_5k_hek_b2="bulkHiChew_5000_2",
               bulkHiChew_500_hek_b1="bulkHiChew_500_1",bulkHiChew_500_hek_b2="bulkHiChew_500_2",
               bulkHiChew_50_hek_b1="bulkHiChew_50_1",bulkHiChew_50_hek_b2="bulkHiChew_50_2",
               bulkpcr_5k_hek_b1="bulkPCR_5000_1",
               bulkpcr_500_hek_b1="bulkPCR_500_1",
               bulkpcr_50_hek_b1="bulkPCR_50_1")
batch_data<- map(batch_f,~fread(paste0(.x,'.mcool.res_500000.eig.cis.vecs.tsv')))

eigen_bin=500000

library(tidyverse)
#merge all data frames in list
df_merge<-batch_data %>% reduce(full_join, by=c('chrom','start'))
df_merge[is.na(df_merge)]<-0

subset_col<-as.vector(2+seq(1,length(names(batch_f)))*4)
#df_merge<-df_merge[df_merge$weight.x>0.023,]

df_E1<-df_merge[, .SD, .SDcols = subset_col]
df_E1[is.na(df_E1)]<-0


colnames(df_E1)<-names(batch_f)
corr <- round(cor(df_E1), 2)
typeof(corr)
library(circlize)

#color_scale =c("royalblue3", "grey", "indianred3",'darkred')
col_fun <- colorRamp2(c(0.4, 0.65,0.9), c("royalblue3", "white",'indianred3'))

typeof(corr)
#corr[corr>0.9]<-0.9
breaks <- seq(0.4, 0.9, length.out = 300)

corr[corr]
library(ComplexHeatmap)
pdf('bulkHiChew_cor_eigen_heatmap.pdf',height = 8,width=8)
Heatmap(corr,column_title = paste0('bulkHiChew eigen correlation'),
        cluster_rows = F,cluster_columns = F,col = col_fun,
        layer_fun = function(j, i, x, y, width, height, fill) {
          # since grid.text can also be vectorized
          grid.text(sprintf("%.2f", pindex(corr, i, j)), x, y, 
                    gp = gpar(fontsize = 10))})
dev.off()







batch_f<- list(biotin_bulk_hek_5Mcell="hek293_ngsHiC_700M",
               bulkHiChew_5k_hek_b1="bulkHiChew_5000_1",bulkHiChew_5k_hek_b2="bulkHiChew_5000_2",
               bulkHiChew_500_hek_b1="bulkHiChew_500_1",bulkHiChew_500_hek_b2="bulkHiChew_500_2",
               bulkHiChew_50_hek_b1="bulkHiChew_50_1",bulkHiChew_50_hek_b2="bulkHiChew_50_2",
               bulkpcr_5k_hek_b1="bulkPCR_5000_1",
               bulkpcr_500_hek_b1="bulkPCR_500_1",
               bulkpcr_50_hek_b1="bulkPCR_50_1")
batch_data<- map(batch_f,~fread(paste0(.x,'.mcool.res_500000.eig.cis.vecs.tsv')))

eigen_bin=500000

library(tidyverse)
#merge all data frames in list
df_merge<-batch_data %>% reduce(full_join, by=c('chrom','start'))
df_merge[is.na(df_merge)]<-0

subset_col<-as.vector(2+seq(1,length(names(batch_f)))*4)
df_E1<-df_merge[, .SD, .SDcols = subset_col]
colnames(df_E1)<-names(batch_f)

df_E1$x<-as.integer(rownames(df_E1))*eigen_bin
df_E1_long <- reshape2::melt(df_E1, id.vars = "x", variable.name = "sample", value.name = "y")
df_E1_long_sub<-df_E1_long

df_E1_long_sub<-df_E1_long[df_E1_long$x>40000000&df_E1_long$x<70000000,]

colors <- c("#E41A1C", "#377EB8", "darkgreen", "#4DAF4A", "#FF7F00", "orange", "#984EA3", "#F781BF", "#999999", "#66C2A5", "#FC8D62")

p2<-ggplot(df_E1_long_sub, aes(x = x, y = y,group = sample)) +
  geom_area(alpha=0.5, fill='grey30',color='grey30') + facet_wrap(~ sample, scales = "free_y",ncol=1) +
  labs(title = "eigen score, chr7:40M-70M, 500k bin",
       x = "Genomic position", y = "eigen") + theme_classic()

pdf('bulkHiChew_cor_eigen_track.pdf',height = 10,width=6)
p2
dev.off()

svg('bulkHiChew_cor_eigen_track.svg',height = 10,width=6)
p2
dev.off()


################################################################################################ ins track
batch_f<- list(biotin_bulk_hek_5Mcell="hek293_ngsHiC_700M",
               bulkHiChew_5k_hek_b1="bulkHiChew_5000_1",bulkHiChew_5k_hek_b2="bulkHiChew_5000_2",
               bulkHiChew_500_hek_b1="bulkHiChew_500_1",bulkHiChew_500_hek_b2="bulkHiChew_500_2",
               bulkHiChew_50_hek_b1="bulkHiChew_50_1",bulkHiChew_50_hek_b2="bulkHiChew_50_2",
               bulkpcr_5k_hek_b1="bulkPCR_5000_1",
               bulkpcr_500_hek_b1="bulkPCR_500_1",
               bulkpcr_50_hek_b1="bulkPCR_50_1")
batch_data<- map(batch_f,~fread(paste0(.x,'.mcool.res_50000.window_150k.ins')))
ins_bin=50000
library(tidyverse)

#merge all data frames in list
df_merge<-batch_data %>% reduce(full_join, by=c('chrom','start'))
df_merge[is.na(df_merge)]<-0
df_merge<-df_merge[complete.cases(df_merge), ]
df_merge<-df_merge[df_merge$is_bad_bin.x!="TRUE"&df_merge$is_bad_bin.x.x!="TRUE"
                   &df_merge$log2_insulation_score_150000.x!=0&df_merge$log2_insulation_score_150000.x.x!=0
                   &df_merge$n_valid_pixels_150000.x>5&df_merge$n_valid_pixels_150000.x.x>5
                   &df_merge$n_valid_pixels_150000.y>5&df_merge$n_valid_pixels_150000.y.y>5,]
#ggplot(df_merge, aes(x = df_merge$n_valid_pixels_150000.x)) + geom_density()
       
#colnames(df_ins)<-c('biotin_bulk_HEK293T','200cell293Tb1','200cell293Tb2','400cellgm12878b1','400cellgm12878b2','4dn_bulk_gm12878')
subset_col<-as.vector(seq(1,length(names(batch_f)))*5)

#subset_col<-as.vector(1+seq(1,length(names(batch_f)))*4)
typeof(df_merge)
df_ins<-df_merge[, .SD, .SDcols = subset_col]
df_ins[is.na(df_ins)]<-0
colnames(df_ins)<-names(batch_f)
#df_ins<-df_ins[df_ins$biotin_bulk_hek_5Mcell!=0&df_ins$bulkHiChew_5k_hek_b1!=0&df_ins$bulkHiChew_5k_hek_b2!=0&df_ins$bulkpcr_50_hek!=0&df_ins$bulkpcr_500_hek!=0,]
corr <- round(cor(df_ins), 2)

col_fun <- colorRamp2(c(0.4, 0.65,0.9), c("royalblue3", "white",'indianred3'))

typeof(corr)
#corr[corr>0.9]<-0.9
breaks <- seq(0.4, 0.9, length.out = 300)

#color_scale =c("royalblue3", "grey", "indianred3",'darkred')
library(ComplexHeatmap)
pdf('bulkHiChew_cor_ins_heatmap.pdf',height = 8,width=8)
Heatmap(corr,column_title = paste0('bulkHiChew ins correlation'),
        cluster_rows = F,cluster_columns = F,col=col_fun,
        layer_fun = function(j, i, x, y, width, height, fill) {
          # since grid.text can also be vectorized
          grid.text(sprintf("%.2f", pindex(corr, i, j)), x, y, 
                    gp = gpar(fontsize = 10))})
dev.off()

#svg('bulkHiChew_cor_ins_heatmap.svg',height = 8,width=8)
#Heatmap(corr,column_title = paste0('bulkHiChew eigen correlation'),cluster_rows = F,cluster_columns = F,col=colorRampPalette(color_scale)(150),
#        layer_fun = function(j, i, x, y, width, height, fill) {
#          # since grid.text can also be vectorized
#          grid.text(sprintf("%.2f", pindex(corr, i, j)), x, y, 
#                    gp = gpar(fontsize = 10))})
#dev.off()





batch_f<- list(biotin_bulk_hek_5Mcell="hek293_ngsHiC_700M",
               bulkHiChew_5k_hek_b1="bulkHiChew_5000_1",bulkHiChew_5k_hek_b2="bulkHiChew_5000_2",
               bulkHiChew_500_hek_b1="bulkHiChew_500_1",bulkHiChew_500_hek_b2="bulkHiChew_500_2",
               bulkHiChew_50_hek_b1="bulkHiChew_50_1",bulkHiChew_50_hek_b2="bulkHiChew_50_2",
               bulkpcr_5k_hek_b1="bulkPCR_5000_1",
               bulkpcr_500_hek_b1="bulkPCR_500_1",
               bulkpcr_50_hek_b1="bulkPCR_50_1")
batch_data<- map(batch_f,~fread(paste0(.x,'.mcool.res_50000.window_150k.ins')))
ins_bin=50000
library(tidyverse)

#merge all data frames in list
df_merge<-batch_data %>% reduce(full_join, by=c('chrom','start'))
df_merge[is.na(df_merge)]<-0
#df_merge<-df_merge[complete.cases(df_merge), ]
#colnames(df_ins)<-c('biotin_bulk_HEK293T','200cell293Tb1','200cell293Tb2','400cellgm12878b1','400cellgm12878b2','4dn_bulk_gm12878')
subset_col<-as.vector(seq(1,length(names(batch_f)))*5)

#subset_col<-as.vector(1+seq(1,length(names(batch_f)))*4)
typeof(df_merge)
df_ins<-df_merge[, .SD, .SDcols = subset_col]
df_ins[is.na(df_ins)]<-0
colnames(df_ins)<-names(batch_f)
df_ins$x<-as.integer(rownames(df_ins))*ins_bin
df_ins_long <- reshape2::melt(df_ins, id.vars = "x", variable.name = "sample", value.name = "y")
#df_ins_long_sub<-df_ins_long[df_ins_long$x>39000000&df_ins_long$x<50000000,]
df_ins_long_sub<-df_ins_long[df_ins_long$x>41000000&df_ins_long$x<43000000,]

# colors <- c("darkred", "darkred", "darkred","darkgreen", "darkgreen", "darkgreen")

p2<-ggplot(df_ins_long_sub, aes(x = x, y = y,group = sample)) +
  geom_area(alpha=0.5, fill='grey30',color='grey30')+ facet_wrap(~ sample, scales = "free_y",ncol=1) +
  labs(title = "insulation score, chr7:41-43M, 50k bin",x = "Genomic position", y = "ins") + theme_classic()

pdf('bulkHiChew_cor_ins_track.pdf',height = 10,width=6)
p2
dev.off()

#svg('bulkHiChew_cor_ins_track.svg',height = 10,width=6)
#p2
#dev.off()




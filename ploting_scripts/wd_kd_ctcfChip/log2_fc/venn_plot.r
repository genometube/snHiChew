library(eulerr)
fit2 <- euler(c("HEK_WT" = 23766,
                "HEK_CTCFKD" = 200,
                "HEK_WT&HEK_CTCFKD" = 30024))
p1<-plot(fit2,fills = c('lightblue','indianred2','lightgrey'),quantities = TRUE,lty = 1:3,labels = list(font = 10,cex=1),shape = "ellipse",
         main=list(label="", cex = 3))

fit2 <- euler(c("HEK_WT" = 26304,
                "HEK_CTCFKD" = 209,
                "HEK_WT&HEK_CTCFKD" = 30901))
p2<-plot(fit2,fills = c('lightblue','indianred2','lightgrey'),quantities = TRUE,lty = 1:3,labels = list(font = 10,cex=1),shape = "ellipse",
         main=list(label="", cex = 3))

pdf('venn_1.pdf')      
p1
dev.off()

pdf('venn_2.pdf')      
p2
dev.off()

###############################################################################################
dpnI_ctcf_rnapII<-read.table('dpnI_ctcf_rnapII.multiinter.xls')

length(which(dpnI_ctcf_rnapII$V5=='1,2,3'))
length(which(dpnI_ctcf_rnapII$V5=='1,2'))

dpnI<-length(which(dpnI_ctcf_rnapII$V5=='1'))
ctcf<-length(which(dpnI_ctcf_rnapII$V5=='2'))
rnapII<-length(which(dpnI_ctcf_rnapII$V5=='3'))
dpnI_ctcf<-length(which(dpnI_ctcf_rnapII$V5=='1,2'))
dpnI_rnapII<-length(which(dpnI_ctcf_rnapII$V5=='1,3'))
ctcf_rnapII<-length(which(dpnI_ctcf_rnapII$V5=='2,3'))
dpnI_ctcf_rnapII<-length(which(dpnI_ctcf_rnapII$V5=='1,2,3'))

fit2 <- euler(c("dpnI" = dpnI,
                "ctcf" = ctcf,
                "rnapII" = rnapII,
                "dpnI&ctcf" = dpnI_ctcf,
                "dpnI&rnapII" = dpnI_rnapII,
                "ctcf&rnapII" = ctcf_rnapII,
                "dpnI&ctcf&rnapII" = dpnI_ctcf_rnapII),shape = "ellipse")
pdf('dpnI_ctcf_rnapII.venn.pdf')      
plot(fit2,fills = c('lightblue','indianred2','lightgrey'),quantities = TRUE,lty = 1:3,labels = list(font = 20,cex=3),shape = "ellipse",
     main=list(label="", cex = 3))

dev.off()
# "rnapII" = rnapII-dpnI_rnapII-ctcf_rnapII+dpnI_ctcf_rnapII

library(UpSetR)
dpnI<-55429
ctcf<-95246
rnapII<-58499
dpnI_ctcf<-36799
dpnI_rnapII<-36097
ctcf_rnapII<-54233
dpnI_ctcf_rnapII<-31540

# Dataset
input <- c(
  c("dpnI" = dpnI,
    "ctcf" = ctcf,
    "rnapII" = rnapII,
    "dpnI&ctcf" = dpnI_ctcf,
    "dpnI&rnapII" = dpnI_rnapII,
    "ctcf&rnapII" = ctcf_rnapII,
    "dpnI&ctcf&rnapII" = dpnI_ctcf_rnapII)
)

pdf('dpnI_ctcf_rnapII.upset.pdf',width=15,height=15)
# Plot
upset(fromExpression(input), 
      nintersects = 7, 
      nsets = 3, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.7, 0.3),
      number.angles = 0, 
      #text.dpnIle = 2.5,
      text.scale = 3,
      point.size = 8, 
      line.size = 1.5,
      set_size.angles=45,
      #dpnIle.sets='identity'
)
dev.off()
?upset





########################################################################

fit2 <- euler(c("dpnI" = dpnI-dpnI_ctcf-dpnI_rnapII+dpnI_ctcf_rnapII,
                "ctcf" = ctcf-dpnI_ctcf-ctcf_rnapII+dpnI_ctcf_rnapII,
                "rnapII" = rnapII-ctcf_rnapII,
                "dpnI&ctcf" = dpnI_ctcf-dpnI_ctcf_rnapII,
                "dpnI&rnapII" = dpnI_rnapII-dpnI_ctcf_rnapII,
                "ctcf&rnapII" = ctcf_rnapII-dpnI_ctcf_rnapII,
                "dpnI&ctcf&rnapII" = dpnI_ctcf_rnapII),shape = "ellipse")
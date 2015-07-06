setwd("D:/ONC-Analysis/users/ash/michelle/michelle_tcga")
load("D:/ONC-Analysis/users/ash/michelle/methylation_analysis/R_image/probe.to.gene.RData")
load("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/r_images/data_LUAD.RData")
#load("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/r_images/LUSC.RData")
load("D:/ONC-Analysis/users/ash/michelle/methylation_analysis/R_image/probesWithin5000TSS.RData")

#load function
source("D:/ONC-Analysis/users/ash/michelle/methylation_analysis/functions/tcga_functions.R")

#Remove probes from X and Y Chromosome
#data<-data[rownames(data) %in% probe.to.gene$probe,]
dim(data)




#get 15m data
stable_hyper <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.stable.hyper_unt_15-unt_1-csc_15.csv", row.names=1)
stable_hypo <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.stable.hypo_unt_15-unt_1-csc_15.csv", row.names=1)
inter_hyper <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.inter.hyper_unt_15-unt_1-csc_15.csv", row.names=1)
inter_hypo <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.inter.hypo_unt_15-unt_1-csc_15.csv", row.names=1)

tcga.stable.hyper <- data[rownames(data) %in% rownames(stable_hyper),]
tcga.stable.hypo <- data[rownames(data) %in% rownames(stable_hypo),]
tcga.inter.hyper <- data[rownames(data) %in% rownames(inter_hyper),]
tcga.inter.hypo <- data[rownames(data) %in% rownames(inter_hypo),]

dim(tcga.stable.hyper)
dim(tcga.stable.hypo )
dim(tcga.inter.hyper)
dim(tcga.inter.hypo)

stable_hyper$Gene <- NULL
stable_hypo$Gene <- NULL
inter_hyper$Gene <- NULL
inter_hypo$Gene <- NULL


#to get what sample is what
#sample <- sample_type("LUAD_450k.csv")
#Plot heatmaps
#LUAD
#heat_tcga(tcga.stable.hyper, plotName="15m-csc-tcga-methylation-heatmap_stable_Hyper_LUAD.jpg", title="Heatmap stable Hypermethylation LUAD")
#heat_tcga(tcga.stable.hypo, plotName="15m-csc-tcga-methylation-heatmap_stable_Hypo_LUAD.jpg", title="Heatmap stable Hypomethylation LUAD")
#heat_tcga(tcga.inter.hyper, plotName="15m-csc-tcga-methylation-heatmap_inter_Hyper_LUAD.jpg", title="Heatmap intermediate Hypermethylation LUAD")
#heat_tcga(tcga.inter.hypo, plotName="15m-csc-tcga-methylation-heatmap_inter_Hypo_LUAD.jpg", title="Heatmap intermediate Hypomethylation LUAD")

#to get what sample is what
#sample <- sample_type("LUSC_450k.csv")

#LUSC
heat_tcga(tcga.stable.hyper, plotName="15m-csc-tcga-methylation-heatmap_stable_Hyper_LUSC.jpg", title="Heatmap stable Hypermethylation LUSC")
heat_tcga(tcga.stable.hypo, plotName="15m-csc-tcga-methylation-heatmap_stable_Hypo_LUSC.jpg", title="Heatmap stable Hypomethylation LUSC")
heat_tcga(tcga.inter.hyper, plotName="15m-csc-tcga-methylation-heatmap_inter_Hyper_LUSC.jpg", title="Heatmap intermediate Hypermethylation LUSC")
heat_tcga(tcga.inter.hypo, plotName="15m-csc-tcga-methylation-heatmap_inter_Hypo_LUSC.jpg", title="Heatmap intermediate Hypomethylation LUSC")




########################################################
#m<- get_row_order_heatmap(tcga.stable.hypo, stable_hypo)

sample<- c("4","3","3","4","3","4","2","3","4")

#LUAD
heat_tcga_col(get_row_order_heatmap(tcga.stable.hyper, stable_hyper), plotName="stable_hyper_15m-LUAD.jpg", title="Heatmap for probes stable Hypermethylated in CSC 15m-LUAD")
heat_tcga_col(get_row_order_heatmap(tcga.stable.hypo, stable_hypo), plotName="stable_hypo_15m-LUAD.jpg", title="Heatmap for probes stable Hypomethylated in CSC 15m-LUAD")
heat_tcga_col(get_row_order_heatmap(tcga.inter.hyper, inter_hyper), plotName="inter_hyper_15m-LUAD.jpg", title="Heatmap for probes intermediate Hypermethylated in CSC 15m-LUAD")
heat_tcga_col(get_row_order_heatmap(tcga.inter.hypo, inter_hypo), plotName="inter_hypo_15m-LUAD.jpg", title="Heatmap for probes intermediate Hypomethylated in CSC 15m-LUAD")

#LUSC
#heat_tcga_col(get_row_order_heatmap(tcga.stable.hyper, stable_hyper), plotName="stable_hyper_15m-LUSC.jpg", title="Heatmap for probes stable Hypermethylated in CSC 15m-LUSC")
#heat_tcga_col(get_row_order_heatmap(tcga.stable.hypo, stable_hypo), plotName="stable_hypo_15m-LUSC.jpg", title="Heatmap for probes stable Hypomethylated in CSC 15m-LUSC")
#heat_tcga_col(get_row_order_heatmap(tcga.inter.hyper, inter_hyper), plotName="inter_hyper_15m-LUSC.jpg", title="Heatmap for probes intermediate Hypermethylated in CSC 15m-LUSC")
#heat_tcga_col(get_row_order_heatmap(tcga.inter.hypo, inter_hypo), plotName="inter_hypo_15m-LUSC.jpg", title="Heatmap for probes intermediate Hypomethylated in CSC 15m-LUSC")


#odata=original data of which heatmap is made. rdata= data that is to be reordered
get_row_order_heatmap <- function (odata, rdata){
  library(gplots)
  test <- heatmap.2(as.matrix(odata))
  z<-odata[rev(test$rowInd), test$colInd]
  reordered <- rdata[rownames(z),,drop=FALSE]
  
  return(reordered)
  
}



















setwd("D:/ONC-Analysis/users/ash/michelle/michelle_tcga")
source("D:/ONC-Analysis/users/ash/michelle/methylation_analysis/functions/tcga_functions.R")

#get 15m data
stable_hyper <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.stable.hyper_unt_15-unt_1-csc_15.csv", row.names=1)
stable_hypo <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.stable.hypo_unt_15-unt_1-csc_15.csv", row.names=1)
inter_hyper <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.inter.hyper_unt_15-unt_1-csc_15.csv", row.names=1)
inter_hypo <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.inter.hypo_unt_15-unt_1-csc_15.csv", row.names=1)

#modify table to be used by heatmap.2

stable_hyper$Gene <- NULL
stable_hypo$Gene <- NULL
inter_hyper$Gene <- NULL
inter_hypo$Gene <- NULL

# unt=3(green) csc=4(blue) 15mcsc=2 (red)
sample<- c("4","3","3","4","3","4","2","3","4")

heat_tcga(stable_hyper, plotName="stable_hyper_15m.jpg", title="Heatmap for probes stable Hypermethylated in CSC 15m")
heat_tcga(stable_hypo, plotName="stable_hypo_15m.jpg", title="Heatmap for probes stable Hypomethylated in CSC 15m")
heat_tcga(inter_hyper, plotName="inter_hyper_15m.jpg", title="Heatmap for probes intermediate Hypermethylated in CSC 15m")
heat_tcga(inter_hypo, plotName="inter_hypo_15m.jpg", title="Heatmap for probes intermediate Hypomethylated in CSC 15m")

################
#combine
c.hyper <- rbind.data.frame(stable_hyper,inter_hyper)
c.hypo <-  rbind.data.frame(stable_hypo,inter_hypo)

c.hyper$Gene <- NULL
c.hypo$Gene <- NULL


c.hyper$X <- NULL
c.hypo$X <- NULL


rownames(c.hyper) <- NULL
rownames(c.hypo) <- NULL


sample<- c("4","3","3","4","3","4","2","3","4")

heat_tcga(c.hyper, plotName="combine_hyper_15m.jpg", title="Heatmap for probes combine Hypermethylated in CSC 15m")
heat_tcga(c.hypo, plotName="combine_hypo_15m.jpg", title="Heatmap for probes combine Hypomethylated in CSC 15m")

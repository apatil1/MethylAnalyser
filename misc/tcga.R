

setwd("D:/ONC-Analysis/users/ash/michelle/michelle_tcga")
load("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/michelle_tcga_luad.RD)

#load function
source("D:/ONC-Analysis/users/ash/michelle/methylation_analysis/functions/tcga_functions.R")

#path to the TCGA files
path <- "D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
type <- "LUAD_450k"

path <- "D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUSC/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
type <- "LUSC_450k"

#Read data from file
data <- read_level3(path, type)

#Remove NA elements
data <- na.omit(data)
head(data)
dim(data)

save.image("LUSC.RData")

###################################################################################

#Remove probes from X and Y Chromosome
data<-data[rownames(data) %in% probe.to.chr$Probe_ID,]
dim(data)

#to get what sample is what
sample <- sample_type("LUSC_450k.csv")

##################################################################################

#get 15m data
stable_hyper <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.stable.hyper_unt_15-unt_1-csc_15.csv")
stable_hypo <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.stable.hypo_unt_15-unt_1-csc_15.csv")
inter_hyper <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.inter.hyper_unt_15-unt_1-csc_15.csv")
inter_hypo <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.inter.hypo_unt_15-unt_1-csc_15.csv")

tcga.stable.hyper <- data[rownames(data) %in% stable_hyper$X,]
tcga.stable.hypo <- data[rownames(data) %in% stable_hypo$X,]
tcga.inter.hyper <- data[rownames(data) %in% inter_hyper$X,]
tcga.inter.hypo <- data[rownames(data) %in% inter_hypo$X,]

dim(tcga.stable.hyper)
dim(tcga.stable.hypo )
dim(tcga.inter.hyper)
dim(tcga.inter.hypo)

#combine michelle and TCGA data
tcga.stable.hyper <- cbind(tcga.stable.hyper, stable_hyper[stable_hyper$X %in% rownames(tcga.stable.hyper), 2:10])
tcga.stable.hypo <- cbind(tcga.stable.hypo, stable_hypo[stable_hypo$X %in% rownames(tcga.stable.hypo), 2:10])
tcga.inter.hyper <- cbind(tcga.inter.hyper, inter_hyper[inter_hyper$X %in% rownames(tcga.inter.hyper), 2:10])
tcga.inter.hypo <- cbind(tcga.inter.hypo, inter_hypo[inter_hypo$X %in% rownames(tcga.inter.hypo), 2:10])

#add 9 samples from michelle data to sample for identification michelle: unt=3(green) csc=4(blue), tumor=1 (black), recurrent tumor=2 (red), normal=0 (white)
sample <- append(sample, c(4,3,3,4,3,4,4,3,4))

#Plot heatmaps
heat_tcga(tcga.stable.hyper, plotName="heatmap_stable_Hyper_LUSC.jpg")
heat_tcga(tcga.stable.hypo, plotName="heatmap_stable_Hypo_LUSC.jpg")
heat_tcga(tcga.inter.hyper, plotName="heatmap_inter_Hyper_LUSC.jpg")
heat_tcga(tcga.inter.hypo, plotName="heatmap_inter_Hypo_LUSC.jpg")





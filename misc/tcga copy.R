

setwd("D:/ONC-Analysis/users/ash/michelle/michelle_tcga")
#load("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/michelle_tcga_luad.RD")

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


##################################################################################

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

#combine michelle and TCGA data
tcga.stable.hyper <- merge(tcga.stable.hyper,stable_hyper ,by="row.names")
tcga.stable.hypo <- merge(tcga.stable.hypo,stable_hypo ,by="row.names")
tcga.inter.hyper <- merge(tcga.inter.hyper,inter_hyper ,by="row.names")
tcga.inter.hypo <- merge(tcga.inter.hypo,inter_hypo ,by="row.names")


#modify table to be used by heatmap.2
rownames(tcga.stable.hyper) <- tcga.stable.hyper$Row.names
tcga.stable.hyper$Row.names <- NULL
tcga.stable.hyper$Gene <- NULL

rownames(tcga.stable.hypo) <- tcga.stable.hypo$Row.names
tcga.stable.hypo$Row.names <- NULL
tcga.stable.hypo$Gene <- NULL

rownames(tcga.inter.hyper) <- tcga.inter.hyper$Row.names
tcga.inter.hyper$Row.names <- NULL
tcga.inter.hyper$Gene <- NULL

rownames(tcga.inter.hypo) <- tcga.inter.hypo$Row.names
tcga.inter.hypo$Row.names <- NULL
tcga.inter.hypo$Gene <- NULL
##################################################################################################################################
#only 15 csc
tcga.stable.hyper <- tcga.stable.hyper[, !colnames(tcga.stable.hyper) %in% c("csc_10", "unt_10" ,"unt_1" , "csc_1",  "unt_6" , "csc_6", "unt_15","csc20" )]
tcga.stable.hypo <- tcga.stable.hypo[, !colnames(tcga.stable.hypo) %in% c("csc_10", "unt_10" ,"unt_1" , "csc_1",  "unt_6" , "csc_6", "unt_15","csc20" )]
tcga.inter.hyper <- tcga.inter.hyper[, !colnames(tcga.inter.hyper) %in% c("csc_10", "unt_10" ,"unt_1" , "csc_1",  "unt_6" , "csc_6", "unt_15","csc20" )]
tcga.inter.hypo <- tcga.inter.hypo[, !colnames(tcga.inter.hypo) %in% c("csc_10", "unt_10" ,"unt_1" , "csc_1",  "unt_6" , "csc_6", "unt_15","csc20" )]

#to get what sample is what
sample <- sample_type("LUAD_450k.csv")
#add 9 samples from michelle data to sample for identification michelle: unt=3(green) csc=4(blue), tumor=1 (black), recurrent tumor=1 (black), normal=0 (white)
sample <- append(sample, c(2))

#only 15 csc and 15 unt
tcga.stable.hyper <- tcga.stable.hyper[, !colnames(tcga.stable.hyper) %in% c("csc_10", "unt_10" ,"unt_1" , "csc_1",  "unt_6" , "csc_6", "csc20" )]
tcga.stable.hypo <- tcga.stable.hypo[, !colnames(tcga.stable.hypo) %in% c("csc_10", "unt_10" ,"unt_1" , "csc_1",  "unt_6" , "csc_6", "csc20" )]
tcga.inter.hyper <- tcga.inter.hyper[, !colnames(tcga.inter.hyper) %in% c("csc_10", "unt_10" ,"unt_1" , "csc_1",  "unt_6" , "csc_6", "csc20" )]
tcga.inter.hypo <- tcga.inter.hypo[, !colnames(tcga.inter.hypo) %in% c("csc_10", "unt_10" ,"unt_1" , "csc_1",  "unt_6" , "csc_6", "csc20" )]

#to get what sample is what
sample <- sample_type("LUAD_450k.csv")
#add 9 samples from michelle data to sample for identification michelle: unt=3(green) csc=4(blue), tumor=1 (black), recurrent tumor=1 (black), normal=0 (white)
sample <- append(sample, c(2,3))



#Plot heatmaps
#LUAD
heat_tcga(tcga.stable.hyper, plotName="15m-csc-tcga-methylation-heatmap_stable_Hyper_LUAD.jpg", title="Heatmap stable Hypermethylation LUAD")
heat_tcga(tcga.stable.hypo, plotName="15m-csc-tcga-methylation-heatmap_stable_Hypo_LUAD.jpg", title="Heatmap stable Hypomethylation LUAD")
heat_tcga(tcga.inter.hyper, plotName="15m-csc-tcga-methylation-heatmap_inter_Hyper_LUAD.jpg", title="Heatmap intermediate Hypermethylation LUAD")
heat_tcga(tcga.inter.hypo, plotName="15m-csc-tcga-methylation-heatmap_inter_Hypo_LUAD.jpg", title="Heatmap intermediate Hypomethylation LUAD")

heat_tcga(tcga.stable.hyper, plotName="15m-csc-unt-tcga-methylation-heatmap_stable_Hyper_LUAD.jpg", title="Heatmap stable Hypermethylation LUAD")
heat_tcga(tcga.stable.hypo, plotName="15m-csc-unt-tcga-methylation-heatmap_stable_Hypo_LUAD.jpg", title="Heatmap stable Hypomethylation LUAD")
heat_tcga(tcga.inter.hyper, plotName="15m-csc-unt-tcga-methylation-heatmap_inter_Hyper_LUAD.jpg", title="Heatmap intermediate Hypermethylation LUAD")
heat_tcga(tcga.inter.hypo, plotName="15m-csc-unt-tcga-methylation-heatmap_inter_Hypo_LUAD.jpg", title="Heatmap intermediate Hypomethylation LUAD")

#LUSC
#heat_tcga(tcga.stable.hyper, plotName="15m-csc-tcga-methylation-heatmap_stable_Hyper_LUSC.jpg", title="Heatmap stable Hypermethylation LUSC")
#heat_tcga(tcga.stable.hypo, plotName="15m-csc-tcga-methylation-heatmap_stable_Hypo_LUSC.jpg", title="Heatmap stable Hypomethylation LUSC")
#heat_tcga(tcga.inter.hyper, plotName="15m-csc-tcga-methylation-heatmap_inter_Hyper_LUSC.jpg", title="Heatmap intermediate Hypermethylation LUSC")
#heat_tcga(tcga.inter.hypo, plotName="15m-csc-tcga-methylation-heatmap_inter_Hypo_LUSC.jpg", title="Heatmap intermediate Hypomethylation LUSC")

#heat_tcga(tcga.stable.hyper, plotName="15m-csc-unt-tcga-methylation-heatmap_stable_Hyper_LUSC.jpg", title="Heatmap stable Hypermethylation LUSC")
#heat_tcga(tcga.stable.hypo, plotName="15m-csc-unt-tcga-methylation-heatmap_stable_Hypo_LUSC.jpg", title="Heatmap stable Hypomethylation LUSC")
#heat_tcga(tcga.inter.hyper, plotName="15m-csc-unt-tcga-methylation-heatmap_inter_Hyper_LUSC.jpg", title="Heatmap intermediate Hypermethylation LUSC")
#heat_tcga(tcga.inter.hypo, plotName="15m-csc-unt-tcga-methylation-heatmap_inter_Hypo_LUSC.jpg", title="Heatmap intermediate Hypomethylation LUSC")

########################
#Random -500 to +1500 TSS

dim(stable_hyper)
dim(stable_hypo )
dim(inter_hyper)
dim(inter_hypo)
#stble_hyper= 431 stable hypo=88 inter hyper= 1207 inter hypo= 608

s <- c(431,88,1207,608)
t <- c("stable hypermethylation", "stable hypomethylation", "intermediate hypermethylation", "intermediate hypomethylation")
dim(data_CpGI.TSS)

for (i in 1:4) {
  tcga.random <- data_CpGI.TSS[sample(c(1:nrow(data_CpGI.TSS)), size=s[i]), ]
  dim(tcga.random)
  
  heat_tcga(tcga.random, plotName=paste("Random LUAD ", t[i],".jpg",sep=""), title=paste("Random LUAD ", t[i],sep=""))
  #heat_tcga(tcga.random, plotName=paste("Random LUSC ", t[i],".jpg",sep=""), title=paste("Random LUSC ", t[i],sep=""))
}
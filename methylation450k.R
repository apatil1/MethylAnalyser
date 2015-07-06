#Load libraries
library(limma)
library(gplots)
library(heatmap.plus)
library(IlluminaHumanMethylation450k.db)
library(minfi)
library(IlluminaHumanMethylation450kmanifest)

#Load images
load("~/Desktop/methylation_analysis/R_image/all.probes.hg19.RData")
load("/Users/ashwini/Desktop/methylation_analysis/R_image/cgi.xy.filter.RData")
load("/Users/ashwini/Desktop/methylation_analysis/R_image/probesWithin5000TSS.RData")
load("~/Desktop/methylation_analysis/R_image/probe.to.gene.RData")
load("~/Desktop/methylation_analysis/R_image/annotation.RData")

#load functions
source("/Users/ashwini/Desktop/methylation_analysis/functions/methylation450k/methylation450k_functions.R")

#Get paths to input files and folders
idat <- "/Users/ashwini/Desktop/michelle_data/idat_files/"
sheet_genelist <- "/Users/ashwini/Desktop/new-annotation/sheet_genelist.csv"
sheet_deltaplot <- "/Users/ashwini/Desktop/new-annotation/sheet_deltaplot.csv"
folder <- "/Users/ashwini/Desktop/new-annotation/"

#Fuction for reading and normalizing data
beta.tab <- data_methylation(idat, folder)


#Retain the probes which are +/-1500 from TSS (loaded from annotation.Rdata image (for file annotation.R))
#Subset the beta.tab table to retain CpG-island probes.
dim(beta.tab)
head(beta.tab)
beta.tab_CpGI <- beta.tab[rownames(beta.tab) %in% cgi.xy.filter, ]
beta.tab_CpGI <- na.omit(beta.tab_CpGI)
dim(beta.tab_CpGI)

########################################################
#Subset the beta.tab table to retain gene body probes.
body <- beta.tab[names(body.probes), ]
body <- na.omit(body)
dim(body)

#Subset the beta.tab table to retain shelf probes.
shelf <- beta.tab[names(shelf.probes), ]
shelf <- na.omit(shelf)
dim(shelf)

#Subset the beta.tab table to retain shore probes.
shore <- beta.tab[names(shore.probes), ]
shore <- na.omit(shore)
dim(shore)
##########################################################

#get +/- 1500 TSS probes
head(probesWithin5000TSS)
tss.1500.probes <- as.character(probesWithin5000TSS[abs(probesWithin5000TSS$Dist.to.TSS) <= 1500, "ProbeName"])
beta.tab_CpGI.PlusMinus1500_TSS <- beta.tab_CpGI[rownames(beta.tab_CpGI) %in% tss.1500.probes,] #these are the CpGI probes within +/-1500 bp from TSS
beta.tab_CpGI.PlusMinus1500_TSS <- na.omit(beta.tab_CpGI.PlusMinus1500_TSS)
beta.tab_CpGI<-data.frame(beta.tab_CpGI)
beta.tab_CpGI.PlusMinus1500_TSS<-data.frame(beta.tab_CpGI.PlusMinus1500_TSS)
dim(beta.tab_CpGI.PlusMinus1500_TSS)

#Write TSS +/-1500 probes to file with gene names
if (length(rownames(beta.tab_CpGI.PlusMinus1500_TSS)) > 0) {
  beta.tab_CpGI.PlusMinus1500_TSS <- cbind.data.frame(beta.tab_CpGI.PlusMinus1500_TSS, "Gene"=rep(NA, times=nrow(beta.tab_CpGI.PlusMinus1500_TSS)))
  for(j in 1:nrow(beta.tab_CpGI.PlusMinus1500_TSS)) {
    probe.id <- rownames(beta.tab_CpGI.PlusMinus1500_TSS)[j]
    gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
    if(length(gene) == 1) { beta.tab_CpGI.PlusMinus1500_TSS$Gene[j] <- as.character(gene) }
  }
  write.csv(beta.tab_CpGI.PlusMinus1500_TSS,paste(folder, "tss_probes.csv",sep="") , quote=F, row.names=T) #row.names set to T to get the probe ids.
}


#Analyze methylation data and create genelist

analysis_methylation(beta.tab_CpGI.PlusMinus1500_TSS, sheet_genelist, folder, "TSS") 
analysis_methylation(body, sheet_genelist, folder, "Body") 
analysis_methylation(shore, sheet_genelist, folder, "Shore") 
analysis_methylation(shelf, sheet_genelist, folder, "Shelf") 

#Make delta plots
delta_plot(beta.tab_CpGI.PlusMinus1500_TSS, sheet_deltaplot, folder, "TSS") 
delta_plot(body, sheet_deltaplot, folder, "Body") 
delta_plot(shelf, sheet_deltaplot, folder, "Shelf") 
delta_plot(shore, sheet_deltaplot, folder, "Shore") 

#save image
save(beta.tab_CpGI,beta.tab_CpGI.PlusMinus1500_TSS,file="michelle_new-annotation.RData")


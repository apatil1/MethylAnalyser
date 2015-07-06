
setwd("D:/ONC-Analysis/users/ash/michelle/LUSC_expression/Expression-Genes/UNC__AgilentG4502A_07_3/Level_1")
path <- "D:/ONC-Analysis/users/ash/michelle/LUSC_expression/Expression-Genes/UNC__AgilentG4502A_07_3/Level_1/"

#create targets file
make_target(path)

library(limma)
targets <- readTargets("targets.txt", sep=" ")
RG <- read.maimages(targets, source="agilent", path="D:/ONC-Analysis/users/ash/michelle/LUSC_expression/Expression-Genes/UNC__AgilentG4502A_07_3/Level_1/")

RG$weights <-  matrix(rep(RG$genes$ControlType,ncol(RG$R)),ncol=ncol(RG$R),byrow=F)
RG$weights[RG$genes$ControlType!=0,] <- 0
RG$weights[RG$genes$ControlType==0,] <- 1

#normalization
#normwithin <-normalizeWithinArrays(RG,method='loess',bc.method='none')
normwithin <-normalizeWithinArrays(RG,method='loess',bc.method='normexp', offset=50)
normbetween <-normalizeBetweenArrays(normwithin,method='Aquantile')

#Remove controls from normwithin/between
normwithin <- normwithin[normbetween$genes$ControlType==0,]
normbetween <- normbetween[normbetween$genes$ControlType==0,]

#make table of gene names and expression values
tcga <- cbind.data.frame(GeneName=normbetween$genes$GeneName, normbetween$M)
tcga.table <- aggregate(. ~  GeneName, data=tcga, median)
colnames(tcga.table) <- c("Gene", targets$Cy5)

head(tcga.table)
dim(tcga.table)

#get up and down genes from 15m expression data and pull the same genes from TCGA data
up_15m <- read.csv("D:/ONC-Analysis/users/ash/michelle/up_MV10-15m.csv")
down_15m <- read.csv("D:/ONC-Analysis/users/ash/michelle/down_MV10-15m.csv")

dim(up_15m)
dim(down_15m)
#make heatmaps of these genes

tcga.table.up <- tcga.table[as.character(tcga.table$Gene) %in% as.character(up_15m$Gene), ]
#tcga.table.down <- tcga.table[tcga.table$Gene %in% down_15m$Gene, ]

up <- merge(up_15m,tcga.table ,by="Gene")
down <- merge(down_15m,tcga.table ,by="Gene")

rownames(up) <- up$Gene
up$Gene <- NULL


rownames(down) <- up$Gene
down$Gene <- NULL


dim(up)
dim(down)
#load function
source("D:/ONC-Analysis/users/ash/michelle/methylation_analysis/functions/tcga_functions.R")

sample <- c("3", rep("1",length(colnames(up))-1))
heat_tcga(up, plotName="up_agilentLUSC.jpg", title="Heatmap up agilent LUSC")

save.image( "tcga_expression_agilent_LUAD.RData")


#function to create target file Agilent TCGA data
make_target <- function(path) {

  FileName <- list.files(path)
  Cy3 <- rep("Mock", length(file))
  Cy5 <- 1:length(FileName)

  target <-cbind(FileName,Cy3,Cy5)
  write.table(target, "Targets.txt", row.names=F)

}
####################################################
up_15m <- read.csv("D:/ONC-Analysis/users/ash/michelle/up_MV10-15m.csv", row.names=1)
down_15m <- read.csv("D:/ONC-Analysis/users/ash/michelle/down_MV10-15m.csv", row.names=1)

dim(up_15m)
dim(down_15m)

up <- tcga.table[tcga.table$Gene %in% rownames(up_15m),]
down <- tcga.table[tcga.table$Gene %in% rownames(down_15m),]

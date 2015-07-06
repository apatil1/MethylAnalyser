
#load function
source("D:/ONC-Analysis/users/ash/michelle/methylation_analysis/functions/tcga_functions.R")

setwd("D:/ONC-Analysis/users/ash/michelle")
#RNAseq V2
#read filenames from folder

#LUSC
path_folder <-"D:/ONC-Analysis/users/ash/michelle/lusc_rnaseq/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/"
path_sdrf <- "D:/ONC-Analysis/users/ash/michelle/lusc_rnaseq/METADATA/UNC__IlluminaHiSeq_RNASeqV2/unc.edu_LUSC.IlluminaHiSeq_RNASeqV2.1.11.0.sdrf.txt"
data <- rnaseqv2_read(path_folder, path_sdrf)

#LUAD
path_folder <-"D:/ONC-Analysis/users/ash/michelle/luad_rnaseq/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/"
path_sdrf <- "D:/ONC-Analysis/users/ash/michelle/luad_rnaseq/METADATA/UNC__IlluminaHiSeq_RNASeqV2/unc.edu_LUAD.IlluminaHiSeq_RNASeqV2.1.14.0.sdrf.txt"
data <- rnaseqv2_read(path_folder, path_sdrf)


#function to reat rnaseq v2 level 3 data
rnaseqv2_read <- function(path_folder,path_sdrf) {
  
  file <- list.files(path_folder)
  
  #get gene.normalized_results file from it
  chars="genes.normalized_results"
  file_names <- file[grepl( chars,file)]
  
  #read sdrf file
  sdrf <- read.table(path_sdrf, sep="\t", header=T)
  
  #read these files
  for (i in 1:length(file_names)) {
    read <- read.table(paste(path_folder, file_names[i], sep=""), sep="\t", header=T, stringsAsFactors=F)
    
    if (i==1) {
      gene <- strsplit(as.character(read$gene_id[1]),"[|]")[[1]][1]
      for (j in 2:length(read$gene_id)) {
        gene <- append(gene, strsplit(as.character(read$gene_id[j]),"[|]")[[1]][1])
      }
      
      data <- cbind.data.frame(gene, read[,2])
      ex <- substr(file_names[i], 9, 44)
      colnames <- c("gene",as.character(sdrf[sdrf$Extract.Name==ex,2][1]))
    }
    else {
      data <- cbind.data.frame(data, read[,2])
      ex <- substr(file_names[i], 9, 44)
      colnames <- append(colnames, as.character(sdrf[sdrf$Extract.Name==ex,2][1]) )
    }
    colnames(data) <- colnames
  }
  
  return (data)
}


# TCGA barcode and parse the barcode to determine tumor and normal sample
colnames(data)

sample <- sample_type(colnames(data))


sample_type <- function(file) {
  
  #f <- read.csv(file)
  
  
  for (i in 2:length(file)) {
    
    #get sample infromation from TCGA bar code
    #ch1 <- strsplit(as.character(f$x[i]), "[.]")[[1]][6]
    ch2 <- substr(file[i], 14, 15)
    
    if (ch2=="11") { #normal sample
      if (i==2) { 
        sample <- 0 
      } else {
        sample <- append(sample, "0")
      }
      
    } else if (ch2=="01") { #primary solid tumor
      if (i==2) { 
        sample <- 1 
      } else {
        sample <- append(sample, "1")
      }
      
    } else if (ch2=="02") { #recurrent solid tumor
      if (i==2) { 
        sample <- 1
      } else {
        sample <- append(sample, "1")
      }
    }
  }
  
  head(sample)
  return(sample)
}






################################
x <- (data[,2:ncol(data)])
x.sub = replace(x, x<1, NA)
x.sub <- cbind.data.frame(gene=data$gene, x.sub)
x.sub.na <- x.sub[complete.cases(x.sub),]
x.log <- log2(x.sub.na[,2:ncol(x.sub.na)])
##################################
#topvaeiable genes

x.sd <- apply(x.log, 1, sd)
range(x.sd)
topvar <- x.log[x.sd >= 2.2, ]
dim(topvar)

#topvar$gene<-NULL
#luad
#heat_tcga_scale(topvar, plotName="scaletopvar_luad.jpg", title="Heatmap LUAD")
#heat_tcga(topvar, plotName="topvar_luad.jpg", title="Heatmap LUAD")

#lusc
heat_tcga_scale(topvar, plotName="scaletopvar_lusc.jpg", title="Heatmap LUSC")
heat_tcga(topvar, plotName="topvar_lusc.jpg", title="Heatmap LUSC")

###############################################
x.log.gene <- cbind.data.frame(gene=x.sub.na$gene, x.log)

up_15m <- read.csv("D:/ONC-Analysis/users/ash/michelle/up_MV10-15m.csv", row.names=1)
down_15m <- read.csv("D:/ONC-Analysis/users/ash/michelle/down_MV10-15m.csv", row.names=1)

tcga_up <- data[data$gene %in% rownames(up_15m),]
tcga_down <- data[data$gene %in% rownames(down_15m),]


tcga_up <-  x.log.gene[ x.log.gene$gene %in% rownames(up_15m),]
tcga_down <-  x.log.gene[ x.log.gene$gene %in% rownames(down_15m),]


tcga_up$gene <- NULL
tcga_down$gene <- NULL

#luaad
heat_tcga_scale(tcga_up, plotName="scale-15m-rnaseqv2-heatmap_tcga_up_LUAD.jpg", title="Heatmap RNAseqV2 tcga up LUAD")
heat_tcga_scale(tcga_down, plotName="scale-15m-rnaseqv2-heatmap_tcga_down_LUAD.jpg", title="Heatmap RNAseqV2 tcga down LUAD")
#heat_tcga(tcga_up, plotName="15m-rnaseqv2-heatmap_tcga_up_LUAD.jpg", title="Heatmap RNAseqV2 tcga up LUAD")
#heat_tcga(tcga_down, plotName="15m-rnaseqv2-heatmap_tcga_down_LUAD.jpg", title="Heatmap RNAseqV2 tcga down LUAD")

#lusc
heat_tcga_scale(tcga_up, plotName="scale-15m-rnaseqv2-heatmap_tcga_up_LUSC.jpg", title="Heatmap RNAseqV2 tcga up LUSC")
heat_tcga_scale(tcga_down, plotName="scale-15m-rnaseqv2-heatmap_tcga_down_LUSC.jpg", title="Heatmap RNAseqV2 tcga down LUSC")
heat_tcga(tcga_up, plotName="15m-rnaseqv2-heatmap_tcga_up_LUSC.jpg", title="Heatmap RNAseqV2 tcga up LUSC")
heat_tcga(tcga_down, plotName="15m-rnaseqv2-heatmap_tcga_down_LUSC.jpg", title="Heatmap RNAseqV2 tcga down LUSC")

###############################################

#########################################


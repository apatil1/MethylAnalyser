

#Function to read level 3 data from TCGA methylation data. Returns table of beta values combining all the samples.
read_level3 <- function(path, type) {
  
  file <- list.files(path) #gets list of files
  #write.csv(file, paste(type,".csv", sep=""))
  
  #reads data from all the files 
  for (i in 1:length(file)) {
  
    read <- read.table(paste(path, file[i], sep=""), sep="\t", header=T, skip=1)
    
    if (i==1) {
      data<-cbind(read[,2])
      rownames(data)<-read[,1]
      col<-paste("s", i, sep="")
    }
    else {
      data<-cbind(data,read[,2])
      col<-append(col,paste("s", i, sep=""))
    }
  }
  
  colnames(data)<-col
  #data <- na.omit(data) #omit all NA elements
  head(data)
  dim(data) 
  
  return (data)

}

#Function to create list of samples
sample_type <- function(file) {
  
  f <- read.csv(file)
  
  
  for (i in 1:length(f$x)) {
    
    #get sample infromation from TCGA bar code
    ch1 <- strsplit(as.character(f$x[i]), "[.]")[[1]][6]
    ch2 <- substr(ch1, 14, 15)
    
    if (ch2=="11") { #normal sample
      if (i==1) { 
        sample <- 0 
      } else {
        sample <- append(sample, "0")
      }
      
    } else if (ch2=="01") { #primary solid tumor
      if (i==1) { 
        sample <- 1 
      } else {
        sample <- append(sample, "1")
      }
      
    } else if (ch2=="02") { #recurrent solid tumor
      if (i==1) { 
        sample <- 1
      } else {
        sample <- append(sample, "1")
      }
    }
  }
  
  head(sample)
  return(sample)
}

#Function to plot heatmap
heat_tcga <- function(x, plotName, title) {  
  
  library(gplots)
  jpeg(plotName, width = 800, height = 800, quality=100)
  heatmap.2(as.matrix(x), col=heat.colors(16), trace="none",labCol=NA,density.info="none", hclustfun = function(x)hclust(x, method = "ward"),labRow=NA, margins=c(10,10), cexRow=0.8, cexCol=1.2, key=T, ColSideColors=sample, main=title)
  dev.off()
  
}
heat_tcga_col <- function(x, plotName, title) {  
  
  library(gplots)
  jpeg(plotName, width = 800, height = 800, quality=100)
  heatmap.2(as.matrix(x), col=heat.colors(16),Rowv=FALSE, trace="none",density.info="none",labRow=NA, margins=c(10,10), cexRow=0.8, cexCol=1.2, key=T, ColSideColors=sample, main=title)
  dev.off()
  
}
heat_tcga_scale <- function(x, plotName, title) {  
  
  library(gplots)
  #breaks for the core of the distribution
  breaks=seq(-4, 4, by=0.2) #41 values
  #now add outliers
  breaks=append(breaks, 10)
  breaks=append(breaks, -10, 0)
  #create colour panel with length(breaks)-1 colours
  mycol <- colorpanel(n=length(breaks)-1,low="green",mid="black",high="red")
  
  jpeg(plotName, width = 800, height = 800, quality=100)
  #scale
  heatmap.2(as.matrix(x), col=mycol,density.info="histogram",labCol=NA, breaks=breaks, trace="none",labRow=NA, hclustfun = function(x)hclust(x, method = "ward"), margins=c(10,10), cexRow=0.8, cexCol=1.2, key=T, ColSideColors=sample, main=title, scale="row")
  dev.off()
  
}

#Functions to create 2 venn diagrams and genelist
venn2_all <- function(x,y, l, filename) {
  
  png(paste(filename,".jpg",sep=""),width=500,height=500)
  venn_hyper <- venndiagram(x, y, unique=T, title=title, labels=l, lines=1, lcol=1, tcol=1, diacol=1, plot=T, type="2", printsub=TRUE)
  dev.off()
  ff(venn_hyper, filename)
}

#Functions to create 3 venn diagrams and genelist
venn3_all <- function(x,y,z, l, filename) {
  
  png(paste(filename,".jpg",sep=""),width=1000,height=1000)
  venn_hyper <- venndiagram(x, y, z, unique=T, title=title, labels=l, lines=1, lcol=1, tcol=1, diacol=1, plot=T, type="3", printsub=TRUE)
  dev.off()
  ff(venn_hyper, filename)
}

#Function makes genelist from venn diagram
ff <- function(venn, title) {
  Dir <- paste(getwd(),"/", sep="")
  
  
  subDir_file <- title
  #Create directory for plots and files generated from the analysis
  dir.create(file.path(Dir, subDir_file), showWarnings = FALSE)
  write(venn$q1,paste(Dir, subDir_file,"/","q1.txt", sep=""))
  write(venn$q2,paste(Dir, subDir_file,"/","q2.txt", sep=""))
  write(venn$q3,paste(Dir, subDir_file,"/","q3.txt", sep=""))
  write(venn$q4,paste(Dir, subDir_file,"/","q4.txt", sep=""))
  write(venn$q5,paste(Dir, subDir_file,"/","q5.txt", sep=""))
  write(venn$q6,paste(Dir, subDir_file,"/","q6.txt", sep=""))
  write(venn$q7,paste(Dir, subDir_file,"/","q7.txt", sep=""))
}

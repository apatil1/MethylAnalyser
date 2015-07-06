read_level3_expressions <- function(path, type) {
  
  file <- list.files(path) #gets list of files
  #write.csv(file, paste(type,".csv", sep=""))
  
  #reads data from all the files 
  for (i in 1:length(file)) {
    bar <- read.table(paste(path, file[i], sep=""), sep="", header=F, nrows=1)
    ch1 <- strsplit(as.character(bar$V3), "[-]")[[1]][4]
    ch2 <- substr(ch1, 1, 2)
    
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
  
  
  head(sample)
  assign("ss",sample, .GlobalEnv)
  
    read <- read.table(paste(path, file[i], sep=""), sep="\t", header=F, skip=2)
    
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

setwd("D:/ONC-Analysis/users/ash/michelle/")
path <- "D:/ONC-Analysis/users/ash/michelle/LUSC_expression/Expression-Genes/UNC__AgilentG4502A_07_3/Level_3/"





data <-read_level3_expressions(path, "lusc_expression")













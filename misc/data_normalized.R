data_normalized <- function(path, folder) {

  #Read files
  file<-list.files(path)

  for (i in 1:length(file)) {
    read<-read.table(paste(path, file[i], sep=""), skip=4)
  
    if (i==1) {
      data<-cbind(read[,2])
      rownames(data)<-read[,1]
      col<-strsplit(file[i],"[.]")[[1]][1]
    }
    else {
      data<-cbind(data,read[,2])
      col<-append(col,strsplit(file[i],"[.]")[[1]][1])
    }
  }
  
  colnames(data)<-col
  head(data)
  dim(data) 
  
  return (data)
}

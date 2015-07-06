sheet <- read.csv("/Users/ashwini/Desktop/michelle_data/final/sheet.csv")

path <- "/Users/ashwini/Desktop/michelle_data/files/"

sub <- "TSS_age.hyper_"
title <-"Hypermethylation"

sub <- "TSS_age.hypo_"
title <-"Hypomethylation"

i<-1
a <- read.csv(paste(path, sub ,sheet$x1[i],"-", sheet$x2[i],".csv", sep=""))
head(a)
j <- 2
b <- read.csv(paste(path,sub,sheet$x1[j],"-", sheet$x2[j],".csv", sep=""))
head(b)
k <-4
c <- read.csv(paste(path,sub,sheet$x1[k],"-", sheet$x2[k],".csv", sep=""))
head(c)

png(paste(title,".jpg",sep=""),width=1000,height=1000)
venn<-venndiagram(x=a$X, y=b$X, z=c$X, unique=T, title=title, labels=c(paste(sheet$x1[i],"-", sheet$x2[i],sep=""), paste(sheet$x1[j],"-", sheet$x2[j],sep=""), paste(sheet$x1[k],"-", sheet$x2[k],sep="")), lines=1, lcol=1, tcol=1, diacol=1, plot=T, type="3", printsub=TRUE)
dev.off()

Dir <-"/Users/ashwini/Desktop/michelle_data/final/venn_diagrams/"

#hyper
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

####################################
sub <- "TSS_age.stable.hyper_"
title <-"Stable_Hypermethylation_probes"

sub <- "TSS_age.stable.hypo_"
title <-"Stable_Hypomethylation_probes"

a <- NULL
b <- NULL
c <- NULL
d <- NULL
i<-1
a <- read.csv(paste(path, sub ,sheet$x1[i],"-", sheet$x2[i], "-", sheet$y1[i],"-", sheet$y2[i],".csv", sep=""))
head(a)
j <- 2
b <- read.csv(paste(path,sub,sheet$x1[j],"-", sheet$x2[j],"-", sheet$y1[j],"-", sheet$y2[j],".csv", sep=""))
head(b)
k <-3
c <- read.csv(paste(path,sub,sheet$x1[k],"-", sheet$x2[k],"-", sheet$y1[k],"-", sheet$y2[k],".csv", sep=""))
head(c)

#m <- 6
#d <- read.csv(paste(path,sub,sheet$x1[m],"-", sheet$x2[m],"-", sheet$y1[m],".csv", sep=""))

png(paste(title,".jpg",sep=""),width=1000,height=1000)
#venn<-venndiagram(x=a$X, y=b$X, z=c$X, w=d$X ,unique=T, title=title, labels=c(paste(sheet$x1[i],"-", sheet$x2[i],"-", sheet$y1[i],sep=""), paste(sheet$x1[j],"-", sheet$x2[j],"-", sheet$y1[j],sep=""), paste(sheet$x1[k],"-", sheet$x2[k],"-", sheet$y1[k],sep=""), paste(sheet$x1[m],"-", sheet$x2[m],"-", sheet$y1[m],sep="")), lines=1, lcol=1, tcol=1, diacol=1, plot=T, type="4", printsub=TRUE)
venn<-venndiagram(x=a$X, y=b$X, z=c$X ,unique=T, title=title, labels=c(paste(sheet$x1[i],"-", sheet$x2[i],"-", sheet$y1[i],"-", sheet$y2[i],sep=""), paste(sheet$x1[j],"-", sheet$x2[j],"-", sheet$y1[j],"-", sheet$y2[j],sep=""), paste(sheet$x1[k],"-", sheet$x2[k],"-", sheet$y1[k],"-", sheet$y2[k],sep="")), lines=1, lcol=1, tcol=1, diacol=1, plot=T, type="3", printsub=TRUE)

dev.off()

Dir <-"/Users/ashwini/Desktop/michelle_data/final/"

#hyper
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


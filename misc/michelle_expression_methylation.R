sheet <- read.csv("/Users/ashwini/Desktop/michelle_data/final/sheet.csv")

#from expression data
path1 <- "/Users/ashwini/Desktop/michelle_data/20140813 Michelle Vaz Agilent GE/up_down/"

colnames(x.agg)

re <- "up_"
re2 <- "down_"

mm <- "-15m"

i <- 8 #6m
i <- 9 #6m
i <- 10 #10m
i <- 11 #15m
up <- read.table(paste(path1, re, colnames(x.agg[i]), mm, ".txt", sep=""))
down <- read.table(paste(path1, re2, colnames(x.agg[i]), mm, ".txt", sep=""))

#########################################################
a <- NULL
b <- NULL

#from methylation data
#hyper
path2 <- "/Users/ashwini/Desktop/michelle_data/final/venn_diagrams/files/"

title <- "Stable up and down regulated genes"
sub <- "TSS_age.stable.hyper_"

j<-3
hyper <- read.csv(paste(path2, sub ,sheet$x1[j],"-", sheet$x2[j], "-", sheet$y1[j],"-", sheet$y2[j],".csv", sep=""))
head(hyper)
dim(hyper)

#hypo
sub2 <- "TSS_age.stable.hypo_"


hypo <- read.csv(paste(path2,sub2,sheet$x1[j],"-", sheet$x2[j],"-", sheet$y1[j],"-", sheet$y2[j],".csv", sep=""))
head(hypo)
dim(hypo)




l1 <- c(paste("Hyper ",sheet$x1[j],"-", sheet$x2[j], "-", sheet$y1[j],"-", sheet$y2[j], sep=""),paste(re, colnames(x.agg[i]), sep=""), paste(re2, colnames(x.agg[i]), sep=""))
l2 <- c(paste("Hypo ",sheet$x1[j],"-", sheet$x2[j], "-", sheet$y1[j],"-", sheet$y2[j], sep=""), paste(re, colnames(x.agg[i]), sep=""), paste(re2, colnames(x.agg[i]), sep=""))

#venn diagram
png(paste(title,".jpg",sep=""),width=1000,height=1000)
par(mfrow=c(1,2))
venn_hyper <- venndiagram(x=hyper$Gene, y=up$V1, z=down$V1 ,unique=T, title=title, labels=l1, lines=1, lcol=1, tcol=1, diacol=1, plot=T, type="3", printsub=TRUE)
venn_hypo <- venndiagram(x=hypo$Gene, y=up$V1, z=down$V1, unique=T, title=title, labels=l2, lines=1, lcol=1, tcol=1, diacol=1, plot=T, type="3", printsub=TRUE)

dev.off()

ff(venn_hyper, "15m-hyper-stable")
ff(venn_hypo, "15m-hypo-stable")

ff <- function(venn, title) {
Dir <-"/Users/ashwini/Desktop/michelle_data/"

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
}
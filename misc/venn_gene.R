sheet <- read.csv("/Users/ashwini/Desktop/michelle_data/final/sheet.csv")
sheet_genelist <- "/Users/ashwini/Desktop/michelle_data/final/sheet.csv"
path <- "/Users/ashwini/Desktop/michelle_data/files/files/"
folder <- "/Users/ashwini/Desktop/michelle_data/files/files/"

x=which(abs(beta.tab_CpGI.PlusMinus1500_TSS$unt_6 -  beta.tab_CpGI.PlusMinus1500_TSS$unt_1) <= 0.05 & abs(beta.tab_CpGI.PlusMinus1500_TSS$unt_10 -  beta.tab_CpGI.PlusMinus1500_TSS$unt_1) <= 0.05 & abs(beta.tab_CpGI.PlusMinus1500_TSS$unt_15 -  beta.tab_CpGI.PlusMinus1500_TSS$unt_1) <= 0.05)
universal_age_stable <- beta.tab_CpGI.PlusMinus1500_TSS[x, ]

x=which(abs(beta.tab_CpGI.PlusMinus1500_TSS$unt_6 -  beta.tab_CpGI.PlusMinus1500_TSS$unt_1) >= 0.05  & abs(beta.tab_CpGI.PlusMinus1500_TSS$unt_6 -  beta.tab_CpGI.PlusMinus1500_TSS$unt_1) <= 0.2  & abs(beta.tab_CpGI.PlusMinus1500_TSS$unt_10 -  beta.tab_CpGI.PlusMinus1500_TSS$unt_1) >= 0.05  & abs(beta.tab_CpGI.PlusMinus1500_TSS$unt_10 -  beta.tab_CpGI.PlusMinus1500_TSS$unt_1) <= 0.2 & abs(beta.tab_CpGI.PlusMinus1500_TSS$unt_15 -  beta.tab_CpGI.PlusMinus1500_TSS$unt_1) >= 0.05 & abs(beta.tab_CpGI.PlusMinus1500_TSS$unt_15 -  beta.tab_CpGI.PlusMinus1500_TSS$unt_1) <= 0.2)
universal_age_inter <- beta.tab_CpGI.PlusMinus1500_TSS[x, ]

analysis_methylation4(universal_age_inter, sheet_genelist, folder, "TSS") 
analysis_methylation4(beta.tab_CpGI.PlusMinus1500_TSS, sheet_genelist, folder, "TSS_all") 
delta_plot(universal_age_inter, sheet_genelist, folder, "TSS_universal") 

sub <- "TSS_age.inter.hyper_"
title <-"inter_Hypermethylation_gene"

sub <- "TSS_age.inter.hypo_"
title <-"inter_Hypomethylation_gene"

a <- NULL
b <- NULL
c <- NULL
d <- NULL
i<-1
a <- read.csv(paste(path, sub ,sheet$x1[i],"-", sheet$x2[i], "-", sheet$y1[i],"-", sheet$y2[i],".csv", sep=""))
head(a)
dim(a)
j <- 2
b <- read.csv(paste(path,sub,sheet$x1[j],"-", sheet$x2[j],"-", sheet$y1[j],"-", sheet$y2[j],".csv", sep=""))
head(b)
dim(b)
k <-3
c <- read.csv(paste(path,sub,sheet$x1[k],"-", sheet$x2[k],"-", sheet$y1[k],"-", sheet$y2[k],".csv", sep=""))
head(c)
dim(c)
#m <- 6
#d <- read.csv(paste(path,sub,sheet$x1[m],"-", sheet$x2[m],"-", sheet$y1[m],".csv", sep=""))

png(paste(title,".jpg",sep=""),width=1000,height=1000)
#venn<-venndiagram(x=a$Gene, y=b$Gene, z=c$Gene, w=d$Gene ,unique=T, title=title, labels=c(paste(sheet$x1[i],"-", sheet$x2[i],"-", sheet$y1[i],sep=""), paste(sheet$x1[j],"-", sheet$x2[j],"-", sheet$y1[j],sep=""), paste(sheet$x1[k],"-", sheet$x2[k],"-", sheet$y1[k],sep=""), paste(sheet$x1[m],"-", sheet$x2[m],"-", sheet$y1[m],sep="")), lines=1, lcol=1, tcol=1, diacol=1, plot=T, type="4", printsub=TRUE)
venn<-venndiagram(x=a$Gene, y=b$Gene, z=c$Gene ,unique=T, title=title, labels=c(paste(sheet$x1[i],"-", sheet$x2[i],"-", sheet$y1[i],"-", sheet$y2[i],sep=""), paste(sheet$x1[j],"-", sheet$x2[j],"-", sheet$y1[j],"-", sheet$y2[j],sep=""), paste(sheet$x1[k],"-", sheet$x2[k],"-", sheet$y1[k],"-", sheet$y2[k],sep="")), lines=1, lcol=1, tcol=1, diacol=1, plot=T, type="3", printsub=TRUE)

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

ven <- venn$q1
type <- "intermediate Hypermethylation gene q1"
plot(ven, type) 

ven <- venn$q2
type <- "intermediate Hypermethylation gene q2"
plot(ven, type) 

ven <- venn$q3
type <- "intermediate Hypermethylation gene q3"
plot(ven, type) 

ven <- venn$q4
type <- "intermediate Hypermethylation gene q4"
plot(ven, type) 

ven <- venn$q5
type <- "intermediate Hypermethylation gene q5"
plot(ven, type) 

ven <- venn$q6
type <- "intermediate Hypermethylation gene q6"
plot(ven, type) 

ven <- venn$q7
type <- "intermediate Hypermethylation gene q7"
plot(ven, type) 

###############
ven <- venn$q1
type <- "intermediate Hypomethylation gene q1"
plot(ven, type) 

ven <- venn$q2
type <- "intermediate Hypomethylation gene q2"
plot(ven, type) 

ven <- venn$q3
type <- "intermediate Hypomethylation gene q3"
plot(ven, type) 

ven <- venn$q4
type <- "intermediate Hypomethylation gene q4"
plot(ven, type) 

ven <- venn$q5
type <- "intermediate Hypomethylation gene q5"
plot(ven, type) 

ven <- venn$q6
type <- "intermediate Hypomethylationgene q6"
plot(ven, type) 

ven <- venn$q7
type <- "intermediate Hypomethylation gene q7"
plot(ven, type) 

plot <- function(ven, type) {
hy<- x.agg[x.agg$GeneName %in% ven, ]
head(hy)
dim(hy)
png(paste(type,".jpg", sep=""), width= 1000, height=1000)
boxplot(hy[,2:11], main=paste("Expression data ", type, sep=""), xlab="Samples", ylab="M values")
dev.off()
}




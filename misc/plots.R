#Load library
library(IlluminaHumanMethylation450kmanifest)
setwd("~/Desktop/new-annotation")
#Load image 
load("~/Desktop/methylation_analysis/R_image/Michelle-1500TSS.RData")
load("~/Desktop/methylation_analysis/R_image/Michelle-all-CPGI.RData")

########################################################################################################################
#Scatter and density plots
########################################################################################################################

#Type1 and Type2 probes plots
IlluminaHumanMethylation450kmanifest
type1<-getProbeInfo(IlluminaHumanMethylation450kmanifest, type = "I")
type2<-getProbeInfo(IlluminaHumanMethylation450kmanifest, type = "II")
type1_cpg<-type1$Name
type2_cpg<-type2$Name

#data type1/2 consist of beta.tab_CpGI.PlusMinus1500_TSS seperated by type1 and type2 probes
data_type1 <- data.frame(beta.tab_CpGI.PlusMinus1500_TSS[rownames(beta.tab_CpGI.PlusMinus1500_TSS) %in% type1_cpg,])
data_type2 <- data.frame(beta.tab_CpGI.PlusMinus1500_TSS[rownames(beta.tab_CpGI.PlusMinus1500_TSS) %in% type2_cpg,])

#+/- 1500 TSS CPGI
colnames(data_type1)


m <- c(2,3,5,8,8, 3,3,3,3,3)
n <- c(1,4,6,7,9, 1,4,6,7,9)

for (i in 1:10) {
  x <- m[i]
  y <- n[i]
  plt(x, y, paste("scatter", "_", colnames(data_type1[x]), "-", colnames(data_type1[y]), ".jpg", sep=""))
}

#plots
plt <- function(x,y,name) {

  png(name, width = 1000, height = 1000)

  plot(data_type1[, x], data_type1[, y], cex.axis =1.6, col=NULL, bg=rgb(0,1,0, alpha=0.5, maxColorValue = 1), pch=21, cex=0.4, xlab=colnames(data_type1[x]), ylab=colnames(data_type1[y]),cex.lab=1.5, cex.axis=1.8)
  points(data_type2[, x], data_type2[, y], col=NULL, bg=rgb(0,0,1, alpha=0.5, maxColorValue = 1), pch=21, cex=0.4)

  legend("topright", inset=c(-0,0), legend=c("Type 1 probes", "Type 2 probes"), col=c("green", "blue"), pt.bg=c("green", "blue"), pch=21, ncol=1, cex=1.5, pt.cex = 1) # inset might need adjusting based on the width of the legend
  title(main="ScatterPlot CpGI +/- 1500 TSS", cex.main=1, col="brown")

  dev.off()
}

#density T1 and T2 all cpgi +/- 1500 TSS probes
colnames(data_type1)

for (i in 1:10) {
  x <- m[i]
  y <- n[i]
  den(x, y, paste("Density_cpgi_1500_tss", "_", colnames(data_type1[x]), "-", colnames(data_type1[y]), ".jpg", sep=""))
}




den <- function(x, y, name) {
  
  png(paste("T1_", name, sep=""), width = 1000, height = 1000)
  
  #density T1 treated, untreated
  plot(density(data_type1[, y]),cex.axis =1.6,col="blue",main="Type 1 +/- 1500 TSS CpG island",cex.lab=1.5)
  lines(density(data_type1[, x]),col="red")
  legend("topright", inset=c(-0,0), legend=c(colnames(data_type1[y]), colnames(data_type1[x])), col=c("blue", "red"), pt.bg=c("blue", "red"), pch=21, ncol=1, cex=1.5, pt.cex = 1)
  dev.off()
  
  #density T2 treated, untreated
  png(paste("T2_", name, sep=""), width = 1000, height = 1000)
  plot(density(data_type2[, y]),col="blue",cex.axis =1.6,main="Type 2 +/- 1500 TSS island",cex.lab=1.5)
  lines(density(data_type2[, x]),col="red")
  legend("topright", inset=c(-0,0), legend=c(colnames(data_type1[y]), colnames(data_type1[x])), col=c("blue", "red"), pt.bg=c("blue", "red"), pch=21, ncol=1, cex=1.5, pt.cex = 1)
  dev.off()
}

#for all cpg island probes
data_type1 <- data.frame(beta.tab_CpGI[rownames(beta.tab_CpGI) %in% type1_cpg,])
data_type2 <- data.frame(beta.tab_CpGI[rownames(beta.tab_CpGI) %in% type2_cpg,])

colnames(data_type1)

for (i in 1:10) {
  x <- m[i]
  y <- n[i]
  plt(x, y, paste("scatter_all_cpgi", "_", colnames(data_type1[x]), "-", colnames(data_type1[y]), ".jpg", sep=""))
}


#plots
plt <- function(x,y,name) {
  
  png(name, width = 1000, height = 1000)
  
  plot(data_type1[, x], data_type1[, y], cex.axis =1.6,col=NULL, bg=rgb(0,1,0, alpha=0.5, maxColorValue = 1), pch=21, cex=0.4, xlab=colnames(data_type1[x]), ylab=colnames(data_type1[y]),cex.lab=1.5)
  points(data_type2[, x], data_type2[, y], col=NULL, bg=rgb(0,0,1, alpha=0.5, maxColorValue = 1), pch=21, cex=0.4)
  
  legend("topright", inset=c(-0,0), legend=c("Type 1 probes", "Type 2 probes"), col=c("green", "blue"), pt.bg=c("green", "blue"), pch=21, ncol=1, cex=1.5, pt.cex = 1) # inset might need adjusting based on the width of the legend
  title(main="ScatterPlot All CpG island probes", cex.main=1, col="brown")
  
  dev.off()
}

#density T1 and T2 all cpgi probes
colnames(data_type1)

for (i in 1:10) {
  x <- m[i]
  y <- n[i]
  den(x, y, paste("Density_all_cpgi", "_", colnames(data_type1[x]), "-", colnames(data_type1[y]), ".jpg", sep=""))
}



den <- function(x, y, name) {

  png(paste("T1_", name, sep=""), width = 1000, height = 1000)

  #density T1 treated, untreated
  plot(density(data_type1[, y]),cex.axis =1.6,col="blue",main="Type 1 all CpG island", cex.lab=1.5)
  lines(density(data_type1[, x]),col="red")
  legend("topright", inset=c(-0,0), legend=c(colnames(data_type1[y]), colnames(data_type1[x])), col=c("blue", "red"), pt.bg=c("blue", "red"), pch=21, ncol=1, cex=1.5, pt.cex = 1)
  dev.off()

  #density T2 treated, untreated
  png(paste("T2_", name, sep=""), width = 1000, height = 1000)
  plot(density(data_type2[, y]),cex.axis =1.6,col="blue",main="Type 2 all CpG island", cex.lab=1.5)
  lines(density(data_type2[, x]),col="red")
  legend("topright", inset=c(-0,0), legend=c(colnames(data_type1[y]), colnames(data_type1[x])), col=c("blue", "red"), pt.bg=c("blue", "red"), pch=21, ncol=1, cex=1.5, pt.cex = 1)
  dev.off()
}

#scatter plot 
for (i in 1:10) {
  x <- m[i]
  y <- n[i]
  scatter(x, y, paste("scatter_probes", "_", colnames(data_type1[x]), "-", colnames(data_type1[y]), ".jpg", sep=""))
}


scatter <- function(x, y, name) {

  png(name, width = 1000, height = 1000)
  
  plot(beta.tab_CpGI[, y], beta.tab_CpGI[, x], cex.axis =1.6,col=NULL, bg=rgb(0,1,0, alpha=0.5, maxColorValue = 1), pch=21, cex=0.4, , xlab=colnames(data_type1[x]), ylab=colnames(data_type1[y]),cex.lab=1.5)
  points(beta.tab_CpGI.PlusMinus1500_TSS[, y], beta.tab_CpGI.PlusMinus1500_TSS[, x], col=NULL, bg=rgb(0,0,1, alpha=0.5, maxColorValue = 1), pch=21, cex=0.4)

  legend("topright", inset=c(0,0), legend=c("All_CpGI", "CpGI +/-1500bp TSS"), col=c("green", "blue"), pt.bg=c("green", "blue"), pch=21, ncol=1, cex=1.5, pt.cex = 1) # inset might need adjusting based on the width of the legend
  title(main="ScatterPlot", cex.main=1, col="brown")

dev.off()
}

#shore
#scatter plot 
for (i in 1:10) {
  x <- m[i]
  y <- n[i]
  scatter(x, y, paste("scatter_probes_shore", "_", colnames(data_type1[x]), "-", colnames(data_type1[y]), ".jpg", sep=""))
}


scatter <- function(x, y, name) {
  
  png(name, width = 1000, height = 1000)
  
  plot(beta.tab_CpGI[, y], beta.tab_CpGI[, x], cex.axis =1.6,col=NULL, bg=rgb(0,1,0, alpha=0.5, maxColorValue = 1), pch=21, cex=0.4, , xlab=colnames(data_type1[x]), ylab=colnames(data_type1[y]),cex.lab=1.5)
  points(beta.tab_CpGI.PlusMinus1500_TSS[, y], beta.tab_CpGI.PlusMinus1500_TSS[, x], col=NULL, bg=rgb(0,0,1, alpha=0.5, maxColorValue = 1), pch=21, cex=0.4)
  
  legend("topright", inset=c(0,0), legend=c("All_CpGI", "shore"), col=c("green", "blue"), pt.bg=c("green", "blue"), pch=21, ncol=1, cex=1.5, pt.cex = 1) # inset might need adjusting based on the width of the legend
  title(main="ScatterPlot", cex.main=1, col="brown")
  
  dev.off()
}


head(beta.tab_CpGI.PlusMinus1500_TSS)
dim(beta.tab_CpGI.PlusMinus1500_TSS)

#pairs plot
png("pairs_plot.jpg", width = 1000, height = 1000)
pairs(beta.tab_CpGI.PlusMinus1500_TSS,pch='.', lower.panel=NULL)
dev.off()

#pairs plot
png("pairs_plot_shore.jpg", width = 1000, height = 1000)
pairs(shore,pch='.', lower.panel=NULL)
dev.off()


#pairs plot
png("pairs_plot_shelf.jpg", width = 1000, height = 1000)
pairs(shelf,pch='.', lower.panel=NULL)
dev.off()

#pairs plot
png("pairs_plot_body.jpg", width = 1000, height = 1000)
pairs(body,pch='.', lower.panel=NULL)
dev.off()


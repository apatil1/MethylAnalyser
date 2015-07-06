#load("")

head(beta.tab_CpGI.PlusMinus1500_TSS)
dim(beta.tab_CpGI.PlusMinus1500_TSS)

# Create a vector of different colrs used for teh plot functions [rang is colors ;-) ]
set.seed(123)
num.colors <- ncol(beta.tab_CpGI.PlusMinus1500_TSS)
temp.rang.rgb <- rgb(red=sample(seq(0,1,length.out=num.colors)), green=sample(seq(0,1,length.out=num.colors)) , blue=sample(seq(0,1,length.out=num.colors)), maxColorValue=1)


# PLOT PCA PLOTS
PCA <- prcomp(t(beta.tab_CpGI.PlusMinus1500_TSS))

pdf("pca.pdf", width = 8, height = 8)

# Add extra space to right of plot area; change clipping to figure
par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)
plot(PCA$x[, c("PC1", "PC2")], xlab=paste("PC1: ",signif((((PCA$sdev)^2)[1]/(sum((PCA$sdev)^2))),2),sep=""), ylab=paste("PC2: ",signif((((PCA$sdev)^2)[2]/(sum((PCA$sdev)^2))),2),sep=""), pch=19, cex=1, bty='L')
for(u in 1:nrow(PCA$x)){
  points(matrix(PCA$x[u,c("PC1", "PC2")], ncol=2), pch=19, col=temp.rang.rgb[u], cex=1)
}

# Add legend to top right, outside plot region
legend("topright", inset=c(-0.3,0), legend=rownames(PCA$x), col=temp.rang.rgb, pch=19, ncol=2, cex=0.7) # inset might need adjusting based on the width of the legend
title(main="PCA plot (gives a snapshot of the variability of samples by\nreducing the variability in the data to two dimensions-\ncloser the points, lesser the vaiability)", sub=i, cex.main=0.7, col.sub="brown")

dev.off()

# PLOT HEATMAP OF VARIABLE PROBES
## Plot density of sd distribution to choose the sd filter for choosing most variable probes for heatmap
beta.tab_CpGI.PlusMinus1500_TSS.sd <- apply(beta.tab_CpGI.PlusMinus1500_TSS, 1, sd)

pdf("density_sd.pdf", width = 8, height = 8)
plot(density(beta.tab_CpGI.PlusMinus1500_TSS.sd), main="")
title(main="Density plot of standard deviations of beta values", cex.main=0.8)
dev.off()

library(heatmap.plus)
heatmap.sd.filter<-0.2
beta.tab_CpGI.PlusMinus1500_TSS.sd.filtered <- beta.tab_CpGI.PlusMinus1500_TSS[which(beta.tab_CpGI.PlusMinus1500_TSS.sd >= heatmap.sd.filter), ]
dim(beta.tab_CpGI.PlusMinus1500_TSS); dim(beta.tab_CpGI.PlusMinus1500_TSS.sd.filtered)
colsidecol <- vector()
for (u in 1:ncol(beta.tab_CpGI.PlusMinus1500_TSS)) {colsidecol[u] <- as.character(u)} #ColSideColors requires colors as characters
#plot.name <- paste(i, "pdf", sep=".")
#plot.name <- paste("Heatmap", plot.name, sep="-")
#plot.name <- paste(plot.dir, plot.name, sep="/")
pdf("heatmap.pdf", width = 8, height = 8)
heatmap(as.matrix(beta.tab_CpGI.PlusMinus1500_TSS.sd.filtered), col=heat.colors(16), trace="none", density.info="none", hclustfun = function(x)hclust(x, method = "ward.D"), margins=c(10,10), cexRow=0.8, cexCol=1.2, key=T, ColSideColors=colsidecol)
dev.off()

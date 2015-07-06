
#Function to read idat files using minfi
data_methylation <- function(baseDir, folder) {
  
  #Minfi QC and normalization with infinum background subtraction and control normalization and swan normalization
  #Read data from idat files
  k<-list.files(baseDir)
  targets <- read.450k.sheet(baseDir)
  sub(baseDir, "", targets$Basename)
  RGset <- read.450k.exp(base = baseDir, targets = targets)
  pd <- pData(RGset)
  
  #Filter probes with minfi probe detection p-value > 0.01
  infi.SOI.filtered.beta <- getBeta(RGset, type = "Illumina")
  total <- length(rownames(infi.SOI.filtered.beta))
  detP <- detectionP(RGset)
  failed <- detP > 0.01
  t.pvalue.col <- grep("TRUE", failed) 
  infi.SOI.filtered.beta[t.pvalue.col] <- NA
  filter_probes_pval <- rownames(infi.SOI.filtered.beta[complete.cases(infi.SOI.filtered.beta),])
  
  #Generate a QC report from minfi package
  qcReport(RGset, sampNames = pd$Sample_Name, sampGroups = pd$Sample_Group, pdf = paste(folder,"qcReport.pdf", sep=""))
  
  #normalization illumina
  MSet.norm <- preprocessIllumina(RGset, bg.correct = TRUE, normalize = "controls", reference = 2)
  
  #swan normalization with Illumina
  Mset.swan.norm <- preprocessSWAN(RGset, MSet.norm)
  
  #Remove probes with detection pval>0.01
  infi.SOI.filtered.beta<-data.frame(getBeta(Mset.swan.norm, type = "Illumina"))
  colnames(infi.SOI.filtered.beta) <- pd$Sample_Group
  infi.SOI.filtered.beta<-infi.SOI.filtered.beta[rownames(infi.SOI.filtered.beta) %in% filter_probes_pval,]
  
  #Remove probes from X and Y Chromosome
  infi.SOI.filtered.beta<-infi.SOI.filtered.beta[rownames(infi.SOI.filtered.beta) %in% all.probes.xy.filter,]
  
  
  
  summary <- paste("Total number of probes :", total, "\nTotal number of probes after filtering minfi probe detection p-value > 0.01 :", length(rownames(infi.SOI.filtered.beta)), "\nNo. of probes removed : ", total - length(rownames(infi.SOI.filtered.beta)), sep="")
  
  #write summary table to file
  write.table(summary, paste(folder,"summary_data.txt",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
  
  return (infi.SOI.filtered.beta)
}


#Function for analysis of methylation data and to generate genelist and delta plots
analysis_methylation <- function(tss, sheet, Dir, type) {
  
  
  #Read the sheet consisting of questions
  dat <- read.csv(sheet, header=T)
  
  summary <- NULL #For summary of the analysis
  rname <- NULL #For rownames of summary of the analysis
  subDir_file <- "files"
  
  #Create directory for plots and files generated from the analysis
  dir.create(file.path(Dir, subDir_file), showWarnings = FALSE)
  
  for(i in 1:length(dat$first)) {
    
    #Get CpG island probes within +/-1500 bp from TSS that get methylated/demethylated by more than 20% (age/proliferation dependent methylation), and then the probes that further get methylated/demethyalted upon CSC exposure.
    age.hypermethylation <- tss[which(tss[colnames(tss)==dat$first[i]] - tss[colnames(tss)==dat$second[i]] >= 0.2), ]
    age.hypomethylation <- tss[which(tss[colnames(tss)==dat$first[i]] - tss[colnames(tss)==dat$second[i]] <= -0.2), ]
    
    #Now, from the probes that get hypermethylated due to ageing in culture, get the probes that get further methylated upon further exposure. And, from the probes that get hypomethylated due to ageing in culture, get the probes that get further demethylated upon further exposure
    age.fur.hypermethylation <- age.hypermethylation[which(age.hypermethylation[,colnames(age.hypermethylation)== dat$further[i]] - age.hypermethylation[,colnames(age.hypermethylation)== dat$first[i]] >= 0.2), ]
    age.fur.hypomethylation <- age.hypomethylation[which(age.hypomethylation[,colnames(age.hypomethylation)== dat$further[i]] - age.hypomethylation[,colnames(age.hypomethylation)== dat$first[i]] <= -0.2), ]
    
    #Get CpG island probes that get methylated specifically upon fur exposure 
    #abs() <= 0.05 will pick probes that do not undergo much changes in methylation upon culturing/ageing.
    age.Stable_fur.hypmeth.probes <- tss[abs(tss[,colnames(tss)==dat$first[i]] - tss[,colnames(tss)==dat$second[i]])  <= 0.05 & tss[,colnames(tss)==dat$further[i]] - tss[,colnames(tss)==dat$first[i]]  >= 0.2, ]
    age.Stable_fur.hypometh.probes <- tss[abs(tss[,colnames(tss)==dat$first[i]] - tss[,colnames(tss)==dat$second[i]]) <= 0.05 & tss[,colnames(tss)==dat$further[i]] - tss[,colnames(tss)==dat$first[i]] <= -0.2, ]
    
    #Get probes that range in methylation between +/-0.05 to +/-0.2 during aging/proliferation but are methylated/demethylated >= 0.2 
    age.intermediate_meth.fur.hypmeth.probes <- tss[abs(tss[,colnames(tss)==dat$first[i]] - tss[,colnames(tss)==dat$second[i]]) >= 0.05 & abs(tss[,colnames(tss)==dat$first[i]] - tss[,colnames(tss)==dat$second[i]]) <= 0.2 & tss[,colnames(tss)==dat$further[i]] - tss[,colnames(tss)==dat$first[i]] >= 0.2, ]
    age.intermediate_meth.fur.hypometh.probes <- tss[abs(tss[,colnames(tss)==dat$first[i]] - tss[,colnames(tss)==dat$second[i]]) >= 0.05 & abs(tss[,colnames(tss)==dat$first[i]] - tss[,colnames(tss)==dat$second[i]]) <= 0.2 & tss[,colnames(tss)==dat$further[i]] - tss[,colnames(tss)==dat$first[i]] <= -0.2, ]
    
    
    #Get gene names of the probes 
    if (length(rownames(age.hypermethylation)) > 0) {
      age.hypermethylation <- cbind.data.frame(age.hypermethylation, "Gene"=rep(NA, times=nrow(age.hypermethylation)))
      for(j in 1:nrow(age.hypermethylation)) {
        probe.id <- rownames(age.hypermethylation)[j]
        gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
        if(length(gene) == 1) { age.hypermethylation$Gene[j] <- as.character(gene) }
      }
      write.csv(age.hypermethylation, paste(Dir, subDir_file,"/", type, "_age.hyper","_",dat$first[i],"-",dat$second[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }
    
    if (length(rownames(age.hypomethylation)) > 0) {
      age.hypomethylation <- cbind.data.frame(age.hypomethylation, "Gene"=rep(NA, times=nrow(age.hypomethylation)))
      for(j in 1:nrow(age.hypomethylation)) {
        probe.id <- rownames(age.hypomethylation)[j]
        gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
        if(length(gene) == 1) { age.hypomethylation$Gene[j] <- as.character(gene) }
      }
      write.csv(age.hypomethylation, paste(Dir, subDir_file,"/", type, "_age.hypo","_",dat$first[i],"-",dat$second[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }
    
    if (length(rownames(age.fur.hypermethylation)) > 0) {
      age.fur.hypermethylation <- cbind.data.frame(age.fur.hypermethylation, "Gene"=rep(NA, times=nrow(age.fur.hypermethylation)))
      for(j in 1:nrow(age.fur.hypermethylation)) {
        probe.id <- rownames(age.fur.hypermethylation)[j]
        gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
        if(length(gene) == 1) { age.fur.hypermethylation$Gene[j] <- as.character(gene) }
      }
      write.csv(age.fur.hypermethylation, paste(Dir, subDir_file,"/", type, "_age.fur.hyper","_",dat$first[i],"-",dat$second[i],"-",dat$further[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }
    
    if (length(rownames(age.fur.hypomethylation)) > 0) {
      age.fur.hypomethylation <- cbind.data.frame(age.fur.hypomethylation, "Gene"=rep(NA, times=nrow(age.fur.hypomethylation)))
      for(j in 1:nrow(age.fur.hypomethylation)) {
        probe.id <- rownames(age.fur.hypomethylation)[j]
        gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
        if(length(gene) == 1) { age.fur.hypomethylation$Gene[j] <- as.character(gene) }
      }
      write.csv(age.fur.hypomethylation, paste(Dir, subDir_file,"/", type, "_age.fur.hypo","_",dat$first[i],"-",dat$second[i],"-",dat$further[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }
    
    if (length(rownames(age.Stable_fur.hypmeth.probes)) > 0) {
      age.Stable_fur.hypmeth.probes <- cbind.data.frame(age.Stable_fur.hypmeth.probes, "Gene"=rep(NA, times=nrow(age.Stable_fur.hypmeth.probes)))
      for(j in 1:nrow(age.Stable_fur.hypmeth.probes)) {
        probe.id <- rownames(age.Stable_fur.hypmeth.probes)[j]
        gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
        if(length(gene) == 1) { age.Stable_fur.hypmeth.probes$Gene[j] <- as.character(gene) }
      }
      write.csv(age.Stable_fur.hypmeth.probes, paste(Dir, subDir_file,"/", type, "_age.stable.hyper","_",dat$first[i],"-",dat$second[i],"-",dat$further[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }
    
    if (length(rownames(age.Stable_fur.hypometh.probes)) > 0) {
      age.Stable_fur.hypometh.probes <- cbind.data.frame(age.Stable_fur.hypometh.probes, "Gene"=rep(NA, times=nrow(age.Stable_fur.hypometh.probes)))
      for(j in 1:nrow(age.Stable_fur.hypometh.probes)) {
        probe.id <- rownames(age.Stable_fur.hypometh.probes)[j]
        gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
        if(length(gene) == 1) { age.Stable_fur.hypometh.probes$Gene[j] <- as.character(gene) }
      }
      write.csv(age.Stable_fur.hypometh.probes, paste(Dir, subDir_file,"/", type, "_age.stable.hypo","_",dat$first[i],"-",dat$second[i],"-",dat$further[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }
    
    if (length(rownames(age.intermediate_meth.fur.hypmeth.probes)) > 0) {
      age.intermediate_meth.fur.hypmeth.probes <- cbind.data.frame(age.intermediate_meth.fur.hypmeth.probes, "Gene"=rep(NA, times=nrow(age.intermediate_meth.fur.hypmeth.probes)))
      for(j in 1:nrow(age.intermediate_meth.fur.hypmeth.probes)) {
        probe.id <- rownames(age.intermediate_meth.fur.hypmeth.probes)[j]
        gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
        if(length(gene) == 1) { age.intermediate_meth.fur.hypmeth.probes$Gene[j] <- as.character(gene) }
      }
      write.csv(age.intermediate_meth.fur.hypmeth.probes, paste(Dir, subDir_file,"/", type, "_age.inter.hyper","_",dat$first[i],"-",dat$second[i],"-",dat$further[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }
    
    if (length(rownames(age.intermediate_meth.fur.hypometh.probes)) > 0) {
      age.intermediate_meth.fur.hypometh.probes <- cbind.data.frame(age.intermediate_meth.fur.hypometh.probes, "Gene"=rep(NA, times=nrow(age.intermediate_meth.fur.hypometh.probes)))
      for(j in 1:nrow(age.intermediate_meth.fur.hypometh.probes)) {
        probe.id <- rownames(age.intermediate_meth.fur.hypometh.probes)[j]
        gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
        if(length(gene) == 1) { age.intermediate_meth.fur.hypometh.probes$Gene[j] <- as.character(gene) }
      }
      write.csv(age.intermediate_meth.fur.hypometh.probes, paste(Dir, subDir_file,"/", type, "_age.inter.hypo","_",dat$first[i],"-",dat$second[i],"-",dat$further[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }
    
    #Add rows to summary file
    summary<-rbind(summary,cbind(length(rownames(age.hypermethylation)), length(rownames(age.hypomethylation)),length(rownames(age.fur.hypermethylation)), length(rownames(age.fur.hypomethylation)), length(rownames(age.Stable_fur.hypmeth.probes)), length(rownames(age.Stable_fur.hypometh.probes)), length(rownames(age.intermediate_meth.fur.hypmeth.probes)), length(rownames(age.intermediate_meth.fur.hypometh.probes))))
    rname<-cbind(rname,paste(dat$first[i] ,"-", dat$second[i], "-", dat$further[i]))
    
  }
  
  #Add colnames and rownames to summary
  colnames(summary)<-cbind("age.hypermethylation", "age.hypomethylation", "age.fur.hypermethylation", "age.fur.hypomethylation", "age.Stable_fur.hypmeth.probes", "age.Stable_fur.hypometh.probes", "age.intermediate_meth.fur.hypmeth.probes","age.intermediate_meth.fur.hypometh.probes")
  rownames(summary)<-rname
  
  #write summary table to file
  write.csv(summary, paste(Dir, type, "_", "summary_probes.csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
  
}


#Function for plotting delta beta plots
delta_plot <- function(tss, sheet, Dir, type) {
  
  #Read the sheet consisting of questions
  dat <- read.csv(sheet, header=T)
  
  subDir_plot <- "Delta-Plots"
  
  #Create directory for plots and files generated from the analysis
  dir.create(file.path(Dir, subDir_plot), showWarnings = FALSE)
  
  for(i in 1:length(dat$x1)) {
    #Delta plot 
    
    #Calculate Dela Beta
    x <- tss[,colnames(tss)==dat$x1[i]] - tss[,colnames(tss)==dat$x2[i]]
    y <- tss[,colnames(tss)==dat$y1[i]] - tss[,colnames(tss)==dat$y2[i]]
    
    #Calculate Delta Beta for probes in the control having Beta values > 0.5
    red <- tss[tss[,colnames(tss)==dat$x2[i]] > 0.5, ]
    point.red.x <- red[,colnames(red)==dat$x1[i]] - red[,colnames(red)==dat$x2[i]]
    point.red.y <- red[,colnames(red)==dat$y1[i]] - red[,colnames(red)==dat$y2[i]]
    
    #Calculate Delta Beta for probes in the control having Beta values < 0.2
    blue <- tss[tss[,colnames(tss)==dat$x2[i]] < 0.2, ]
    point.blue.x <- blue[,colnames(blue)==dat$x1[i]] - blue[,colnames(blue)==dat$x2[i]]
    point.blue.y <- blue[,colnames(blue)==dat$y1[i]] - blue[,colnames(blue)==dat$y2[i]]
    
    png(paste(Dir, subDir_plot,"/", type, "_",dat$x1[i],"-",dat$x2[i],"-",dat$y1[i], "-",dat$y2[i],".jpg",sep=""), width = 1000, height = 1000)
    
    #plot the Delta Beta values
    plot(y~x, cex=0.4, cex.axis =1.8 ,cex.lab=1.5,ylim=range(-1,1),xlim=range(-1,1), col=rgb(0.6,0.6,0.6, alpha=0.4), pch=19, xlab=paste(dat$x1[i],"-",dat$x2[i],sep=""), ylab=paste(dat$y1[i],"-",dat$y2[i],sep=""))
    points(point.red.x, point.red.y, col=NULL, bg=rgb(0,0,0.702, alpha=0.4, maxColorValue = 1), pch=21, cex=0.4)
    points(point.blue.x, point.blue.y, col=NULL, bg=rgb(1,0,0, alpha=0.4, maxColorValue = 1), pch=21, cex=0.4)
    
    #Add lines to x and y axis
    abline(v=-0.05, col="green"); abline(v=0.05, col="green"); abline(v=-0.2, col="pink"); abline(v=0.2, col="pink")
    abline(h=0.2, col="pink"); abline(h=-0.2, col="pink")
    
    #Add title to the plot
    title(main=paste( type, "_", dat$x1[i],"-",dat$x2[i]," vs ",dat$y1[i],"-",dat$y2[i],sep=""), cex.main=1, col="brown")
    
    #Add legend to the plot
    legend("topright", inset=c(0,0), legend=c("All probes", paste(dat$x2[i]," > 0.5",sep=""), paste(dat$x2[i]," < 0.2",sep="")), col=c("grey","darkblue","red"), pt.bg=c("grey","darkblue","red"), pch=21, ncol=1, cex=2) # inset might need adjusting based on the width of the legend
    
    #Add counts of probes in each block in the plot
    text(-0.3,0.6,length(x[(x< -0.2) & (y > 0.2)]),cex=1.5)
    text(-0.3,0,length(x[(x< -0.2) & ((y> -0.2) & (y< 0.2))]),cex=1.5)
    text(-0.3,-0.6,length(x[(x < -0.2)  & (y < -0.2)]),cex=1.5)
    text(-0.1,-0.6,length(x[(x < -0.05) & (x > -0.2)  & (y < -0.2)]),cex=1.5)
    text(-0.1,0.6,length(x[(x < -0.05) & (x > -0.2)  & (y > 0.2)]),cex=1.5)
    text(-0.1,0,length(x[(x < -0.05) & (x > -0.2) & (y < 0.2) & (y > -0.2)]),cex=1.5)
    text(-0,0,length(x[(x < 0.05)  & (x > -0.05) & (y < 0.2) & (y > -0.2)]),cex=1.5)
    text(-0,0.6,length(x[(x < 0.05)  & (x > -0.05) & (y > 0.2)]),cex=1.5)
    text(-0,-0.6,length(x[(x < 0.05)  & (x > -0.05) & (y < -0.2)]),cex=1.5)
    text(0.1,-0.6,length(x[(x < 0.2)  & (x > 0.05) & (y < -0.2)]),cex=1.5)
    text(0.1,0.6,length(x[(x < 0.2)  & (x > 0.05) & (y > 0.2)]),cex=1.5)
    text(0.1,0,length(x[(x < 0.2)  & (x > 0.05) & (y < 0.2) & (y > -0.2)]),cex=1.5)
    text(0.3,0,length(x[(x > 0.2) & (y < 0.2) & (y > -0.2)]),cex=1.5)
    text(0.3,0.6,length(x[(x > 0.2) & (y > 0.2)]),cex=1.5)
    text(0.3,-0.6,length(x[(x > 0.2) & (y < -0.2) ]),cex=1.5)
    
    dev.off()
  }
}
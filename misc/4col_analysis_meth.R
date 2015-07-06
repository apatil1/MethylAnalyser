
sheet <- "/Users/ashwini/Desktop/michelle_data/final/sheet.csv"
type <- "TSS"
folder <- "/Users/ashwini/Desktop/new-annotation/"
analysis_methylation4(beta.tab_CpGI.PlusMinus1500_TSS, sheet, folder, type)
#Function for analysis of methylation data and to generate genelist and delta plots


analysis_methylation4 <- function(tss, sheet, Dir, type) {
  
  
  #Read the sheet consisting of questions
  dat <- read.csv(sheet, header=T)
  
  summary <- NULL #For summary of the analysis
  rname <- NULL #For rownames of summary of the analysis
  subDir_file <- "files"
  
  #Create directory for plots and files generated from the analysis
  dir.create(file.path(Dir, subDir_file), showWarnings = FALSE)
  
  for(i in 1:length(dat$x1)) {
    
    
    #Get CpG island probes that get methylated specifically upon fur exposure 
    #abs() <= 0.05 will pick probes that do not undergo much changes in methylation upon culturing/ageing.
    age.Stable_fur.hypmeth.probes <- tss[abs(tss[,colnames(tss)==dat$x1[i]] - tss[,colnames(tss)==dat$x2[i]])  <= 0.05 & tss[,colnames(tss)==dat$y1[i]] - tss[,colnames(tss)==dat$y2[i]]  >= 0.2, ]
    age.Stable_fur.hypometh.probes <- tss[abs(tss[,colnames(tss)==dat$x1[i]] - tss[,colnames(tss)==dat$x2[i]]) <= 0.05 & tss[,colnames(tss)==dat$y1[i]] - tss[,colnames(tss)==dat$y2[i]] <= -0.2, ]
    
   
    
    if (length(rownames(age.Stable_fur.hypmeth.probes)) > 0) {
      age.Stable_fur.hypmeth.probes <- cbind.data.frame(age.Stable_fur.hypmeth.probes, "Gene"=rep(NA, times=nrow(age.Stable_fur.hypmeth.probes)))
      for(j in 1:nrow(age.Stable_fur.hypmeth.probes)) {
        probe.id <- rownames(age.Stable_fur.hypmeth.probes)[j]
        gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
        if(length(gene) == 1) { age.Stable_fur.hypmeth.probes$Gene[j] <- as.character(gene) }
      }
      write.csv(age.Stable_fur.hypmeth.probes, paste(Dir, subDir_file,"/", type, "_age.stable.hyper","_",dat$x1[i],"-",dat$x2[i],"-",dat$y1[i],"-",dat$y2[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }
    
    if (length(rownames(age.Stable_fur.hypometh.probes)) > 0) {
      age.Stable_fur.hypometh.probes <- cbind.data.frame(age.Stable_fur.hypometh.probes, "Gene"=rep(NA, times=nrow(age.Stable_fur.hypometh.probes)))
      for(j in 1:nrow(age.Stable_fur.hypometh.probes)) {
        probe.id <- rownames(age.Stable_fur.hypometh.probes)[j]
        gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
        if(length(gene) == 1) { age.Stable_fur.hypometh.probes$Gene[j] <- as.character(gene) }
      }
      write.csv(age.Stable_fur.hypometh.probes, paste(Dir, subDir_file,"/", type, "_age.stable.hypo","_",dat$x1[i],"-",dat$x2[i],"-",dat$y1[i],"-",dat$y2[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }
    
    
    #Get probes that range in methylation between +/-0.05 to +/-0.2 during aging/proliferation but are methylated/demethylated >= 0.2 
    age.intermediate_meth.fur.hypmeth.probes <- tss[abs(tss[,colnames(tss)==dat$x1[i]] - tss[,colnames(tss)==dat$x2[i]]) >= 0.05 & abs(tss[,colnames(tss)==dat$x1[i]] - tss[,colnames(tss)==dat$x2[i]]) <= 0.2 & tss[,colnames(tss)==dat$y1[i]] - tss[,colnames(tss)==dat$y2[i]] >= 0.2, ]
    age.intermediate_meth.fur.hypometh.probes <- tss[abs(tss[,colnames(tss)==dat$x1[i]] - tss[,colnames(tss)==dat$x2[i]]) >= 0.05 & abs(tss[,colnames(tss)==dat$x1[i]] - tss[,colnames(tss)==dat$x2[i]]) <= 0.2 & tss[,colnames(tss)==dat$y1[i]] - tss[,colnames(tss)==dat$y2[i]] <= -0.2, ]
    
    if (length(rownames(age.intermediate_meth.fur.hypmeth.probes)) > 0) {
      age.intermediate_meth.fur.hypmeth.probes <- cbind.data.frame(age.intermediate_meth.fur.hypmeth.probes, "Gene"=rep(NA, times=nrow(age.intermediate_meth.fur.hypmeth.probes)))
      for(j in 1:nrow(age.intermediate_meth.fur.hypmeth.probes)) {
        probe.id <- rownames(age.intermediate_meth.fur.hypmeth.probes)[j]
        gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
        if(length(gene) == 1) { age.intermediate_meth.fur.hypmeth.probes$Gene[j] <- as.character(gene) }
      }
      write.csv(age.intermediate_meth.fur.hypmeth.probes, paste(Dir, subDir_file,"/", type, "_age.inter.hyper","_",dat$x1[i],"-",dat$x2[i],"-",dat$y1[i],"-",dat$y2[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }
    
    if (length(rownames(age.intermediate_meth.fur.hypometh.probes)) > 0) {
      age.intermediate_meth.fur.hypometh.probes <- cbind.data.frame(age.intermediate_meth.fur.hypometh.probes, "Gene"=rep(NA, times=nrow(age.intermediate_meth.fur.hypometh.probes)))
      for(j in 1:nrow(age.intermediate_meth.fur.hypometh.probes)) {
        probe.id <- rownames(age.intermediate_meth.fur.hypometh.probes)[j]
        gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
        if(length(gene) == 1) { age.intermediate_meth.fur.hypometh.probes$Gene[j] <- as.character(gene) }
      }
      write.csv(age.intermediate_meth.fur.hypometh.probes, paste(Dir, subDir_file,"/", type, "_age.inter.hypo","_",dat$x1[i],"-",dat$x2[i],"-",dat$y1[i],"-",dat$y2[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }
    
    
  } 
    
  
  
 
}


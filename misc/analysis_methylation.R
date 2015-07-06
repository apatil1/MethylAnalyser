
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
        gene <- probe.to.genesymbol[[probe.id]]
        if(length(gene) == 1) { age.hypermethylation$Gene[j] <- gene }
      }
      write.csv(age.hypermethylation, paste(Dir, subDir_file,"/", type, "_age.hyper","_",dat$first[i],"-",dat$second[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }

    if (length(rownames(age.hypomethylation)) > 0) {
      age.hypomethylation <- cbind.data.frame(age.hypomethylation, "Gene"=rep(NA, times=nrow(age.hypomethylation)))
      for(j in 1:nrow(age.hypomethylation)) {
        probe.id <- rownames(age.hypomethylation)[j]
        gene <- probe.to.genesymbol[[probe.id]]
        if(length(gene) == 1) { age.hypomethylation$Gene[j] <- gene }
      }
      write.csv(age.hypomethylation, paste(Dir, subDir_file,"/", type, "_age.hypo","_",dat$first[i],"-",dat$second[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }

    if (length(rownames(age.fur.hypermethylation)) > 0) {
      age.fur.hypermethylation <- cbind.data.frame(age.fur.hypermethylation, "Gene"=rep(NA, times=nrow(age.fur.hypermethylation)))
      for(j in 1:nrow(age.fur.hypermethylation)) {
        probe.id <- rownames(age.fur.hypermethylation)[j]
        gene <- probe.to.genesymbol[[probe.id]]
        if(length(gene) == 1) { age.fur.hypermethylation$Gene[j] <- gene }
      }
      write.csv(age.fur.hypermethylation, paste(Dir, subDir_file,"/", type, "_age.fur.hyper","_",dat$first[i],"-",dat$second[i],"-",dat$further[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }

    if (length(rownames(age.fur.hypomethylation)) > 0) {
      age.fur.hypomethylation <- cbind.data.frame(age.fur.hypomethylation, "Gene"=rep(NA, times=nrow(age.fur.hypomethylation)))
      for(j in 1:nrow(age.fur.hypomethylation)) {
        probe.id <- rownames(age.fur.hypomethylation)[j]
        gene <- probe.to.genesymbol[[probe.id]]
        if(length(gene) == 1) { age.fur.hypomethylation$Gene[j] <- gene }
      }
      write.csv(age.fur.hypomethylation, paste(Dir, subDir_file,"/", type, "_age.fur.hypo","_",dat$first[i],"-",dat$second[i],"-",dat$further[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }

    if (length(rownames(age.Stable_fur.hypmeth.probes)) > 0) {
      age.Stable_fur.hypmeth.probes <- cbind.data.frame(age.Stable_fur.hypmeth.probes, "Gene"=rep(NA, times=nrow(age.Stable_fur.hypmeth.probes)))
      for(j in 1:nrow(age.Stable_fur.hypmeth.probes)) {
        probe.id <- rownames(age.Stable_fur.hypmeth.probes)[j]
        gene <- probe.to.genesymbol[[probe.id]]
        if(length(gene) == 1) { age.Stable_fur.hypmeth.probes$Gene[j] <- gene }
      }
      write.csv(age.Stable_fur.hypmeth.probes, paste(Dir, subDir_file,"/", type, "_age.stable.hyper","_",dat$first[i],"-",dat$second[i],"-",dat$further[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }

    if (length(rownames(age.Stable_fur.hypometh.probes)) > 0) {
      age.Stable_fur.hypometh.probes <- cbind.data.frame(age.Stable_fur.hypometh.probes, "Gene"=rep(NA, times=nrow(age.Stable_fur.hypometh.probes)))
      for(j in 1:nrow(age.Stable_fur.hypometh.probes)) {
        probe.id <- rownames(age.Stable_fur.hypometh.probes)[j]
        gene <- probe.to.genesymbol[[probe.id]]
        if(length(gene) == 1) { age.Stable_fur.hypometh.probes$Gene[j] <- gene }
      }
      write.csv(age.Stable_fur.hypometh.probes, paste(Dir, subDir_file,"/", type, "_age.stable.hypo","_",dat$first[i],"-",dat$second[i],"-",dat$further[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }

    if (length(rownames(age.intermediate_meth.fur.hypmeth.probes)) > 0) {
      age.intermediate_meth.fur.hypmeth.probes <- cbind.data.frame(age.intermediate_meth.fur.hypmeth.probes, "Gene"=rep(NA, times=nrow(age.intermediate_meth.fur.hypmeth.probes)))
      for(j in 1:nrow(age.intermediate_meth.fur.hypmeth.probes)) {
        probe.id <- rownames(age.intermediate_meth.fur.hypmeth.probes)[j]
        gene <- probe.to.genesymbol[[probe.id]]
        if(length(gene) == 1) { age.intermediate_meth.fur.hypmeth.probes$Gene[j] <- gene }
      }
      write.csv(age.intermediate_meth.fur.hypmeth.probes, paste(Dir, subDir_file,"/", type, "_age.inter.hyper","_",dat$first[i],"-",dat$second[i],"-",dat$further[i],".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
    }

    if (length(rownames(age.intermediate_meth.fur.hypometh.probes)) > 0) {
      age.intermediate_meth.fur.hypometh.probes <- cbind.data.frame(age.intermediate_meth.fur.hypometh.probes, "Gene"=rep(NA, times=nrow(age.intermediate_meth.fur.hypometh.probes)))
      for(j in 1:nrow(age.intermediate_meth.fur.hypometh.probes)) {
        probe.id <- rownames(age.intermediate_meth.fur.hypometh.probes)[j]
        gene <- probe.to.genesymbol[[probe.id]]
        if(length(gene) == 1) { age.intermediate_meth.fur.hypometh.probes$Gene[j] <- gene }
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

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
  infi.SOI.filtered.beta<-infi.SOI.filtered.beta[rownames(infi.SOI.filtered.beta) %in% probe.to.chr$Probe_ID,]
  
  
  
  summary <- paste("Total number of probes :", total, "\nTotal number of probes after filtering minfi probe detection p-value > 0.01 :", length(rownames(infi.SOI.filtered.beta)), "\nNo. of probes removed : ", total - length(rownames(infi.SOI.filtered.beta)), sep="")
  
  #write summary table to file
  write.table(summary, paste(folder,"summary_data.txt",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.

  return (infi.SOI.filtered.beta)
}
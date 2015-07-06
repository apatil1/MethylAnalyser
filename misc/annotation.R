############################################################################################################################################################
#ANNOTATION DATA
#Illumina 450K probe annotation using IlluminaHumanMethylation450k.db
############################################################################################################################################################

library(IlluminaHumanMethylation450k.db)

probe.loc.in.gene <- as.list(IlluminaHumanMethylation450kPROBELOCATION) #Map a probe ID with its location in the gene
CpGloc <- as.list(IlluminaHumanMethylation450kCPGILOCATION) #IlluminaHumanMethylation450kCPGILOCATION maps between Illumina probe IDs and the UCSC CpG island features, if any, with which they are associated
island.probes <- unlist(sapply(c(1:length(CpGloc)), FUN=function(i) {return(grep("island" , CpGloc[i], ignore.case=T, value=T))})) #gets island probes
shelf.probes <- unlist(sapply(c(1:length(CpGloc)), FUN=function(i) {return(grep("shelf" , CpGloc[i], ignore.case=T, value=T))})) #gets shelf probes
shore.probes <- unlist(sapply(c(1:length(CpGloc)), FUN=function(i) {return(grep("shore" , CpGloc[i], ignore.case=T, value=T))})) #gets shelf probes

# Get postion of probes relative to gene structure
unique(sapply(probe.loc.in.gene, FUN=function(i){
  strsplit(i, ":" )[[1]][2]
}))# gets the various positional information relative to the gene structure for the probes in the annotaiton data.

tss.probes <- unlist(sapply(c(1:length(probe.loc.in.gene)), FUN=function(i) {return(grep("TSS" , probe.loc.in.gene[i], ignore.case=T, value=T))})) #gets probes that are near TSS (200 or 1500 bp near TSS)
body.probes <- unlist(sapply(c(1:length(probe.loc.in.gene)), FUN=function(i) {return(grep("Body" , probe.loc.in.gene[i], ignore.case=T, value=T))})) #gets probes that are in the body of gene
island.tss.probes <- intersect(names(island.probes), names(tss.probes)) #gets probes near TSS and in CpG island
names(island.tss.probes) <- island.tss.probes # this is done because various plot function below require the input as names vectors. 
length(names(island.probes)); length(names(tss.probes)); length(island.tss.probes)

# Map infinium probe names to gene symbols
x <- IlluminaHumanMethylation450kSYMBOL
mapped_probes <- mappedkeys(x)
probe.to.genesymbol <- as.list(x[mapped_probes])
head(probe.to.genesymbol)

# Map infinium probes to starting position of the gene (is this TSS?)
x <- IlluminaHumanMethylation450kCHRLOC
mapped_probes <- mappedkeys(x)
probe.to.GeneStart <- as.list(x[mapped_probes])
head(probe.to.GeneStart)

# Map infinium probe names to CpG coordinate
x <- IlluminaHumanMethylation450kCPGCOORDINATE
mapped_probes <- mappedkeys(x)
probe.to.CpG_coordinate <- as.list(x[mapped_probes])
head(probe.to.CpG_coordinate)

# Map infinium probe names to chromosomes
x <- IlluminaHumanMethylation450kCHR37
mapped_probes <- mappedkeys(x)
probe.to.chr <- as.data.frame(x[mapped_probes])
head(probe.to.chr)

# Remove  X, Y and MULTI chromosome probes
probe.to.chr <- probe.to.chr[ !probe.to.chr$Chromosome_37 %in% c("X", "Y", ""), ] #removes probes corresponding to X, Y, MULTI and "" chromosomes

unique(probe.to.chr$Chromosome_37)

#save image
save.image("annotation.RData")

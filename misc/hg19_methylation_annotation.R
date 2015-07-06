## load the library
library(FDb.InfiniumMethylation.hg19)
## list the contents that are loaded into memory
ls(FDb.InfiniumMethylation.hg19)
## show the db object that is loaded by calling its name
FDb.InfiniumMethylation.hg19
## extract features for use in constructing SummarizedExperiments
## or comparing chip features against other data (e.g. ChIP-seq)
InfiniumMethylation <- features(FDb.InfiniumMethylation.hg19)

show(InfiniumMethylation)
all.data<-as.data.frame(InfiniumMethylation)
dim(all.data)
all.probes <-all.data[all.data$seqnames %in% c("chr16", "chr1",  "chr15", "chr9",  "chr3",  "chr12", "chr11", "chr8",  "chr2",  "chr19", "chr17", "chr18", "chr4",  "chr10", "chr21", "chr6",  "chr14", "chr7",  "chr13",  "chr5",  "chr22", "chr20"),]
all.probes.xy.filter <- rownames(all.probes)
save(all.probes.xy.filter ,file="all.probes.hg19.RData")

## Example: probes that overlap Irizarrys HMM CpG islands
data(hg19.islands)
CGI.probes <- subsetByOverlaps(InfiniumMethylation, hg19.islands)
cgi.island<-(as.data.frame(CGI.probes))

dim(cgi.island)
cgi.xy.filter <-cgi.island[cgi.island$seqnames %in% c("chr16", "chr1",  "chr15", "chr9",  "chr3",  "chr12", "chr11", "chr8",  "chr2",  "chr19", "chr17", "chr18", "chr4",  "chr10", "chr21", "chr6",  "chr14", "chr7",  "chr13",  "chr5",  "chr22", "chr20"),]
cgi.xy.filter <- rownames(cgi.xy.filter)

save(cgi.xy.filter ,file="cgi.xy.filter.RData")

## Same as above, but now for "shores"
hg19.shores <- c(flank(hg19.islands, 2000, start=TRUE),
                 flank(hg19.islands, 2000, start=FALSE))
shore.probes <- subsetByOverlaps(InfiniumMethylation, hg19.shores)
head(shore.probes)
tail(shore.probes)
shore <- as.data.frame(shore.probes)
dim(shore)
## The same logic works for overlapping probes with other data.
shelf2000-5000
body
#get shelf, body




######################################################################################
#For gene symbol to probe mapping, I'd use the newer annotation package IlluminaHumanMethylation450kanno.ilmn12.hg19.

#source("http://bioconductor.org/biocLite.R")
#biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)

 

#The documentation and methods are not great just now, but you can extract a DataFrame from it like this:

anno <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other
names(anno)

probe.to.gene <- cbind(probe=rownames(anno),gene=anno$UCSC_RefGene_Name)
head(probe.to.gene)
dim(probe.to.gene)
probe.to.gene<- data.frame(probe.to.gene)
save(probe.to.gene ,file="probe.to.gene.RData")
######################################################################################

s
setwd("D:/ONC-Analysis/users/ash/michelle/michelle_tcga")

#load function
source("D:/ONC-Analysis/users/ash/michelle/methylation_analysis/functions/tcga_functions.R")
source("D:/ONC-Analysis/users/ash/michelle/methylation_analysis/functions/venn.R")

load("D:/ONC-Analysis/users/ash/michelle/methylation_analysis/R_image/cgi.xy.filter.RData")
load("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/r_images/data_LUAD.RData")
#load("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/r_images/LUSC.RData")
load("D:/ONC-Analysis/users/ash/michelle/methylation_analysis/R_image/probesWithin5000TSS.RData")

#Remove probes from X and Y Chromosome
#data<-data[rownames(data) %in% probe.to.chr$Probe_ID,]
dim(data)

#to get what sample is what
sample <- sample_type("LUAD_450k.csv")

#Subset the data table to retain CpG-island probes.
dim(data)
head(data)
data_CpGI <- data[rownames(data) %in% cgi.xy.filter , ]
dim(data_CpGI)

#get -500 to 1500 TSS probes
head(probesWithin5000TSS)
tss.1500.probes <- as.character(probesWithin5000TSS[ (probesWithin5000TSS$Dist.to.TSS <= 1500) & (probesWithin5000TSS$Dist.to.TSS >= -500), "ProbeName"])
data_CpGI.TSS <- data_CpGI[rownames(data_CpGI) %in% tss.1500.probes,] #these are the CpGI probes within +/-1500 bp from TSS
data_CpGI.TSS<-data.frame(data_CpGI.TSS)

dim(data_CpGI.TSS)
head(data_CpGI.TSS)

save( data_CpGI.TSS, file="LUAD_cpgi_tss.RData")
#############################################################################

setwd("D:/ONC-Analysis/users/ash/michelle/michelle_tcga")

#load function
source("D:/ONC-Analysis/users/ash/michelle/methylation_analysis/functions/tcga_functions.R")
source("D:/ONC-Analysis/users/ash/michelle/methylation_analysis/functions/venn.R")


load("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD_cpgi_tss.RData")
#load("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUSC_cpgi_tss.RData")
load("D:/ONC-Analysis/users/ash/michelle/methylation_analysis/R_image/probe.to.gene.RData")
#to get what sample is what
sample <- sample_type("LUAD_450k.csv")

#Differntially methylated probes
normal <- data_CpGI.TSS[,grep("0", sample)]
tumor <- data_CpGI.TSS[,grep("1", sample)]

# calculate means of each class for each probe
normal.m <- apply(normal,1,mean,na.rm=T)
tumor.m <- apply(tumor,1,mean,na.rm=T)

head(normal.m)
length(normal.m)
class(normal.m)
head(tumor.m)
length(tumor.m)

# get fold changes 
delta_beta <- tumor.m - normal.m #delta beta 
#fold <- tumor.m/normal.m #fold change of 1 would mean no change

head(delta_beta)
length(delta_beta)
#head(fold)
#length(fold)

###############################################################################

library(limma)
design <-cbind(Grp1=1,Grp2vs1=as.numeric(sample))

fit <- lmFit(data_CpGI.TSS,design)
fit <- eBayes(fit)
z <- topTable(fit,coef=ncol(design))

head(fit$p.value) 
dim(fit$p.value)
#################################################################################

#p value < 0.01
pval.filter <- fit$p.value[fit$p.value[,2] < 0.01,] 
dim(pval.filter)

#hypermethylation
#delta beta > 0.2 and beta normal < 0.2
delta_beta.hyper <- data.frame(delta_beta[delta_beta > 0.2]) 
normal.m.hyper <- data.frame(normal.m[normal.m < 0.2]) 

#probes with p value < 0.01, delta beta > 0.2 and beta normal < 0.2 
tcga.hyper.probes <- rownames(delta_beta.hyper)[( rownames(delta_beta.hyper) %in% rownames(pval.filter)) & (rownames(delta_beta.hyper) %in% rownames(normal.m.hyper))]
tcga.hyper.gene <- probe.to.gene$gene[probe.to.gene$probe %in% tcga.hyper.probes]

#hypomethylation
#delta beta < -0.2 and beta normal > 0.5
delta_beta.hypo <- data.frame(delta_beta[delta_beta < -0.2]) 
normal.m.hypo <- data.frame(normal.m[normal.m > 0.5])

#probes with p value < 0.01, delta beta < 0.2 and beta normal > 0.5 
tcga.hypo.probes <- rownames(delta_beta.hypo)[( rownames(delta_beta.hypo) %in% rownames(pval.filter)) & (rownames(delta_beta.hypo) %in% rownames(normal.m.hypo))]
tcga.hypo.gene <- probe.to.gene$gene[probe.to.gene$probe %in% tcga.hypo.probes]


save( tcga.hyper.probes,tcga.hypo.probes, tcga.hyper.gene, tcga.hypo.gene, file="LUAD_tcga-list.RData")
#######################################################################################


#load("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD_tcga-list.RData")
load("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUSC_tcga-list.RData")

setwd("D:/ONC-Analysis/users/ash/michelle/michelle_tcga")
#get 15m data
stable_hyper <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.stable.hyper_unt_15-unt_1-csc_15.csv")
stable_hypo <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.stable.hypo_unt_15-unt_1-csc_15.csv")
inter_hyper <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.inter.hyper_unt_15-unt_1-csc_15.csv")
inter_hypo <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.inter.hypo_unt_15-unt_1-csc_15.csv")

####################################################################################
#venn diagram gene
title <- "Venn diagram for TCGA LUAD and 15 month CSC treated methylation data "

l <- c("stable hypermethylated genes", "TCGA Hypermethylated genes")
venn2_all(stable_hyper$Gene, tcga.hyper.gene, l, "LUAD_stable_hyper_gene")

l <- c("stable hypomethylated genes", "TCGA Hypomethylated genes")
venn2_all(stable_hypo$Gene, tcga.hypo.gene, l, "LUAD_stable_hypo_gene")

l <- c("intermediate hypermethylated genes", "TCGA Hypermethylated genes")
venn2_all(inter_hyper$Gene, tcga.hyper.gene, l, "LUAD_inter_hyper_gene")

l <- c("intermediate hypomethylated genes", "TCGA Hypomethylated genes")
venn2_all(inter_hypo$Gene, tcga.hypo.gene, l, "LUAD_inter_hypo_gene")

###########################################################
c.hyper <- c(as.character(stable_hyper$Gene),as.character(inter_hyper$Gene))
c.hypo <- c(as.character(stable_hypo$Gene),as.character(inter_hypo$Gene))

title <- "combine inter stable genes "
l <- c("hypermethylated genes", "LUAD TCGA Hypermethylated genes")
venn2_all(c.hyper, luad.hyper, l, "combine LUAD_hyper_gene")

l <- c("hypomethylated genes", "LUAD TCGA Hypomethylated genes")
venn2_all(c.hypo, luad.hypo, l, "combineLUAD_hypo_gene")

c.hyper <- c(as.character(stable_hyper$Gene),as.character(inter_hyper$Gene))
c.hypo <- c(as.character(stable_hypo$Gene),as.character(inter_hypo$Gene))
l <- c("hypermethylated genes", "LUSC TCGA Hypermethylated genes")
venn2_all(c.hyper, lusc.hyper, l, "combineLUSC_hyper_gene")

l <- c("hypomethylated genes", "LUSC TCGA Hypomethylated genes")
venn2_all(c.hypo, lusc.hypo, l, "combineLUSC_hypo_gene")

###################################

#venn diagram gene LUSC
title <- "Venn diagram for TCGA LUSC and 15 month CSC treated methylation data "

l <- c("stable hypermethylated genes", "TCGA Hypermethylated genes")
venn2_all(stable_hyper$Gene, tcga.hyper.gene, l, "LUSC_stable_hyper_gene")

l <- c("stable hypomethylated genes", "TCGA Hypomethylated genes")
venn2_all(stable_hypo$Gene, tcga.hypo.gene, l, "LUSC_stable_hypo_gene")

l <- c("intermediate hypermethylated genes", "TCGA Hypermethylated genes")
venn2_all(inter_hyper$Gene, tcga.hyper.gene, l, "LUSC_inter_hyper_gene")

l <- c("intermediate hypomethylated genes", "TCGA Hypomethylated genes")
venn2_all(inter_hypo$Gene, tcga.hypo.gene, l, "LUSC_inter_hypo_gene")
####################################################################################

#bivalent genes
file.bi <- read.table("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/ES_chromatin genes.txt", sep="", header=T)
head(file.bi)
dim(file.bi)
colnames(file.bi)
col_bi <- grep("GEOandKu_MACS.ES.K4K27.genes", colnames(file.bi))

for (i in 1:length(col_bi)) {
  if(i==1) {
    bivalent.genes <- as.character(file.bi[,col_bi[i]])
  }else {
    bivalent.genes <- append(bivalent.genes, as.character(file.bi[,col_bi[i]]))
  }
}

length(bivalent.genes)
head(c)

#########################################
#combine and bivalent
c.hyper <- c(as.character(stable_hyper$Gene),as.character(inter_hyper$Gene))
c.hypo <- c(as.character(stable_hypo$Gene),as.character(inter_hypo$Gene))

title <- "combine bivalent genes "
l <- c("hypermethylated genes", "bivalent genes")
venn2_all(c.hyper, bivalent.genes, l, "combine hyper bivalent gene")

l <- c("hypomethylated genes", "bivalent genes")
venn2_all(c.hypo, bivalent.genes, l, "combine hypo bivalent gene")

#########################################################################################
#venn diagram bivalent gene to michelle data
title <- "Venn diagram for bivalent genes and 15 month CSC treated methylation data "

l <- c( "bivalent genes","stable hypermethylated genes")
venn2_all( bivalent.genes,stable_hyper$Gene, l, "bivalent_stable_hyper_gene")

l <- c( "bivalent genes","stable hypomethylated genes")
venn2_all( bivalent.genes,stable_hypo$Gene, l, "bivalent_stable_hypo_gene")

l <- c( "bivalent genes", "intermediate hypermethylated genes")
venn2_all( bivalent.genes,inter_hyper$Gene, l, "bivalent_inter_hyper_gene")

l <- c( "bivalent genes", "intermediate hypomethylated genes")
venn2_all( bivalent.genes,inter_hypo$Gene, l, "bivalent_inter_hypo_gene")

###########################################################################################

#venn diagram bivalent gene to TCGA LUAD
title <- "Venn diagram for TCGA LUAD and bivalent genes "

l <- c("bivalent genes", "TCGA Hypermethylated genes")
venn2_all(bivalent.genes, tcga.hyper.gene, l, "Bivalent_LUAD_hyper_gene")

l <- c("bivalent genes", "TCGA Hypomethylated genes")
venn2_all(bivalent.genes, tcga.hypo.gene, l, "Bivalent_LUAD_hypo_gene")

###########################################################################################

#venn diagram bivalent gene to TCGA LUSC
title <- "Venn diagram for TCGA LUSC and bivalent genes "

l <- c("bivalent genes", "TCGA Hypermethylated genes")
venn2_all(bivalent.genes, tcga.hyper.gene, l, "Bivalent_LUSC_hyper_gene")

l <- c("bivalent genes", "TCGA Hypomethylated genes")
venn2_all(bivalent.genes, tcga.hypo.gene, l, "Bivalent_LUSC_hypo_gene")

####################################################################

#venn diagram for bivalent genes, michelle's genes and TCGA LUSC genes

title <- "Venn diagram for TCGA LUSC and bivalent genes and 15 month CSC treated methylation data "
l <- c("bivalent genes", "TCGA Hypermethylated genes", "stable hypermethylated genes" )
venn3_all(bivalent.genes, tcga.hyper.gene, stable_hyper$Gene, l, "TCGA_LUSC_Bivalent_stable_hyper_gene")

l <- c("bivalent genes", "TCGA Hypermethylated genes", "intermediate hypermethylated genes" )
venn3_all(bivalent.genes, tcga.hyper.gene, inter_hyper$Gene, l, "TCGA_LUSC_Bivalent_inter_hyper_gene")

l <- c("bivalent genes", "TCGA Hypomethylated genes", "intermediate hypomethylated genes" )
venn3_all(bivalent.genes, tcga.hypo.gene, inter_hypo$Gene, l, "TCGA_LUSC_Bivalent_inter_hypo_gene")

l <- c("bivalent genes", "TCGA Hypomethylated genes", "stable hypomethylated genes" )
venn3_all(bivalent.genes, tcga.hypo.gene, stable_hypo$Gene, l, "TCGA_LUSC_Bivalent_stable_hypo_gene")
##################################################################################

#venn diagram for bivalent genes, michelle's genes and TCGA LUAD genes

title <- "Venn diagram for TCGA LUAD and bivalent genes and 15 month CSC treated methylation data "
l <- c("bivalent genes", "TCGA Hypermethylated genes", "stable hypermethylated genes" )
venn3_all(bivalent.genes, tcga.hyper.gene, stable_hyper$Gene, l, "TCGA_LUAD_Bivalent_stable_hyper_gene")

l <- c("bivalent genes", "TCGA Hypermethylated genes", "intermediate hypermethylated genes" )
venn3_all(bivalent.genes, tcga.hyper.gene, inter_hyper$Gene, l, "TCGA_LUAD_Bivalent_inter_hyper_gene")

l <- c("bivalent genes", "TCGA Hypomethylated genes", "intermediate hypomethylated genes" )
venn3_all(bivalent.genes, tcga.hypo.gene, inter_hypo$Gene, l, "TCGA_LUAD_Bivalent_inter_hypo_gene")

l <- c("bivalent genes", "TCGA Hypomethylated genes", "stable hypomethylated genes" )
venn3_all(bivalent.genes, tcga.hypo.gene, stable_hypo$Gene, l, "TCGA_LUAD_Bivalent_stable_hypo_gene")
######################################################################################

#venn diagram for LUAD genes TCGA LUSC genes
luad.hyper <- tcga.hyper.gene
luad.hypo <- tcga.hypo.gene
lusc.hyper <- tcga.hyper.gene
lusc.hypo <- tcga.hypo.gene

title <- "Venn diagram for TCGA LUAD and LUSC genes"
l <- c("LUAD Hypermethylated genes", "LUSC Hypermethylated genes"  )
venn2_all(luad.hyper, lusc.hyper, l, "TCGA_luad_lusc_hyper_gene")

l <- c("LUAD Hypomethylated genes", "LUSC Hypomethylated genes" )
venn2_all(luad.hypo, lusc.hypo, l, "TCGA_luad_lusc_hypo_gene")

#####################################################################################

#hyper
h.q1.common <- read.table("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/TCGA_luad_lusc_hyper_gene/q1.txt")
h.q2.luad <- read.table("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/TCGA_luad_lusc_hyper_gene/q2.txt")
h.q3.lusc <-  read.table("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/TCGA_luad_lusc_hyper_gene/q3.txt")

#hypo
o.q1.common <- read.table("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/TCGA_luad_lusc_hypo_gene/q1.txt")
o.q2.luad <- read.table("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/TCGA_luad_lusc_hypo_gene/q2.txt")
o.q3.lusc <-  read.table("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/TCGA_luad_lusc_hypo_gene/q3.txt")

#LUAD/LUSC hyper unique stable/inter hyper
#LUAD/LUSC hyper common stable/inter hyper
title <- "common hyper"
l <- c("stable Hypermethylated 15m genes", "common hyper"  )
venn2_all(stable_hyper$Gene, h.q1.common$V1, l, "stable common_hyper")

title <- "LUAD hyper"
l <- c("stable Hypermethylated 15m genes", "unique LUAD hyper"  )
venn2_all(stable_hyper$Gene, h.q2.luad$V1, l, "stable unique LUAD hyper")

title <- "LUSC hyper"
l <- c("stable Hypermethylated 15m genes", "unique LUSC hyper"  )
venn2_all(stable_hyper$Gene, h.q3.lusc$V1, l, "stable unique LUSC hyper")

#inter
title <- "inter common hyper"
l <- c("inter Hypermethylated 15m genes", "common hyper"  )
venn2_all(inter_hyper$Gene, h.q1.common$V1, l, "inter common_hyper")

title <- "inter LUAD hyper"
l <- c("inter Hypermethylated 15m genes", "unique LUAD hyper"  )
venn2_all(inter_hyper$Gene, h.q2.luad$V1, l, "inter unique LUAD hyper")

title <- "inter LUSC hyper"
l <- c("inter Hypermethylated 15m genes", "unique LUSC hyper"  )
venn2_all(inter_hyper$Gene, h.q3.lusc$V1, l, "inter unique LUSC hyper")

dim(stable_hyper)
dim(stable_hypo )
dim(inter_hyper)
dim(inter_hypo)

#same for hypo

title <- "common hypo"
l <- c("stable hypomethylated 15m genes", "common hypo"  )
venn2_all(stable_hypo$Gene, o.q1.common$V1, l, "stable common_hypo")

title <- "LUAD hypo"
l <- c("stable hypomethylated 15m genes", "unique LUAD hypo"  )
venn2_all(stable_hypo$Gene, o.q2.luad$V1, l, "stable unique LUAD hypo")

title <- "LUSC hypo"
l <- c("stable hypomethylated 15m genes", "unique LUSC hypo"  )
venn2_all(stable_hypo$Gene, o.q3.lusc$V1, l, "stable unique LUSC hypo")

#inter
title <- "inter common hypo"
l <- c("inter hypomethylated 15m genes", "common hypo"  )
venn2_all(inter_hypo$Gene, o.q1.common$V1, l, "inter common_hypo")

title <- "inter LUAD hypo"
l <- c("inter hypomethylated 15m genes", "unique LUAD hypo"  )
venn2_all(inter_hypo$Gene, o.q2.luad$V1, l, "inter unique LUAD hypo")

title <- "inter LUSC hypo"
l <- c("inter hypomethylated 15m genes", "unique LUSC hypo"  )
venn2_all(inter_hypo$Gene, o.q3.lusc$V1, l, "inter unique LUSC hypo")


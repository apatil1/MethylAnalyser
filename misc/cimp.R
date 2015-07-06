setwd("D:/ONC-Analysis/users/ash/michelle/michelle_tcga")
#get 15m data
stable_hyper <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.stable.hyper_unt_15-unt_1-csc_15.csv")
stable_hypo <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.stable.hypo_unt_15-unt_1-csc_15.csv")
inter_hyper <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.inter.hyper_unt_15-unt_1-csc_15.csv")
inter_hypo <- read.csv("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/LUAD/files/TSS_age.inter.hypo_unt_15-unt_1-csc_15.csv")

cimp <-read.table("supplTable_CIMP.txt", header=T)
head(cimp)
dim(cimp)

source("D:/ONC-Analysis/users/ash/michelle/methylation_analysis/functions/venn.R")
#venn diagram gene LUSC
title <- "Venn diagram for methylation and cimp data "

l <- c("stable hypermethylated genes", "cimp")
venn2_all(stable_hyper$Gene, cimp$Gene_Symbol, l, "cimp_stable_hyper_gene")

l <- c("stable hypomethylated genes", "cimp")
venn2_all(stable_hypo$Gene,cimp$Gene_Symbol, l, "cimp_stable_hypo_gene")

l <- c("intermediate hypermethylated genes", "cimp")
venn2_all(inter_hyper$Gene, cimp$Gene_Symbol, l, "cimp_inter_hyper_gene")

l <- c("intermediate hypomethylated genes", "cimp")
venn2_all(inter_hypo$Gene, cimp$Gene_Symbol, l, "cimp_inter_hypo_gene")


#combine

c.hyper <- c(as.character(stable_hyper$Gene),as.character(inter_hyper$Gene))
c.hypo <- c(as.character(stable_hypo$Gene),as.character(inter_hypo$Gene))

title <- "combine inter stable genes cimp "
l <- c("hypermethylated genes", "cimp")
venn2_all(c.hyper, cimp$Gene_Symbol, l, "combine-cimp_hyper_gene")

l <- c("hypomethylated genes", "cimp")
venn2_all(c.hypo, cimp$Gene_Symbol, l, "combine-cimp_hypo_gene")

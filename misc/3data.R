rownames(pval.filter)
probes <- data_CpGI.TSS[rownames(data_CpGI.TSS) %in% tcga.hyper.probes, ]

heat_tcga(probes, plotName="0.3hyper_tcga_luad.jpg", title="Heatmap 0.3hyper_tcga_luad ")

#LUSC
title <- "Venn diagram for TCGA LUSC and 15 month CSC treated methylation data "

l <- c("stable hypermethylated genes", "TCGA Hypermethylated genes")
venn2_all(stable_hyper$Gene, tcga.hyper.gene, l, "0.3LUSC_stable_hyper_gene")

l <- c("stable hypomethylated genes", "TCGA Hypomethylated genes")
venn2_all(stable_hypo$Gene, tcga.hypo.gene, l, "0.3LUSC_stable_hypo_gene")

l <- c("intermediate hypermethylated genes", "TCGA Hypermethylated genes")
venn2_all(inter_hyper$Gene, tcga.hyper.gene, l, "0.3LUSC_inter_hyper_gene")

l <- c("intermediate hypomethylated genes", "TCGA Hypomethylated genes")
venn2_all(inter_hypo$Gene, tcga.hypo.gene, l, "0.3LUSC_inter_hypo_gene")

#LUAD
title <- "Venn diagram for TCGA LUAD and 15 month CSC treated methylation data "

l <- c("stable hypermethylated genes", "TCGA Hypermethylated genes")
venn2_all(stable_hyper$Gene, tcga.hyper.gene, l, "0.3LUAD_stable_hyper_gene")

l <- c("stable hypomethylated genes", "TCGA Hypomethylated genes")
venn2_all(stable_hypo$Gene, tcga.hypo.gene, l, "0.3LUAD_stable_hypo_gene")

l <- c("intermediate hypermethylated genes", "TCGA Hypermethylated genes")
venn2_all(inter_hyper$Gene, tcga.hyper.gene, l, "0.3LUAD_inter_hyper_gene")

l <- c("intermediate hypomethylated genes", "TCGA Hypomethylated genes")
venn2_all(inter_hypo$Gene, tcga.hypo.gene, l, "0.3LUAD_inter_hypo_gene")
######################################################################

#venn diagram for LUAD genes TCGA LUSC genes
luad.hyper <- tcga.hyper.gene
luad.hypo <- tcga.hypo.gene
#lusc.hyper <- tcga.hyper.gene
#lusc.hypo <- tcga.hypo.gene

title <- "Venn diagram for TCGA LUAD and LUSC genes"
l <- c("LUAD Hypermethylated genes", "LUSC Hypermethylated genes"  )
venn2_all(luad.hyper, lusc.hyper, l, "0.3TCGA_luad_lusc_hyper_gene")

l <- c("LUAD Hypomethylated genes", "LUSC Hypomethylated genes" )
venn2_all(luad.hypo, lusc.hypo, l, "0.3TCGA_luad_lusc_hypo_gene")

#####################################################################################

#hyper
h.q1.common <- read.table("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/0.3TCGA_luad_lusc_hyper_gene/q1.txt")
h.q2.luad <- read.table("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/0.3TCGA_luad_lusc_hyper_gene/q2.txt")
h.q3.lusc <-  read.table("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/0.3TCGA_luad_lusc_hyper_gene/q3.txt")

#hypo
o.q1.common <- read.table("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/0.3TCGA_luad_lusc_hypo_gene/q1.txt")
o.q2.luad <- read.table("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/0.3TCGA_luad_lusc_hypo_gene/q2.txt")
o.q3.lusc <-  read.table("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/0.3TCGA_luad_lusc_hypo_gene/q3.txt")

#LUAD/LUSC hyper unique stable/inter hyper
#LUAD/LUSC hyper common stable/inter hyper
title <- "common hyper"
l <- c("stable Hypermethylated 15m genes", "common hyper"  )
venn2_all(stable_hyper$Gene, h.q1.common$V1, l, "0.3stable common_hyper")

title <- "LUAD hyper"
l <- c("stable Hypermethylated 15m genes", "unique LUAD hyper"  )
venn2_all(stable_hyper$Gene, h.q2.luad$V1, l, "0.3stable unique LUAD hyper")

title <- "LUSC hyper"
l <- c("stable Hypermethylated 15m genes", "unique LUSC hyper"  )
venn2_all(stable_hyper$Gene, h.q3.lusc$V1, l, "0.3stable unique LUSC hyper")

#inter
title <- "inter common hyper"
l <- c("inter Hypermethylated 15m genes", "common hyper"  )
venn2_all(inter_hyper$Gene, h.q1.common$V1, l, "0.3inter common_hyper")

title <- "inter LUAD hyper"
l <- c("inter Hypermethylated 15m genes", "unique LUAD hyper"  )
venn2_all(inter_hyper$Gene, h.q2.luad$V1, l, "0.3inter unique LUAD hyper")

title <- "inter LUSC hyper"
l <- c("inter Hypermethylated 15m genes", "unique LUSC hyper"  )
venn2_all(inter_hyper$Gene, h.q3.lusc$V1, l, "0.3inter unique LUSC hyper")

dim(stable_hyper$Gene)
dim(stable_hypo )
dim(inter_hyper)
dim(inter_hypo)

#same for hypo


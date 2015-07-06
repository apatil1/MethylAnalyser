source("/Users/ashwini/Desktop/methylation_analysis/functions/venn.R")
source("/Users/ashwini/Desktop/methylation_analysis/functions/tcga_functions.R")

setwd("~/Desktop/new-annotation/files")
#intermediate hypo
a <- read.csv("TSS_age.inter.hypo_unt_15-unt_1-csc_15-csc_1.csv")
b <- read.csv("TSS_age.inter.hypo_unt_10-unt_1-csc_10-csc_1.csv")
c <- read.csv("TSS_age.inter.hypo_unt_6-unt_1-csc_6-csc_1.csv")

setwd("~/Desktop/new-annotation/hari's method")
title <- "methylation inter hypo"
l <- c("TSS_age.inter.hypo_unt_15-unt_1-csc_15-csc_1", "TSS_age.inter.hypo_unt_10-unt_1-csc_10-csc_1" , "TSS_age.inter.hypo_unt_6-unt_1-csc_6-csc_1")
venn3_all(a$Gene, b$Gene, c$Gene, l, "methylation inter hypo")

setwd("~/Desktop/new-annotation/files")
#intermediate hyper
a <- read.csv("TSS_age.inter.hyper_unt_15-unt_1-csc_15-csc_1.csv")
b <- read.csv("TSS_age.inter.hyper_unt_10-unt_1-csc_10-csc_1.csv")
c <- read.csv("TSS_age.inter.hyper_unt_6-unt_1-csc_6-csc_1.csv")

setwd("~/Desktop/new-annotation/hari's method")
title <- "methylation inter hyper"
l <- c("TSS_age.inter.hyper_unt_15-unt_1-csc_15-csc_1", "TSS_age.inter.hyper_unt_10-unt_1-csc_10-csc_1" , "TSS_age.inter.hyper_unt_6-unt_1-csc_6-csc_1")
venn3_all(a$Gene, b$Gene, c$Gene, l, "methylation inter hyper")


setwd("~/Desktop/new-annotation/files")
#stable hypo
a <- read.csv("TSS_age.stable.hypo_unt_15-unt_1-csc_15-csc_1.csv")
b <- read.csv("TSS_age.stable.hypo_unt_10-unt_1-csc_10-csc_1.csv")
c <- read.csv("TSS_age.stable.hypo_unt_6-unt_1-csc_6-csc_1.csv")

setwd("~/Desktop/new-annotation/hari's method")
title <- "methylation stable hypo"
l <- c("TSS_age.stable.hypo_unt_15-unt_1-csc_15-csc_1", "TSS_age.stable.hypo_unt_10-unt_1-csc_10-csc_1" , "TSS_age.stable.hypo_unt_6-unt_1-csc_6-csc_1")
venn3_all(a$Gene, b$Gene, c$Gene, l, "methylation stable hypo")



setwd("~/Desktop/new-annotation/files")
#stable hyper
a <- read.csv("TSS_age.stable.hyper_unt_15-unt_1-csc_15-csc_1.csv")
b <- read.csv("TSS_age.stable.hyper_unt_10-unt_1-csc_10-csc_1.csv")
c <- read.csv("TSS_age.stable.hyper_unt_6-unt_1-csc_6-csc_1.csv")

setwd("~/Desktop/new-annotation/hari's method")
title <- "methylation stable hyper"
l <- c("TSS_age.stable.hyper_unt_15-unt_1-csc_15-csc_1", "TSS_age.stable.hyper_unt_10-unt_1-csc_10-csc_1" , "TSS_age.stable.hyper_unt_6-unt_1-csc_6-csc_1")
venn3_all(a$Gene, b$Gene, c$Gene, l, "methylation stable hyper")


###########################################
#hyper

setwd("~/Desktop/new-annotation/files")
#stable hyper
a <- read.csv("TSS_age.stable.hyper_unt_15-unt_1-csc_15-csc_1.csv")
b <- read.csv("TSS_age.stable.hyper_unt_10-unt_1-csc_10-csc_1.csv")
c <- read.csv("TSS_age.stable.hyper_unt_6-unt_1-csc_6-csc_1.csv")
setwd("~/Desktop/new-annotation/files")
#intermediate hyper
ai <- read.csv("TSS_age.inter.hyper_unt_15-unt_1-csc_15-csc_1.csv")
bi <- read.csv("TSS_age.inter.hyper_unt_10-unt_1-csc_10-csc_1.csv")
ci <- read.csv("TSS_age.inter.hyper_unt_6-unt_1-csc_6-csc_1.csv")

ac<- c(as.character(a$Gene),as.character(ai$Gene))
bc<- c(as.character(b$Gene),as.character(bi$Gene))
cc<- c(as.character(c$Gene),as.character(ci$Gene))

setwd("~/Desktop/new-annotation/hari's method")
title <- "methylation combine"
l <- c("combine.hyper_unt_15-unt_1-csc_15-csc_1", "combine.hyper_unt_10-unt_1-csc_10-csc_1" , "combine.hyper_unt_6-unt_1-csc_6-csc_1")
venn3_all(ac, bc, cc, l, "methylation combine hyper")
############################################################

#hypo

setwd("~/Desktop/new-annotation/files")
#stable hypo
a <- read.csv("TSS_age.stable.hypo_unt_15-unt_1-csc_15-csc_1.csv")
b <- read.csv("TSS_age.stable.hypo_unt_10-unt_1-csc_10-csc_1.csv")
c <- read.csv("TSS_age.stable.hypo_unt_6-unt_1-csc_6-csc_1.csv")
setwd("~/Desktop/new-annotation/files")
#intermediate hypo
ai <- read.csv("TSS_age.inter.hypo_unt_15-unt_1-csc_15-csc_1.csv")
bi <- read.csv("TSS_age.inter.hypo_unt_10-unt_1-csc_10-csc_1.csv")
ci <- read.csv("TSS_age.inter.hypo_unt_6-unt_1-csc_6-csc_1.csv")

ac<- c(as.character(a$Gene),as.character(ai$Gene))
bc<- c(as.character(b$Gene),as.character(bi$Gene))
cc<- c(as.character(c$Gene),as.character(ci$Gene))

setwd("~/Desktop/new-annotation/hari's method")
title <- "methylation combine"
l <- c("combine.hypo_unt_15-unt_1-csc_15-csc_1", "combine.hypo_unt_10-unt_1-csc_10-csc_1" , "combine.hypo_unt_6-unt_1-csc_6-csc_1")
venn3_all(ac, bc, cc, l, "methylation combine hypo")

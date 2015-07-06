#heatmaps

load("D:/ONC-Analysis/users/ash/michelle/michelle_tcga/r_images/michelle.RData")
length(tcga.hyper.probes)
length(tcga.hypo.probes)

#LUAD
dim(beta.tab)
head(beta.tab)

#michelle data
hyper.michelle <- beta.tab[rownames(beta.tab) %in% tcga.hyper.probes, ]
hypo.michelle <- beta.tab[rownames(beta.tab) %in% tcga.hypo.probes, ]

head(hyper.michelle)
dim(hyper.michelle)

head(hypo.michelle)
dim(hypo.michelle)

head(data_CpGI.TSS)
dim(data_CpGI.TSS)

#tcga data
tcga.hyper.data <- data_CpGI.TSS[rownames(data_CpGI.TSS) %in% tcga.hyper.probes, ]
tcga.hypo.data <- data_CpGI.TSS[rownames(data_CpGI.TSS) %in% tcga.hypo.probes, ]


#combine data
hyper_all <-  merge(tcga.hyper.data, hyper.michelle ,by="row.names")
hypo_all <-  merge(tcga.hypo.data, hypo.michelle ,by="row.names")

rownames(hyper_all) <- hyper_all$Row.names
hyper_all$Row.names <- NULL

rownames(hypo_all) <- hypo_all$Row.names
hypo_all$Row.names <- NULL

#only 15m
hyper_all <- hyper_all[, !colnames(hyper_all) %in% c("csc_10", "unt_10" ,"unt_1" , "csc_1",  "unt_6" , "csc_6", "unt_15","csc20" )]
hypo_all <- hypo_all[, !colnames(hypo_all) %in% c("csc_10", "unt_10" ,"unt_1" , "csc_1",  "unt_6" , "csc_6", "unt_15","csc20" )]


head(hyper_all)
dim(hyper_all)

head(hypo_all)
dim(hypo_all)


#add 9 samples from michelle data to sample for identification michelle: unt=3(green) csc=4(blue), tumor=1 (black), recurrent tumor=1 (black), normal=0 (white)
sample <- append(sample, c(2))

#heatmap
#heat_tcga(hyper_all, plotName="0.3only_15mcsc_opp_heatmap_Hyper_LUAD.jpg", title="Heatmap for probes Hypermethylation LUAD and CSC 15m")
#heat_tcga(hypo_all, plotName="0.3only_15mcsc_opp_heatmap_Hypo_LUAD.jpg", title="Heatmap for probes Hypomethylation LUAD and CSC 15m")

heat_tcga(hyper_all, plotName="0.3only_15mcsc_opp_heatmap_Hyper_LUSC.jpg", title="Heatmap for probes Hypermethylation LUSC and CSC 15m")
heat_tcga(hypo_all, plotName="0.3only_15mcsc_opp_heatmap_Hypo_LUSC.jpg", title="Heatmap for probes Hypomethylation LUSC and CSC 15m")


#################
sample<- c("4","3","3","4","3","4","2","3","4")
#heat_tcga(hyper.michelle, plotName="opp_TCGA_heatmap_Hyper_LUSC.jpg", title="Heatmap for probes Hypermethylation in TCGA LUSC from experiment")
#heat_tcga(hypo.michelle, plotName="opp_TCGA_heatmap_Hypo_LUSC.jpg", title="Heatmap for probes Hypomethylation in TCGA LUSC from experiment")

heat_tcga(hyper.michelle, plotName="opp_TCGA_heatmap_Hyper_LUAD.jpg", title="Heatmap for probes Hypermethylation in TCGA LUAD from experiment")
heat_tcga(hypo.michelle, plotName="opp_TCGA_heatmap_Hypo_LUAD.jpg", title="Heatmap for probes Hypomethylation in TCGA LUAD from experiment")



######
#only tcga
#to get what sample is what
sample <- sample_type("LUAD_450k.csv")

#heat_tcga(tcga.hyper.data, plotName="opp_onlyTCGA_heatmap_Hyper_LUSC.jpg", title="Heatmap for probes Hypermethylation in TCGA LUSC")
#heat_tcga(tcga.hypo.data, plotName="opp_onlyTCGA_heatmap_Hypo_LUSC.jpg", title="Heatmap for probes Hypomethylation in TCGA LUSC")

heat_tcga(tcga.hyper.data, plotName="opp_onlyTCGA_heatmap_Hyper_LUAD.jpg", title="Heatmap for probes Hypermethylation in TCGA LUAD")
heat_tcga(tcga.hypo.data, plotName="opp_onlyTCGA_heatmap_Hypo_LUAD.jpg", title="Heatmap for probes Hypomethylation in TCGA LUAD")








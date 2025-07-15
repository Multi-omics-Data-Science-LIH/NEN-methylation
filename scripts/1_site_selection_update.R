####################################################
##           Extract the unique sites list        ##
##   for most of the tissues(except. minor ones   ##
##    Yue.zhang@lih.lu          2022March         ##
####################################################



# Read-in all the .csv files of the Localizaiton based differential expression


## Set env.
rm(list=ls())
gc()



# set folder
path =".../nen_project/medcom" 
# dir.create(path)
setwd(path)
getwd()




### discard the previous method, use another way to better select the cmps needed

appendix <- read.csv(".../data/original_reprocess/differential_methylation_data/diffMethTable_site_cmp1.csv", header = T)
colon <- read.csv(".../data/original_reprocess/differential_methylation_data/diffMethTable_site_cmp2.csv", header = T)
# Duodenum <- read.csv(".../data/original_reprocess/differential_methylation_data/diffMethTable_site_cmp3.csv", header = T)
Ileum <- read.csv(".../data/original_reprocess/differential_methylation_data/diffMethTable_site_cmp4.csv", header = T)
# Liver <- read.csv(".../data/original_reprocess/differential_methylation_data/diffMethTable_site_cmp5.csv", header = T)
lung <- read.csv(".../data/original_reprocess/differential_methylation_data/diffMethTable_site_cmp6.csv", header = T)
lymph <- read.csv(".../data/original_reprocess/differential_methylation_data/diffMethTable_site_cmp7.csv", header = T)
pancreas <- read.csv(".../data/original_reprocess/differential_methylation_data/diffMethTable_site_cmp8.csv", header = T)
# Papilla <- read.csv(".../data/original_reprocess/differential_methylation_data/diffMethTable_site_cmp9.csv", header = T)
# Pleura <- read.csv(".../data/original_reprocess/differential_methylation_data/diffMethTable_site_cmp10.csv", header = T)
rectum <- read.csv(".../data/original_reprocess/differential_methylation_data/diffMethTable_site_cmp11.csv", header = T)
skin <- read.csv(".../data/original_reprocess/differential_methylation_data/diffMethTable_site_cmp12.csv", header = T)
# Soft_tissue <- read.csv(".../data/original_reprocess/differential_methylation_data/diffMethTable_site_cmp13.csv", header = T)
stomach <- read.csv(".../data/original_reprocess/differential_methylation_data/diffMethTable_site_cmp14.csv", header = T)

# First to extra the sites with an adj.fdr<0.05

appendix_fdr <- appendix[appendix$diffmeth.p.adj.fdr < 0.05, "cgid"]
colon_fdr <- colon[colon$diffmeth.p.adj.fdr < 0.05, "cgid"]
rectum_fdr <- rectum[rectum$diffmeth.p.adj.fdr < 0.05, "cgid"]
Ileum_fdr <- Ileum[Ileum$diffmeth.p.adj.fdr < 0.05, "cgid"]
stomach_fdr <- stomach[stomach$diffmeth.p.adj.fdr < 0.05, "cgid"]
lung_fdr <- lung[lung$diffmeth.p.adj.fdr < 0.05, "cgid"]
lymph_fdr <- lymph[lymph$diffmeth.p.adj.fdr < 0.05, "cgid"]
pancreas_fdr <- pancreas[pancreas$diffmeth.p.adj.fdr < 0.05, "cgid"]
skin_fdr <- skin[skin$diffmeth.p.adj.fdr < 0.05, "cgid"]

# Create the list for setdiff

listinput <- list(appendix_fdr,colon_fdr,rectum_fdr,Ileum_fdr,
                  stomach_fdr,lung_fdr,lymph_fdr,pancreas_fdr,
                  skin_fdr)

ob <- list()
for(i in 1:length(listinput)){
  a <- setdiff(listinput[[i]], unlist(listinput[-i]))
  ob[[i]] <- a
  
}
ob <- setNames(ob,c("appendix_fdr","colon_fdr","rectum_fdr","Ileum_fdr",
                    "stomach_fdr","lung_fdr","lymph_fdr",
                    "pancreas_fdr","skin_fdr"))


# In ob, unique sites are selected
# Go back to big table to exact the features

# for each table/Location, only keep the unique ones


appendix_uni <- appendix[appendix$cgid %in% ob[[1]],]
appendix_uni <- appendix_uni[order(abs(appendix_uni$mean.diff), decreasing = TRUE),]  
appendix_top <- appendix_uni[1:500,"cgid"]


colon_uni <- colon[colon$cgid %in% ob[[2]],]
colon_uni <- colon_uni[order(abs(colon_uni$mean.diff), decreasing = TRUE),]  
colon_top <- colon_uni[1:500,"cgid"]


rectum_uni <- rectum[rectum$cgid %in% ob[[3]],]
rectum_uni <- rectum_uni[order(abs(rectum_uni$mean.diff), decreasing = TRUE),]  
rectum_top <- rectum_uni[1:500,"cgid"]


Ileum_uni <- Ileum[Ileum$cgid %in% ob[[4]],]
Ileum_uni <- Ileum_uni[order(abs(Ileum_uni$mean.diff), decreasing = TRUE),]  
Ileum_top <- Ileum_uni[1:500,"cgid"]


stomach_uni <- stomach[stomach$cgid %in% ob[[5]],]
stomach_uni <- stomach_uni[order(abs(stomach_uni$mean.diff), decreasing = TRUE),]  
stomach_top <- stomach_uni[1:500,"cgid"]


lung_uni <- lung[lung$cgid %in% ob[[6]],]
lung_uni <- lung_uni[order(abs(lung_uni$mean.diff), decreasing = TRUE),]  
lung_top <- lung_uni[1:500,"cgid"]


lymph_uni <- lymph[lymph$cgid %in% ob[[7]],]
lymph_uni <- lymph_uni[order(abs(lymph_uni$mean.diff), decreasing = TRUE),]  
lymph_top <- lymph_uni[1:500,"cgid"]


pancreas_uni <- pancreas[pancreas$cgid %in% ob[[8]],]
pancreas_uni <- pancreas_uni[order(abs(pancreas_uni$mean.diff), decreasing = TRUE),]  
pancreas_top <- pancreas_uni[1:500,"cgid"]

skin_uni <- skin[skin$cgid %in% ob[[9]],]
skin_uni <- skin_uni[order(abs(skin_uni$mean.diff), decreasing = TRUE),]  
skin_top <- skin_uni[1:500,"cgid"]

# get the union of those combined with the original 5k and unique
library(RnBeads)
library(DecompPipeline)

rnb.original_processed_updated_anno_012025_hg19 <- readRDS(".../data/rnb.original_processed_updated_anno_012025_hg19.RDS")

# select a subset of CpGs sites used for MeDeCom analysis
cg_subset <- prepare.CG.subsets(rnb.set = rnb.original_processed_updated_anno_012025_hg19,
                                marker.selection = "var",
                                n.markers = 5000)
names(cg_subset)

meth_5k <- meth(rnb.original_processed_updated_anno_012025_hg19, row.names=TRUE, i = cg_subset$var)
meth_5k <- na.omit(meth_5k)
metht <- t(meth_5k)
id_5k <- colnames(metht)


# toplist reduce

bigfur <-Reduce(union, list(appendix_top, colon_top,lung_top,rectum_top,Ileum_top,
                            stomach_top,lymph_top,pancreas_top,skin_top,id_5k))

sample(1:length(rownames(appendix)), 3000)
random_sites <- appendix[sample(1:length(rownames(appendix)), 3000), "cgid"]

bigfur_random <- union(bigfur,random_sites)
# 12291 CpGs

saveRDS(bigfur_random,file = "selected12k_sites_update", compress = F)

save.image(file='sites_selection_update.RData')



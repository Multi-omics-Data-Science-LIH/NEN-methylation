#########################################################
##   Liver/NEN project RF/other ML data preparation    ##
##         Yue.zhang@lih.lu      Jan. 2022             ##
##        Test with small number of samples            ##
#########################################################
rm(list=ls())
gc()



# set folder
path =".../nen_project/medcom"   
# dir.create(path)
setwd(path)
getwd()

library(tidyverse)
library(RnBeads)
library(DecompPipeline)
library(ggcorrplot)


set.seed(33)
## Load the medecom result file from Decomppipeline.R

md.res <- readRDS(file=".../data/md.res_12k_val_012025.rds")
plotParameters(md.res, K=10,lambdaScale=0.00001)

lmcs<-getLMCs(md.res, K=10, lambda=0.00001)





# load the Rnbeads set from previous analysis(preprocessing)
# rnb.set.prepro <- load.rnb.set(path = ".../data/rnb_prepro.zip",
#                                temp.dir = "/.../data/temp_extraction/")

rnb.set.prepro <- readRDS(".../data/rnb.original_processed_updated_anno_012025_hg19.RDS")

# select a subset of CpGs sites used for MeDeCom analysis
cg_subset <- prepare.CG.subsets(rnb.set = rnb.set.prepro,
                                marker.selection = "var",
                                n.markers = 5000)
names(cg_subset)
str(cg_subset)

meth.data <- meth(rnb.set.prepro, row.names = T, i = cg_subset$var)
meth.data.t <- as.data.frame(t(as.data.frame(meth.data)))


sample.annotation <- ".../data/Clean_up_annotation_updated2025.csv"
annotation <- read.csv(sample.annotation, header = TRUE)
annotation$index <- 1:198

meth.data.t$index <- 1:198

meth.data.ml <- left_join(annotation, meth.data.t, by = "index")

write.csv(meth.data.ml, ".../data/NEN_ml_5k_org.csv", row.names = F)



####################################################
# load the Rnbeads set from previous analysis(preprocessing)
# Extract the methylation table with row names
meth_table <- meth(rnb.set.prepro, row.names = TRUE)


# Update the rnb.set.prepro object to include only the rows present in the cleaned methylation table
rnb.set.prepro <- remove.sites(rnb.set.prepro, !complete.cases(meth_table))
meth_table <- meth_table[complete.cases(meth_table),]
# Get proportion for LMC-10
## ((We choosed LMC10, because it shows more correlation with LUMP.))
prop <- getProportions(md.res, K=10, lambda=0.00001)
colnames(prop) <- annotation$ID


prop <- as.data.frame(prop)
prop_t <- t(prop)
prop_t <- as.data.frame(prop_t)

#redefine row and column names
rownames(prop_t) <- colnames(prop)
colnames(prop_t) <- rownames(prop)

# re-assign an index column for left_join

prop_t$ID <- as.integer(rownames(prop_t))

df_plotting <- inner_join(annotation, prop_t, by = "ID")
lump <- rnb.execute.lump(rnb.set.prepro)
df_plotting$lump <- lump

corrdf <- cor(df_plotting[13:23])
p9 <- ggcorrplot(corrdf, type = "lower",
                 lab = TRUE,  method = "circle")

print(p9)



#############################################################
## 12k + validation CpGs

cg_subset_new <- md.res@parameters$GROUP_LISTS$var
# Extract methylation data for the randomly selected CpGs
meth.data_val <- meth(rnb.set.prepro, row.names = TRUE, i = cg_subset_new)
dim(meth.data_val)


### OR use this secend option to create meth table based on 12k and validation data (Just to be sure we have same CpGs)
rnb.set_validation_processed  <- readRDS(".../data/rnb.set_validation_processed.RDS")


selected_sites <- readRDS(".../data/original_reprocess/selected12k_sites_update")  # Sites' names, need to be converted to indices

# Combine the two datasets: 12k + validation CpGs

# Extract CpG site IDs from rnb.set_validation_processed
cpg_sites_rnb <- rownames(rnb.set_validation_processed@sites)
cpg_sites_cleaned <- gsub("_.*", "", cpg_sites_rnb)

# Ensure selected_sites is unique
selected_sites <- unique(selected_sites)

# Combine the CpG sites: union or intersect
# Union: Combine all unique CpGs from both datasets
combined_cpgs <- unique(c(cpg_sites_cleaned, selected_sites))

# Intersect: Only keep common CpGs between the two datasets
common_cpgs <- intersect(cpg_sites_cleaned, selected_sites)

# Step 1: Find the common CpG sites
common_cpgs_found <- intersect(rownames(meth_table), common_cpgs)

# Display the number of common CpGs found
cat("Number of common CpGs found:", length(common_cpgs_found), "\n")

# (Optional) Step 2: Subset meth_table to include only the common CpGs
meth_table_common <- meth_table[common_cpgs_found, ]

# Verify the subset
dim(meth_table_common) 

###########################################################
# Calculate the variance for each row
row_variances <- apply(meth.data_val, 1, var)

# Order the rows by variance in decreasing order and get the top 5000 indices
top5000_indices <- order(row_variances, decreasing = TRUE)[1:5000]

# Subset the original data to include only the top 5000 most variable sites
meth.data_val_top5000 <- meth.data_val[top5000_indices, ]

# Optionally, you can check the dimensions to confirm
dim(meth.data_val_top5000)  # Should return 5000 x 199

common_for_med <- intersect(common_cpgs, rownames(meth.data_val_top5000))
cat("Number of common CpGs found:", length(common_for_med), "\n")
###################################################

LMC3 <- prop_t$V3
meth.data <- as.data.frame(meth.data_val_top5000)



## Residual function

residuals <- function(meth, pheno, adjustment_columns=c()){
  
  residual_values <- apply(meth, 1, function(x) {
    pheno$x <- x
    form <- as.formula(paste0("x ~ ", paste(adjustment_columns, collapse="+")))
    mod<-lm(form, data = pheno)
    pval<-pf(summary(mod)$fstatistic[1], summary(mod)$fstatistic[2], summary(mod)$fstatistic[3], lower.tail = FALSE)
    if(pval<0.01){
      stats::residuals(mod)
    }else{
      scale(x, scale = FALSE)
    }
  })
  return(residual_values)
  
}

reka_test <- residuals(meth = meth.data, df_plotting, adjustment_columns=c("LMC3"))


#  reka test
reka_test <- as.data.frame(reka_test)
# Format the matrix and pass to python
reka_test$index <- 1:198
meth.data.res.ml <- left_join(annotation,reka_test, by = "index")
meth.data.res.ml.nona <- meth.data.res.ml[, colSums(is.na(meth.data.res.ml))==0]

write.csv(meth.data.res.ml.nona, ".../data/NEN_ml_5k_res_test_update.csv", row.names = F)


############################
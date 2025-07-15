rm(list=ls())
gc()


# set folder
path =".../nen_project/medcom"          # Needs to be modified
# dir.create(path)
setwd(path)
getwd()



library(RnBeads)
library(R.utils)
library(bigstatsr)
library(DecompPipeline)
library(RnBeads.hg38)



# Load the prepared rnb.set
rnb.original_processed_updated_anno_012025_hg19 <- readRDS(".../data/rnb.original_processed_updated_anno_012025_hg19.RDS")

# Load the prepared rnb.set
# rnb.set.prepro <- load.rnb.set(path = ".../data/rnb_prepro.zip",
#                                temp.dir = ".../data/temp_extraction/")

rnb.set_validation_processed <- readRDS(".../data/rnb.set_validation_processed.RDS")


selected_sites <- readRDS(".../data/selected12k_sites_update")  # Sites' names, need to be converted to indices

# Combine the two datasets: 12k + validation CpGs

# Extract CpG site IDs from rnb.set_validation_processed
cpg_sites_rnb <- rownames(rnb.set_validation_processed@sites)
cpg_sites_cleaned <- gsub("_.*", "", cpg_sites_rnb)

# Ensure selected_sites is unique
selected_sites <- unique(selected_sites)

# Combine the CpG sites: union or intersect
# Union: Combine all unique CpGs from both datasets
# combined_cpgs <- unique(c(cpg_sites_cleaned, selected_sites))

# Intersect: Only keep common CpGs between the two datasets
common_cpgs <- intersect(cpg_sites_cleaned, selected_sites)

# ## Remove NA values 
# # Extract the methylation table with row names
meth_table <- meth(rnb.original_processed_updated_anno_012025_hg19, row.names = TRUE)
# # Remove rows with missing values in the methylation table
# #clean_meth_table <- meth_table[complete.cases(meth_table), ]
# 
# # Update the rnb.set.prepro object to include only the rows present in the cleaned methylation table
rnb.set.prepro <- remove.sites(rnb.original_processed_updated_anno_012025_hg19, !complete.cases(meth_table))
meth_table <- meth_table[complete.cases(meth_table),]

# selected sites index
selected_sites_index <- which(rownames(meth_table) %in% common_cpgs)
selected_sites_index <- list(var = selected_sites_index)  # To keep the same format as the prepare_cg() function output

## check for NA values

# Check the number of NA values in the methylation matrix of rnb.set.prepro
na_count_rnb_set <- sum(is.na(meth(rnb.set.prepro)))
cat("Number of NA values in rnb.set.prepro:", na_count_rnb_set, "\n")

# Check the number of NA values in meth_table
na_count_meth_table <- sum(is.na(meth_table))
cat("Number of NA values in meth_table:", na_count_meth_table, "\n")

na_count_selected_sites <- sum(is.na(selected_sites_index$var))
cat("Number of NA values in selected_sites_index$var:", na_count_selected_sites, "\n")

## MEDECOME

md.res <- start.medecom.analysis(
  rnb.set = rnb.set.prepro,
  cg.groups = selected_sites_index,
  Ks = 7:14,
  lambda.grid = c(0, 10^-(2:5)),
  factorviz.outputs = TRUE,
  analysis.name = "NEN_12k_EPIC_update",
  cores = 10
)


# Define the path

md.res_12k_val <- md.res
save_path <- ".../data/md.res_12k_val_012025.rds"

# Save the MeDeComSet object
saveRDS(md.res_12k_val, file = save_path)

# Confirm that the file is saved
cat("MeDeComSet object saved to:", save_path)


rm(list=ls())
gc()



# set folder
path =".../nen_project/medcom"   
# dir.create(path)
setwd(path)
getwd()



## All cpgs
library(tidyverse)
library(RnBeads)
library(DecompPipeline)
library(xgboost)

# load the Rnbeads set from previous analysis(preprocessing)
# rnb.set.prepro <- load.rnb.set(path = ".../data/rnb_prepro.zip",
#                                temp.dir = ".../data/temp_extraction/")
rnb.set.prepro <- readRDS(".../data/rnb.original_processed_updated_anno_012025_hg19.RDS")

md.res_12k_val <- readRDS(file=".../data/md.res_12k_val_012025.rds")
sample.annotation <- ".../data/Clean_up_annotation_updated2025.csv"



## xgboost model 

## 1. Extract the subset of CpGs from the MeDeComSet
cg_subset_new <- md.res_12k_val@parameters$GROUP_LISTS$var
str(cg_subset_new)
# This should be the vector of CpG indices (or IDs) used in MeDeCom analysis.

## 2. Get the methylation data from the RnBeads object
#    using the same CpG indices so it matches the MeDeCom subset.

# Set a seed for reproducibility
set.seed(123)  


# Extract methylation data for the randomly selected CpGs
meth.data_val <- meth(rnb.set.prepro, row.names = TRUE, i = cg_subset_new)
# Transpose the methylation data
meth.data_val.t <- as.data.frame(t(as.data.frame(meth.data_val)))

# Inspect the structure
str(meth.data_val.t)



annotation <- read.csv(sample.annotation, header = TRUE)
annotation$index <- 1:198
meth.data_val.t$index <- 1:198


meth.data_val.ml <- left_join(annotation, meth.data_val.t, by = "index")

# Prepare the data
# Ensure that 'Class' is a factor
meth.data_val.t$Class <- as.factor(meth.data_val.ml$P_grouping)

# Split data into features and target
features <- meth.data_val.t[, !(names(meth.data_val.t) %in% c("index", "Class"))]
target <- meth.data_val.t$Class

x <- as.matrix(features)
y <- as.numeric(target) - 1
########
# Check for NA values
# Check if there are any NA values in the dataset
any_na <- any(is.na(x))
cat("Are there any NA values in the dataset? ", any_na, "\n")

# Count the total number of NA values
total_na <- sum(is.na(x))
cat("Total number of NA values: ", total_na, "\n")

# Identify the rows and columns containing NA values
na_rows <- which(apply(x, 1, function(row) any(is.na(row))))
na_cols <- which(apply(x, 2, function(col) any(is.na(col))))

cat("Rows with NA values: ", na_rows, "\n")
cat("Columns with NA values: ", na_cols, "\n")
#############################

#Remove NA
# Remove columns with any NA values
# x_clean <- x[, colSums(is.na(x)) == 0]
# 
# # Check the dimensions of the cleaned dataset
# cat("Original dataset dimensions: ", dim(x), "\n")
# cat("Cleaned dataset dimensions: ", dim(x_clean), "\n")



# complete_cases <- complete.cases(features)
# features <- features[complete_cases, ]
# target <- target[complete_cases]


# Prepare data for XGBoost
# x <- as.matrix(x_clean)
# y <- as.numeric(target) - 1  # XGBoost requires numeric labels starting from 0

# Train a Random Forest model
set.seed(123)
xgb_model <- xgboost(data = x, label = y, objective = "multi:softprob",
                     num_class = length(levels(target)),
                     nrounds = 800, verbose = 0)

# Get feature importance
importance_matrix <- xgb.importance(feature_names = colnames(x), model = xgb_model)
important_cpgs_xgboost <- importance_matrix$Feature

# Verify the number of selected CpGs
length(important_cpgs_xgboost)  

# Extract methylation data for the important CpGs
# Get all row names (CpG site IDs) from the RnBSet
all_cpgs <- sites(rnb.set.prepro)

# Match the important CpGs with available CpGs
matched_indices <- match(important_cpgs_xgboost, rownames(all_cpgs))

# Remove NA values (for CpGs not found in the dataset)
matched_indices <- matched_indices[!is.na(matched_indices)]

# Step 1: Extract methylation data for the matched CpGs
meth.data_val_selected <- meth(rnb.set.prepro, i = matched_indices)

# Step 2: Transpose the data so rows = samples and columns = CpG sites
meth.data_val.t <- as.data.frame(t(as.data.frame(meth.data_val_selected)))

# Step 3: Rename columns with CpG names from important_cpgs_xgboost
colnames(meth.data_val.t) <- important_cpgs_xgboost

# Step 4: Add an index column to match annotation rows
meth.data_val.t$index <- 1:nrow(meth.data_val.t)

# Step 5: Merge with the annotation data on "index"
meth.data_val.ml <- left_join(annotation, meth.data_val.t, by = "index")

# Step 6: Save the merged data frame as a CSV file
write.csv(meth.data_val.ml, 
          file = ".../data/NEN_ml_xgboosting_update2.csv",
          row.names = FALSE)


####################################
# ============================================================================
# 1. Load Libraries

library(pheatmap)

# ============================================================================
# 2. Load Data & Annotations
# ============================================================================
# Read proportions & annotation for the validation cohort
annotation            <- read.csv(".../data/Clean_up_annotation_updated2025.csv", header = TRUE)
md.res_val_12k        <- readRDS(".../data/md.res_12k_val_012025.rds")

# ============================================================================
# 3. Compute & Plot Proportions (Discovery Data)
# ============================================================================
# Subset & order annotation for plotting
anno_plot <- annotation %>% 
  select(Localization, Primary, NEN.type, P_grouping)
rownames(anno_plot) <- annotation$ID
anno_plot <- anno_plot[order(anno_plot$P_grouping), ]
# ============================================================================
# 4. Prepare Validation Cohort Annotation & Combine Proportions
# ============================================================================
rnb.validation <- readRDS(".../data/rnb.set_validation_processed.RDS")

anno_plot_validation <- data.frame(
  Localization = rep("liver", nrow(rnb.validation@pheno)),
  Primary      = rep("liver", nrow(rnb.validation@pheno)),
  NEN.type     = rep("", nrow(rnb.validation@pheno)),
  P_grouping   = rep("CUP", nrow(rnb.validation@pheno))
)
proportions_val_LMC10 <- readRDS(".../data/validation/proportions_val_LMC10.RDS")
rownames(anno_plot_validation) <- colnames(proportions_val_LMC10)

# ============================================================================
# 5. Read xg CpGs & Extract Relevant IDs
# ============================================================================
cg_list <- important_cpgs_xgboost

# ============================================================================
# 6. Extract Validation Methylation Data from rnb.validation
# ============================================================================
beta_matrix <- rnb.validation@meth.sites
cpg_ids_raw <- rownames(rnb.validation@sites)
cpg_ids     <- sub("_.*", "", cpg_ids_raw)

# Match and subset CpGs
idx           <- match(cg_list, cpg_ids)
valid         <- !is.na(idx)
idx           <- idx[valid]
selected_cpgs <- cg_list[valid]

meth_table_val <- beta_matrix[idx, ]
rownames(meth_table_val) <- selected_cpgs
colnames(meth_table_val) <- rnb.validation@pheno$txt_PATHOLOGIE_NUMMER
dim(meth_table_val)
meth_table_val <- as.data.frame(t(meth_table_val))
anno_plot_validation$ID <- rownames(anno_plot_validation)
meth_table_val$ID <- rownames(meth_table_val)
meth_data_val_ml <- left_join(anno_plot_validation, meth_table_val, by = "ID")

# Export to CSV
write.csv(meth_data_val_ml,
          ".../data/NEN_ml_xgb_val_update.csv",
          row.names = FALSE)

# ============================================================================
# 1. Load Libraries
# ============================================================================
library(tidyverse)
library(RnBeads)
library(DecompPipeline)
library(pheatmap)

# ============================================================================
# 2. Load Data & Annotations
# ============================================================================
# Read proportions & annotation for the validation cohort
#proportions_val_LMC9  <- readRDS(".../data/validation/proportions_val_LMC9.RDS")
proportions_val_LMC10 <- readRDS(".../data/validation/proportions_val_LMC10.RDS")
annotation            <- read.csv(".../data/Clean_up_annotation_updated2025.csv", header = TRUE)
md.res_val_12k        <- readRDS(".../data/md.res_12k_val_012025.rds")

# ============================================================================
# 3. Compute & Plot Proportions (Discovery Data)
# ============================================================================
prop <- getProportions(md.res_val_12k, K = 10, lambda = 0.00001)
colnames(prop) <- annotation$ID
prop <- as.data.frame(prop)

# Subset & order annotation for plotting
anno_plot <- annotation %>% 
  select(Localization, Primary, NEN.type, P_grouping)
rownames(anno_plot) <- annotation$ID
prop   <- prop[, order(anno_plot$P_grouping)]
anno_plot <- anno_plot[order(anno_plot$P_grouping), ]

# Plot heatmap
pheatmap(prop,
         annotation = anno_plot["P_grouping", drop = FALSE],
         cluster_cols = FALSE,
         cluster_rows = FALSE)
dev.off()  # Close graphics device

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
rownames(anno_plot_validation) <- colnames(proportions_val_LMC10)

# Combine discovery and validation proportions/annotations
prop_val     <- cbind(prop, proportions_val_LMC10)
anno_plot_val <- rbind(anno_plot, anno_plot_validation)

# ============================================================================
# 5. Read 5000 CpGs & Extract Relevant IDs
# ============================================================================
df_5000_cpgs <- read_csv(".../data/NEN_ml_5k_res_test_update.csv",
                         show_col_types = FALSE)
cg_list <- colnames(df_5000_cpgs)[grepl("^cg", colnames(df_5000_cpgs))]
cg_vector <- as.vector(cg_list)  # (optional)

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
str(meth_table_val)

# ============================================================================
# 7. Prepare Proportion Data for Plotting
# ============================================================================
# Subset columns (199:214) from combined proportions and transpose
prop_val_phe <- as.data.frame(prop_val[, 199:214])
prop_t <- as.data.frame(t(prop_val_phe))

# Reassign row/column names
rownames(prop_t) <- colnames(prop_val_phe)
colnames(prop_t) <- rownames(prop_val_phe)

# Add ID for merging
prop_t$ID <- rownames(prop_t)
anno_plot_validation$ID <- rownames(anno_plot_validation)
df_plotting <- inner_join(anno_plot_validation, prop_t, by = "ID")

# ============================================================================
# 8. Calculate Residuals for Methylation Data
# ============================================================================
# Convert methylation data to data frame
meth.data <- as.data.frame(meth_table_val)

# Define a function to calculate residuals with optional adjustment
residuals_function <- function(meth, pheno, adjustment_columns = c()){
  apply(meth, 1, function(x) {
    pheno$x <- x
    form <- as.formula(paste0("x ~ ", paste(adjustment_columns, collapse = " + ")))
    mod <- lm(form, data = pheno)
    pval <- pf(summary(mod)$fstatistic[1],
               summary(mod)$fstatistic[2],
               summary(mod)$fstatistic[3],
               lower.tail = FALSE)
    if (pval < 0.01) {
      stats::residuals(mod)
    } else {
      scale(x, scale = FALSE)
    }
  })
}

# Adjust for LMC3 (assumed available in prop_t)
LMC3 <- prop_t$V3
reka_test <- residuals_function(meth = meth.data, pheno = df_plotting, adjustment_columns = c("LMC3"))
reka_test <- as.data.frame(reka_test)

# ============================================================================
# 9. Merge Residuals with Validation Annotation & Export Results
# ============================================================================
# Add index columns for merging
reka_test$index <- 1:nrow(reka_test)
anno_plot_validation$index <- 1:nrow(anno_plot_validation)

meth_data_val_ml <- left_join(anno_plot_validation, reka_test, by = "index")
# Remove columns with any NA values
meth_data_val_ml_nona <- meth_data_val_ml[, colSums(is.na(meth_data_val_ml)) == 0]

# Export to CSV
write.csv(meth_data_val_ml_nona,
          ".../data/NEN_ml_5k_res_val_update.csv",
          row.names = FALSE)

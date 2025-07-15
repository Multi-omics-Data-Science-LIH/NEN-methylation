library(plotly)
library(tidyverse)
library(RnBeads)
library(DecompPipeline)
library(RColorBrewer)
library(MeDeCom)
library(ComplexHeatmap)
library(pheatmap)
library(ggcorrplot)

my_cols <- c("Merkel cell carcinoma"="#00008B",
             "NEN appendix"="#003D00",
             "NEN colorectal"	= "#00FFFF",
             "NEN gastroduodenal"	= "#00E300",
             "NEN ileum"	= "#007D00",
             "NEN liver CUP"	= "#FF0000",
             "NEN liver metastasis"= "#FFA500",
             "NEN pancreas"= "#FF00FF",
             "Pulmonal NEC"= "#0D8CFF",
             "Pulmonary carcinoid"= "#0000FF")

# 0. Medcom 12k+Val

md.res_val_12k <- readRDS(file=".../data/md.res_12k_val_012025.rds")

annotation <- read.csv(".../data/Clean_up_annotation_updated2025.csv", header = TRUE)

# 1. Get proportions and assign sample IDs as column names
prop <- getProportions(md.res_val_12k, K = 10, lambda = 0.00001)
colnames(prop) <- annotation$ID

prop <- as.data.frame(prop)

anno_plot <- annotation[,c("Localization", "Primary", "NEN.type", "P_grouping")]

# 2. Create an annotation data frame using the P_grouping column.
rownames(anno_plot) <- annotation$ID
prop <- prop[, order(anno_plot$P_grouping)]
anno_plot <- anno_plot[order(anno_plot$P_grouping),]

# 3. Set up annotation colors
anno_colors <- list(P_grouping = my_cols)

# 4. Save the heatmap as a PDF file with specified dimensions
pdf(".../results/P_grouping_LMC10.pdf",
    width = 10, height = 5)
pheatmap(prop,
         show_colnames = FALSE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = FALSE,
         annotation_col = anno_plot[, "P_grouping", drop = FALSE],
         annotation_colors = anno_colors)
dev.off()
####################################
# Proportion heatmap (LMC2 removed)

# 1a. Label each row as LMC1, LMC2, … LMC10
rownames(prop) <- paste0("LMC", seq_len(nrow(prop)))

# 1b. Remove LMC2
prop <- prop[ !rownames(prop) %in% "LMC2", ]

# 2. Prepare annotation
anno_plot <- annotation[, c("Localization", "Primary", "NEN.type", "P_grouping")]
rownames(anno_plot) <- annotation$ID

# 2a. Order samples by P_grouping
prop <- prop[, order(anno_plot$P_grouping)]
anno_plot <- anno_plot[order(anno_plot$P_grouping), ]

# 3. Set up annotation colors
anno_colors <- list(P_grouping = my_cols)


mat <- as.matrix(prop)

annotation_col <- anno_plot[ colnames(mat), "P_grouping", drop = FALSE ]
# (this gives a one-column data.frame, rownames = sample IDs)
heat_colors <- colorRampPalette( brewer.pal(9, "YlOrRd") )(100)

# Now open your PDF and plot
pdf(".../results/P_grouping_LMC10_non_lmc2.pdf",
    width = 10, height = 5)
pheatmap(
  mat,
  color = heat_colors,
  show_colnames     = FALSE,
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,
  border_color      = NA,
  annotation_col    = annotation_col,
  annotation_colors = anno_colors,
  main              = "Proportion heatmap (LMC2 removed)"
)
dev.off()

####################################

# Get proportion for LMC-9

# 1. Get proportions and assign sample IDs as column names
prop_lmc9 <- getProportions(md.res_val_12k, K = 9, lambda = 0.00001)
colnames(prop_lmc9) <- annotation$ID

prop_lmc9 <- as.data.frame(prop_lmc9)

anno_plot <- annotation[,c("Localization", "Primary", "NEN.type", "P_grouping")]
# 2. Create an annotation data frame using the P_grouping column.
rownames(anno_plot) <- annotation$ID
prop_lmc9 <- prop_lmc9[, order(anno_plot$P_grouping)]
anno_plot <- anno_plot[order(anno_plot$P_grouping),]

# 3. Set up annotation colors
anno_colors <- list(P_grouping = my_cols)

# 4. Save the heatmap as a PDF file with specified dimensions
pdf(".../results/P_grouping_LMC9.pdf",
    width = 10, height = 5)
pheatmap(prop_lmc9,
         show_colnames = FALSE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = FALSE,
         annotation_col = anno_plot[, "P_grouping", drop = FALSE],
         annotation_colors = anno_colors)
dev.off()

######################################

## boxplots LMC10
LMCs_P_groups <- cbind(anno_plot, t(prop))

for (var in c(paste0("LMC", c(1,2:10)))){
  
  
  LMCs_P_groups$LMC <- LMCs_P_groups[,var]
  g <- ggplot(LMCs_P_groups) + geom_boxplot(aes(x=P_grouping, y=LMC, color=P_grouping), outlier.shape = NA, position = position_nudge(x=-0.1), width=0.3, lwd=1.2)+
    geom_point(aes(x=P_grouping, y=LMC, color=P_grouping), position = position_nudge(x=0.2),size = 2)+
    theme_minimal()+
    scale_color_manual(values = my_cols)+
    theme(axis.title = element_blank(),
          #panel.grid.major.x = element_blank(),
          panel.border = element_rect(colour = "darkgrey", fill=NA),
          panel.grid.minor = element_blank(),
          legend.background = element_blank(),
          legend.direction="horizontal",
          legend.position = "top",
          plot.title = element_text(size = 15),
          plot.subtitle = element_text(size = 10),
          axis.text.y = element_text(size = (14)),
          axis.text.x = element_blank(),
          plot.caption = element_text(size = 8, color = "grey70", hjust = 0))
  ggsave(paste0(".../results/lmc10/", var, "_boxplot.pdf"), width = 9, height = 4)}

##################

## boxplots LMC9
LMCs_P_groups <- cbind(anno_plot, t(prop_lmc9 ))

for (var in c(paste0("LMC", c(1,2:10)))){
  
  
  LMCs_P_groups$LMC <- LMCs_P_groups[,var]
  g <- ggplot(LMCs_P_groups) + geom_boxplot(aes(x=P_grouping, y=LMC, color=P_grouping), outlier.shape = NA, position = position_nudge(x=-0.1), width=0.3, lwd=1.2)+
    geom_point(aes(x=P_grouping, y=LMC, color=P_grouping), position = position_nudge(x=0.2),size = 2)+
    theme_minimal()+
    scale_color_manual(values = my_cols)+
    theme(axis.title = element_blank(),
          #panel.grid.major.x = element_blank(),
          panel.border = element_rect(colour = "darkgrey", fill=NA),
          panel.grid.minor = element_blank(),
          legend.background = element_blank(),
          legend.direction="horizontal",
          legend.position = "top",
          plot.title = element_text(size = 15),
          plot.subtitle = element_text(size = 10),
          axis.text.y = element_text(size = (14)),
          axis.text.x = element_blank(),
          plot.caption = element_text(size = 8, color = "grey70", hjust = 0))
  ggsave(paste0(".../results/", var, "_boxplot.pdf"), width = 9, height = 4)}

###########################################

# Calculate lump for Medcom 12k+Val
rnb.original_processed_updated_anno_012025_hg19 <- readRDS(".../data/rnb.original_processed_updated_anno_012025_hg19.RDS")
lump <- rnb.execute.lump(rnb.original_processed_updated_anno_012025_hg19)

### correlation between LUMP and LMC10
# Get proportion for LMC-10

prop <- getProportions(md.res_val_12k, K=10, lambda=0.00001)
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
df_plotting$lump <- lump


corrdf <- cor(df_plotting[13:23])
p2 <- ggcorrplot(corrdf, type = "lower",
                 lab = TRUE,  method = "circle")

ggsave(
  filename = paste0(".../results/test/lump_lmc10.pdf"),
  plot = p2,
  width = 9,
  height = 4
)


### correlation between LUMP and LMC9

# Get proportion for LMC-9

prop_lmc9 <- getProportions(md.res_val_12k, K=9, lambda=0.00001)
colnames(prop_lmc9) <- annotation$ID

prop_lmc9 <- as.data.frame(prop_lmc9)
prop_t <- t(prop_lmc9)
prop_t <- as.data.frame(prop_t)

#redefine row and column names
rownames(prop_t) <- colnames(prop_lmc9)
colnames(prop_t) <- rownames(prop_lmc9)

# re-assign an index column for left_join

prop_t$ID <- as.integer(rownames(prop_t))

df_plotting <- inner_join(annotation, prop_t, by = "ID")
df_plotting$lump <- lump


corrdf <- cor(df_plotting[13:22])
p9 <- ggcorrplot(corrdf, type = "lower",
                 lab = TRUE,  method = "circle")

ggsave(
  filename = paste0(".../results/test/lump_lmc9.pdf"),
  plot = p9,
  width = 9,
  height = 4
)

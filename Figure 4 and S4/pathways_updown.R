# Load necessary libraries
source("functions.R")
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(factoextra)
library(clusterProfiler)
library(org.Mm.eg.db)
library(pathview)
library(Rtsne)
library(pheatmap)
library(RColorBrewer)
library(viridis)

set.seed(1234)
# Read the count data
counts_data <- read.csv("bigomics/counts.csv")

# Creating the metadata
sample_names <- colnames(counts_data)[-1]  # exclude the gene ID column
conditions <- c(rep("Ctrl", 5), rep("Ctrl_IGF1", 3), rep("LPHC_IGF1", 5), rep("LPHC", 5))
coldata <- data.frame(row.names=sample_names, condition=conditions)

# Renaming the levels
coldata$condition <- factor(coldata$condition, levels = c("Ctrl", "Ctrl_IGF1", "LPHC_IGF1", "LPHC"))
levels(coldata$condition) <- c("Ctrl", "Ctrl_IGF1", "LPHC_IGF1", "LPHC")

# Create the DESeqDataSet object again
dds <- DESeqDataSetFromMatrix(countData = counts_data[,-1],  # Exclude the gene ID column
                              colData = coldata,
                              design = ~ condition) #general groups, no interesection IGF1 / diet

# Normalize the data
dds <- DESeq(dds)
res <- results(dds)

## DEGS ####
# Differential expression for Standard Diet with IGF1 vs. Standard Diet without IGF1
res_ctrl <- results(dds, contrast=c("condition", "Ctrl_IGF1", "Ctrl"))
rownames(res_ctrl) <- counts_data$X

# Differential expression for LPHC Diet with IGF1 vs. LPHC Diet without IGF1
res_lphc <- results(dds, contrast=c("condition", "LPHC_IGF1", "LPHC"))
rownames(res_lphc) <- counts_data$X


## From Deepseek Coder
# Filter significant genes
padj_threshold <- 0.05
log2fc_threshold <- 1

# For Ctrl_IGF1 vs Ctrl
sig_genes_ctrl <- res_ctrl %>%
  as.data.frame() %>%
  filter(padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold) %>%
  rownames_to_column(var = "gene_id")

# save
sig_genes_ctrl$symbol <- mapIds(org.Mm.eg.db, keys=sig_genes_ctrl$gene_id, keytype="ENSEMBL", column="SYMBOL")
write.csv(sig_genes_ctrl, file="review/sig_degs_Ctrl_IGF1_vs_Ctrl.csv", row.names=TRUE)

# For LPHC_IGF1 vs LPHC
sig_genes_lphc <- res_lphc %>%
  as.data.frame() %>%
  filter(padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold) %>%
  rownames_to_column(var = "gene_id")

# save
sig_genes_lphc$symbol <- mapIds(org.Mm.eg.db, keys=sig_genes_lphc$gene_id, keytype="ENSEMBL", column="SYMBOL")
write.csv(sig_genes_lphc, file="review/sig_degs_LPHC_IGF1-vs_LPHC.csv", row.names=TRUE)


# Split up-regulated and down-regulated genes
# For up-regulated genes in Ctrl_IGF1 vs Ctrl
up_genes_ctrl <- sig_genes_ctrl %>% filter(log2FoldChange > 0)
up_genes_ctrl_entrez <- mapIds(org.Mm.eg.db, keys=up_genes_ctrl$gene_id, keytype="ENSEMBL", column="ENTREZID")
up_genes_ctrl_entrez <- up_genes_ctrl_entrez[!is.na(up_genes_ctrl_entrez)]

# For down-regulated genes in Ctrl_IGF1 vs Ctrl
down_genes_ctrl <- sig_genes_ctrl %>% filter(log2FoldChange < 0)
down_genes_ctrl_entrez <- mapIds(org.Mm.eg.db, keys=down_genes_ctrl$gene_id, keytype="ENSEMBL", column="ENTREZID")
down_genes_ctrl_entrez <- down_genes_ctrl_entrez[!is.na(down_genes_ctrl_entrez)]

# For up-regulated genes in LPHC_IGF1 vs LPHC
up_genes_lphc <- sig_genes_lphc %>% filter(log2FoldChange > 0)
up_genes_lphc_entrez <- mapIds(org.Mm.eg.db, keys=up_genes_lphc$gene_id, keytype="ENSEMBL", column="ENTREZID")
up_genes_lphc_entrez <- up_genes_lphc_entrez[!is.na(up_genes_lphc_entrez)]

# For down-regulated genes in LPHC_IGF1 vs LPHC
down_genes_lphc <- sig_genes_lphc %>% filter(log2FoldChange < 0)
down_genes_lphc_entrez <- mapIds(org.Mm.eg.db, keys=down_genes_lphc$gene_id, keytype="ENSEMBL", column="ENTREZID")
down_genes_lphc_entrez <- down_genes_lphc_entrez[!is.na(down_genes_lphc_entrez)]




##### Perform enrichGO GSEA analysis ####
# For up-regulated genes in Ctrl_IGF1 vs Ctrl
ego_up_ctrl <- enrichGO(gene = up_genes_ctrl_entrez,
                        OrgDb = org.Mm.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)
ego_up_ctrl_smpld <- simplify(ego_up_ctrl, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")


# For down-regulated genes in Ctrl_IGF1 vs Ctrl
ego_down_ctrl <- enrichGO(gene = down_genes_ctrl_entrez,
                          OrgDb = org.Mm.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)
ego_down_ctrl_smpld <- simplify(ego_down_ctrl, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")

# For up-regulated genes in LPHC_IGF1 vs LPHC
ego_up_lphc <- enrichGO(gene = up_genes_lphc_entrez,
                        OrgDb = org.Mm.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)
ego_up_lphc_smpld <- simplify(ego_up_lphc, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")

# For down-regulated genes in LPHC_IGF1 vs LPHC
ego_down_lphc <- enrichGO(gene = down_genes_lphc_entrez,
                          OrgDb = org.Mm.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)
ego_down_lphc_smpld <- simplify(ego_down_lphc, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")

# Write the results of the up-regulated genes in Ctrl_IGF1 vs Ctrl to a CSV file
write.csv(ego_up_ctrl, file="review/ego_up_ctrl.csv", row.names=TRUE)
write.csv(ego_up_ctrl_smpld, file="review/ego_up_ctrl_smpld.csv", row.names=TRUE)

# Write the results of the down-regulated genes in Ctrl_IGF1 vs Ctrl to a CSV file
write.csv(ego_down_ctrl, file="review/ego_down_ctrl.csv", row.names=TRUE)
write.csv(ego_down_ctrl_smpld, file="review/ego_down_ctrl_smpld.csv", row.names=TRUE)

# Write the results of the up-regulated genes in LPHC_IGF1 vs LPHC to a CSV file
write.csv(ego_up_lphc, file="review/ego_up_lphc.csv", row.names=TRUE)
write.csv(ego_up_lphc_smpld, file="review/ego_up_lphc_smpld.csv", row.names=TRUE)

# Write the results of the down-regulated genes in LPHC_IGF1 vs LPHC to a CSV file
write.csv(ego_down_lphc, file="review/ego_down_lphc.csv", row.names=TRUE)
write.csv(ego_down_lphc_smpld, file="review/ego_down_lphc_smpld.csv", row.names=TRUE)


# Visualization for up-regulated genes in Ctrl_IGF1 vs Ctrl
plot_up_ctrl <- dotplot(ego_up_ctrl, showCategory = 10) + ggtitle("Up-regulated genes in Ctrl_IGF1 vs Ctrl")
ggsave(paste0("review/Plots/Up_Ctrl_IGF1_vs_Ctrl.png"), plot = plot_up_ctrl, dpi = 300, width = 6, height = 5)

# Visualization for down-regulated genes in Ctrl_IGF1 vs Ctrl
plot_down_ctrl <- dotplot(ego_down_ctrl, showCategory = 20) + ggtitle("Down-regulated genes in Ctrl_IGF1 vs Ctrl")
ggsave(paste0("review/Plots/Down_Ctrl_IGF1_vs_Ctrl.png"), plot = plot_down_ctrl, dpi = 300, width = 6, height = 5)

# Visualization for up-regulated genes in LPHC_IGF1 vs LPHC
plot_up_lphc <- dotplot(ego_up_lphc, showCategory = 10) + ggtitle("Up-regulated genes in LPHC_IGF1 vs LPHC")
ggsave(paste0("review/Plots/Up_LPHC_IGF1_vs_LPHC.png"), plot = plot_up_lphc, dpi = 300, width = 6, height = 5)

# Visualization for down-regulated genes in LPHC_IGF1 vs LPHC
plot_down_lphc <- dotplot(ego_down_lphc, showCategory = 20) + ggtitle("Down-regulated genes in LPHC_IGF1 vs LPHC")
ggsave(paste0("review/Plots/Down_LPHC_IGF1_vs_LPHC.png"), plot = plot_down_lphc, dpi = 300, width = 6, height = 5)

# Visualization for up-regulated genes in Ctrl_IGF1 vs Ctrl
plot_up_ctrl_smpld <- dotplot(ego_up_ctrl_smpld, showCategory = 10) + ggtitle("Up-regulated genes in Ctrl_IGF1 vs Ctrl (Simplified)")
ggsave(paste0("review/Plots/Up_Ctrl_IGF1_vs_Ctrl_smpld.png"), plot = plot_up_ctrl_smpld, dpi = 300, width = 6, height = 5)

# Visualization for down-regulated genes in Ctrl_IGF1 vs Ctrl
plot_down_ctrl_smpld <- dotplot(ego_down_ctrl_smpld, showCategory = 10) + ggtitle("Down-regulated genes in Ctrl_IGF1 vs Ctrl (Simplified)")
ggsave(paste0("review/Plots/Down_Ctrl_IGF1_vs_Ctrl_smpld.png"), plot = plot_down_ctrl_smpld, dpi = 300, width = 6, height = 5)

# Visualization for up-regulated genes in LPHC_IGF1 vs LPHC
plot_up_lphc_smpld <- dotplot(ego_up_lphc_smpld, showCategory = 10) + ggtitle("Up-regulated genes in LPHC_IGF1 vs LPHC (Simplified)")
ggsave(paste0("review/Plots/Up_LPHC_IGF1_vs_LPHC_smpld.png"), plot = plot_up_lphc_smpld, dpi = 300, width = 6, height = 5)

# Visualization for down-regulated genes in LPHC_IGF1 vs LPHC
plot_down_lphc_smpld <- dotplot(ego_down_lphc_smpld, showCategory = 10) + ggtitle("Down-regulated genes in LPHC_IGF1 vs LPHC (Simplified)")
ggsave(paste0("review/Plots/Down_LPHC_IGF1_vs_LPHC_smpld.png"), plot = plot_down_lphc_smpld, dpi = 300, width = 6, height = 5)

# KEGG enrichment analysis for up-regulated genes in Ctrl_IGF1 vs Ctrl
kegg_up_ctrl <- enrichKEGG(gene = up_genes_ctrl_entrez, 
                           organism = 'mmu', # Mouse
                           pvalueCutoff = 0.05)

# Write the results of the KEGG enrichment analysis for up-regulated genes in Ctrl_IGF1 vs Ctrl to a CSV file
write.csv(kegg_up_ctrl, file="review/kegg_up_ctrl.csv", row.names=TRUE)

# Visualization for KEGG enrichment analysis for up-regulated genes in Ctrl_IGF1 vs Ctrl
plot_kegg_up_ctrl <- dotplot(kegg_up_ctrl, showCategory = 10) + ggtitle("KEGG Enrichment - Up-regulated genes in Ctrl_IGF1 vs Ctrl")
ggsave(paste0("review/Plots/KEGG_Up_Ctrl_IGF1_vs_Ctrl.png"), plot = plot_kegg_up_ctrl, dpi = 300, width = 6, height = 5)

# KEGG enrichment analysis for down-regulated genes in Ctrl_IGF1 vs Ctrl
kegg_down_ctrl <- enrichKEGG(gene = down_genes_ctrl_entrez, 
                             organism = 'mmu', # Mouse
                             pvalueCutoff = 0.05)


# Write the results of the KEGG enrichment analysis for down-regulated genes in Ctrl_IGF1 vs Ctrl to a CSV file
write.csv(kegg_down_ctrl, file="review/kegg_down_ctrl.csv", row.names=TRUE)

# Visualization for KEGG enrichment analysis for down-regulated genes in Ctrl_IGF1 vs Ctrl
plot_kegg_down_ctrl <- dotplot(kegg_down_ctrl, showCategory = 10) + ggtitle("KEGG Enrichment - Down-regulated genes in Ctrl_IGF1 vs Ctrl")
ggsave(paste0("review/Plots/KEGG_Down_Ctrl_IGF1_vs_Ctrl.png"), plot = plot_kegg_down_ctrl, dpi = 300, width = 6, height = 5)

# KEGG enrichment analysis for up-regulated genes in LPHC_IGF1 vs LPHC
kegg_up_lphc <- enrichKEGG(gene = up_genes_lphc_entrez, 
                           organism = 'mmu', # Mouse
                           pvalueCutoff = 0.05)

# Write the results of the KEGG enrichment analysis for up-regulated genes in LPHC_IGF1 vs LPHC to a CSV file
write.csv(kegg_up_lphc, file="review/kegg_up_lphc.csv", row.names=TRUE)

# Visualization for KEGG enrichment analysis for up-regulated genes in LPHC_IGF1 vs LPHC
plot_kegg_up_lphc <- dotplot(kegg_up_lphc, showCategory = 10) + ggtitle("KEGG Enrichment - Up-regulated genes in LPHC_IGF1 vs LPHC")
ggsave(paste0("review/Plots/KEGG_Up_LPHC_IGF1_vs_LPHC.png"), plot = plot_kegg_up_lphc, dpi = 300, width = 6, height = 5)

# KEGG enrichment analysis for down-regulated genes in LPHC_IGF1 vs LPHC
kegg_down_lphc <- enrichKEGG(gene = down_genes_lphc_entrez, 
                             organism = 'mmu', # Mouse
                             pvalueCutoff = 0.05)

# Write the results of the KEGG enrichment analysis for down-regulated genes in LPHC_IGF1 vs LPHC to a CSV file
write.csv(kegg_down_lphc, file="review/kegg_down_lphc.csv", row.names=TRUE)

# Visualization for KEGG enrichment analysis for down-regulated genes in LPHC_IGF1 vs LPHC
plot_kegg_down_lphc <- dotplot(kegg_down_lphc, showCategory = 10) + ggtitle("KEGG Enrichment - Down-regulated genes in LPHC_IGF1 vs LPHC")
ggsave(paste0("review/Plots/KEGG_Down_LPHC_IGF1_vs_LPHC.png"), plot = plot_kegg_down_lphc, dpi = 300, width = 6, height = 5)


#### Common gene Ctrl and LPHC anaylsis ####
up_common_genes <- intersect(up_genes_ctrl_entrez,up_genes_lphc_entrez)
down_common_genes <- intersect(down_genes_ctrl_entrez,down_genes_lphc_entrez)

# EnrichGO for up-regulated common genes
ego_up_common <- enrichGO(gene = up_common_genes,
                          OrgDb = org.Mm.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

# Simplify the GO enrichment results for up-regulated common genes
ego_up_common_smpld <- simplify(ego_up_common, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")

# Write the results of the GO enrichment analysis for up-regulated common genes to a CSV file
write.csv(ego_up_common, file="ego_up_common.csv", row.names=TRUE)
write.csv(ego_up_common_smpld, file="ego_up_common_smpld.csv", row.names=TRUE)

# Visualization for up-regulated common genes (non-simplified)
plot_up_common <- dotplot(ego_up_common, showCategory = 10) + ggtitle("Up-regulated common genes between Ctrl and LPHC")
ggsave(paste0("review/Plots/Up_Common_Genes.png"), plot = plot_up_common, dpi = 300, width = 6, height = 5)

# Visualization for up-regulated common genes (simplified)
plot_up_common_smpld <- dotplot(ego_up_common_smpld, showCategory = 10) + ggtitle("Up-regulated common genes between Ctrl and LPHC (Simplified)")
ggsave(paste0("review/Plots/Up_Common_Genes_smpld.png"), plot = plot_up_common_smpld, dpi = 300, width = 6, height = 5)

# EnrichGO for down-regulated common genes
ego_down_common <- enrichGO(gene = down_common_genes,
                            OrgDb = org.Mm.eg.db,
                            keyType = "ENTREZID",
                            ont = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)

# Simplify the GO enrichment results for down-regulated common genes
ego_down_common_smpld <- simplify(ego_down_common, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")

# Write the results of the GO enrichment analysis for down-regulated common genes to a CSV file
write.csv(ego_down_common, file="ego_down_common.csv", row.names=TRUE)
write.csv(ego_down_common_smpld, file="ego_down_common_smpld.csv", row.names=TRUE)

# Visualization for down-regulated common genes (non-simplified)
plot_down_common <- dotplot(ego_down_common, showCategory = 10) + ggtitle("Down-regulated common genes between Ctrl and LPHC")
ggsave(paste0("review/Plots/Down_Common_Genes.png"), plot = plot_down_common, dpi = 300, width = 6, height = 5)

# Visualization for down-regulated common genes (simplified)
plot_down_common_smpld <- dotplot(ego_down_common_smpld, showCategory = 10) + ggtitle("Down-regulated common genes between Ctrl and LPHC (Simplified)")
ggsave(paste0("review/Plots/Down_Common_Genes_smpld.png"), plot = plot_down_common_smpld, dpi = 300, width = 6, height = 5)


##### EnrichKEGG for up-regulated common genes ####
kegg_up_common <- enrichKEGG(gene = up_common_genes,
                             organism = 'mmu', # Mouse
                             pvalueCutoff = 0.05)

# Write the results of the KEGG enrichment analysis for up-regulated common genes to a CSV file
write.csv(kegg_up_common, file="kegg_up_common.csv", row.names=TRUE)

# Visualization for up-regulated common genes (non-simplified)
plot_kegg_up_common <- dotplot(kegg_up_common, showCategory = 10) + ggtitle("KEGG Enrichment - Up-regulated common genes between Ctrl and LPHC")
ggsave(paste0("review/Plots/KEGG_Up_Common_Genes.png"), plot = plot_kegg_up_common, dpi = 300, width = 6, height = 5)

# EnrichKEGG for down-regulated common genes
kegg_down_common <- enrichKEGG(gene = down_common_genes,
                               organism = 'mmu', # Mouse
                               pvalueCutoff = 0.05)

# Write the results of the KEGG enrichment analysis for down-regulated common genes to a CSV file
write.csv(kegg_down_common, file="kegg_down_common.csv", row.names=TRUE)

# Visualization for down-regulated common genes (non-simplified)
plot_kegg_down_common <- dotplot(kegg_down_common, showCategory = 10) + ggtitle("KEGG Enrichment - Down-regulated common genes between Ctrl and LPHC")
ggsave(paste0("review/Plots/KEGG_Down_Common_Genes.png"), plot = plot_kegg_down_common, dpi = 300, width = 6, height = 5)


##### heatmap common genes #####
## z-score for heatmap DEG gene in modules

# Renaming the levels
coldata$condition <- factor(coldata$condition, levels = c("Ctrl", "LPHC", "Ctrl_IGF1", "LPHC_IGF1"))
levels(coldata$condition) <- c("Ctrl", "LPHC", "Ctrl_IGF1", "LPHC_IGF1")

# Filter out lowly expressed genes
keep_genes <- rowSums(counts_data > 10) > (0.5 * ncol(counts_data))
expr_matrix_filtered <- counts_data[keep_genes, ]
rownames(expr_matrix_filtered) <- expr_matrix_filtered$X

# Normalization using DESeq2 (as already done)
dds <- DESeqDataSetFromMatrix(countData = expr_matrix_filtered[,-1], colData = coldata, design = ~ condition)
dds <- estimateSizeFactors(dds)
vst_transform <- vst(dds, blind=FALSE)  # Directly get VST transformation
datExpr <- as.data.frame(assay(vst_transform))
datExpr$gene <- rownames(datExpr)

common_genes <- data.frame(entrez= c(up_common_genes,down_common_genes))
common_genes$ensembl <- mapIds(org.Mm.eg.db, list_genes$entrez, keytype = "ENTREZID", column = "ENSEMBL")

exp_common_genes <- datExpr %>% dplyr::filter(gene %in% common_genes$ensembl)

exp_common_genes$Symbol <- mapIds(org.Mm.eg.db, rownames(exp_common_genes), keytype = "ENSEMBL", column = "SYMBOL")
exp_common_genes <- na.omit(exp_common_genes)
rownames(exp_common_genes) <- exp_common_genes$Symbol
exp_common_genes <- subset(exp_common_genes, select = - c(gene,Symbol))

all.group <- as.vector(coldata$condition)

#Declare z-score function
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# Function to compute percentile transformation
percentile_transformation <- function(logCPM_matrix) {
  
  # Calculate the percentile rank for each value in the matrix
  percentile_matrix <- t(apply(logCPM_matrix, 1, function(x) {
    ranks <- rank(x, na.last = "keep", ties.method = "average")
    percentiles <- (ranks - 1) / (length(x) - 1)
    return(percentiles)
  }))
  
  return(percentile_matrix)
}

#calculate z-score
mean_common_genes <- apply(exp_common_genes, 1, function(x)by(x, all.group, mean))
mean_exp_common_genes_zscore <- t(apply(mean_common_genes, 1, cal_z_score))

# Perform percentile transformation
percentile_matrix <- percentile_transformation(exp_common_genes)
mean_percentile_matrix <- apply(percentile_matrix, 1, function(x)by(x, all.group, mean))
mean_exp_common_genes_zscore <- mean_percentile_matrix

# Order the matrix based on the custom order vector
custom_order <- c("Ctrl", "LPHC", "Ctrl_IGF1", "LPHC_IGF1")
mean_common_genes <- mean_common_genes[match(custom_order, rownames(mean_common_genes)), ]
mean_exp_common_genes_zscore <- mean_exp_common_genes_zscore[match(custom_order, rownames(mean_exp_common_genes_zscore)), ]

# heatmap of gene expression from the selected module
group <- data.frame(DEGs = c(rep("Upregulated", length(up_common_genes)),rep("Downregulated",length(down_common_genes))), gene = c(up_common_genes,down_common_genes))
group$gene <- mapIds(org.Mm.eg.db, group$gene, keytype = "ENTREZID", column = "SYMBOL")
group <- group %>% dplyr::filter(gene %in% colnames(mean_exp_common_genes_zscore))

# select top 20 of each
# not ordered 
top20_group <- data.frame(rbind(head(group[group$DEGs == "Upregulated",],20),head(group[group$DEGs == "Downregulated",],20)))

# top slopes
# Calculate the slopes for each gene
# Function to calculate the slope (coefficient) from linear regression
calculate_slope <- function(gene_exp) {
  model <- lm(gene_exp ~ seq_along(gene_exp))
  return(coef(model)[2])  # Return the slope (second coefficient)
}

# Calculate the slopes for each gene
slopes <- apply(mean_common_genes, 2, calculate_slope)
# Convert slopes to a named vector
# Convert slopes to a data frame
slopes_df <- data.frame(Slope = slopes, Gene = colnames(mean_common_genes))

# Extract the top 20 positive slopes and top 20 negative slopes
top_positive_slopes <- head(slopes_df[order(slopes_df$Slope, decreasing = TRUE), ], 20)
top_negative_slopes <- head(slopes_df[order(slopes_df$Slope, decreasing = FALSE), ], 20)
# Combine the top slopes into a single data frame
top_slopes_df <- rbind(top_positive_slopes, top_negative_slopes)
top_slopes_df$group <- c(rep("Upregulated",20),rep("Downregulated",20))

#cleanup groups
row.names(group) <- group[,2]
colnames(group)[1] <- "group"
group <- subset(group, select = - 2)

row.names(top20_group) <- top20_group[,2]
colnames(top20_group)[1] <- "group"
top20_group <- subset(top20_group, select = - 2)

row.names(top_slopes_df) <- top_slopes_df[,2]
colnames(top_slopes_df)[3] <- "group"
top_slopes_df <- subset(top_slopes_df, select = - c(1:2))


my_colour = list("group" = 
                   c("Upregulated" = "red", 
                     "Downregulated" = "darkblue"))

pdf("review/Plots/Heatmap_top20_updown_deg_common_genes_p_transform.pdf",onefile=T,width=13,height=6)
mybreaks <- c(
  #seq(min(mean_kid_modGenes), -0.01, length.out=50), #_zscore
  #seq(0, max(mean_kid_modGenes),length.out=50)
  seq(0, 1,length.out=100)
) 
color <- viridis(100) #colorRampPalette(c("blue","white","red"))(100)
pheatmap(mean_exp_common_genes_zscore[, colnames(mean_exp_common_genes_zscore) %in% rownames(top_slopes_df)], #mean_exp_common_genes_zscore, #mean_exp_common_genes_zscore[, colnames(mean_exp_common_genes_zscore) %in% rownames(top20_group)],
         cluster_rows = F,
         cluster_cols = T,
         #show_colnames = TRUE,
         #show_rownames = TRUE,
         main = "Percentile-score of mean gene expression",
         color = color,
         #breaks = mybreaks,
         #cutree_col = 4, 
         annotation_col = top_slopes_df, #top20_group, #group, 
         annotation_colors = my_colour, 
         annotation_legend = TRUE,
         border_color = NA,
         cellheight = 15,
         cellwidth = 12,
         fontsize_col = 12
         #filename = "heatmap.pdf"
)
dev.off()


#####################################################
## Intersection analysis with diet and IGF1
diets <- c(rep("Standard", 5), rep("Standard", 3), rep("LowProteinHighCarbs", 5), rep("LowProteinHighCarbs", 5))
coldata$diet <- factor(diets)

treatment <- c(rep("Veh", 5), rep("IGF1", 3), rep("IGF1", 5), rep("Veh", 5))
coldata$treatment <- factor(treatment)

# Update design formula to include diet and interaction term
dds <- DESeqDataSetFromMatrix(countData = counts_data[,-1], colData = coldata, design = ~ diet + treatment + diet:treatment)
dds <- DESeq(dds)

# Results for interaction terms
resultsNames(dds)
res_interaction <- results(dds, name="dietStandard.treatmentVeh")

# Plotting MDS
library(limma)
library(Mfuzz)
rld <- rlog(dds)
data <- plotMDS(as.matrix(assay(rld)), col=c("red", "blue", "green", "purple")[coldata$condition])
print(data)

# Prepare the data for clustering
# Prepare the data for Mfuzz clustering
eset <- new("ExpressionSet", exprs = assay(rld))
eset <- filter.NA(eset, thres = 0.25)
eset <- fill.NA(eset, mode = "mean")
eset <- standardise(eset)


# Prepare the data for Mfuzz clustering
eset <- new("ExpressionSet", exprs = assay(rld))

# Remove rows with NA, NaN, or Inf values
eset <- eset[complete.cases(exprs(eset)) & is.finite(rowSums(exprs(eset))), ]

# Filter genes with low expression
eset <- filter.std(eset, min.std = 0)

# Fill NA values with mean
eset <- fill.NA(eset, mode = "mean")

# Standardize the data
eset <- standardise(eset)


# Perform the clustering
num.clusters <- 12  # Adjust based on your analysis
membership <- 1.25  # Fuzzification parameter, adjust based on your data
cl <- mfuzz(eset, centers = num.clusters, m = membership)

# Plotting clusters
pdf("Cluster_plots.pdf")
mfuzz.plot2(eset, cl = cl, mfrow = c(4, 3), time.labels = colnames(eset))
dev.off()

###############################################################
### test zone


## Clustering and Visualization
# Normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
# Assign gene IDs to rownames
rownames(normalized_counts) <- counts_data$X

# t-SNE for visualization
tsne_result <- Rtsne(t(normalized_counts), perplexity = 5)
tsne_df <- data.frame(tsne_result$Y, condition = coldata$condition)

ggplot(tsne_df, aes(x = X1, y = X2, color = condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "t-SNE plot of samples", x = "t-SNE dimension 1", y = "t-SNE dimension 2")

# Heatmap of top significant genes
top_genes <- unique(c(sig_genes_ctrl$gene_id, sig_genes_lphc$gene_id))
top_genes_counts <- normalized_counts[top_genes, ]

heatmap(as.matrix(top_genes_counts), col = colorRampPalette(c("blue", "white", "red"))(100),
        scale = "row", margins = c(10, 5), labRow = NA, labCol = coldata$condition)

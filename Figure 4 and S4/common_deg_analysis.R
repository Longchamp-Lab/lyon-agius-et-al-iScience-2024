# Load necessary libraries
source("functions.R")
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(pathview)
library(pheatmap)
library(edgeR)
library(viridis)

# Read the degs data from ctrl igf1 vs ctrl and lphc igf1 vs lhpc
deg_ctrl_csv <- read.csv("deg_ctrl.csv")
deg_lphc_csv <- read.csv("deg_lphc.csv")

# Identify common genes (assuming gene names are in the first column)
common_genes <- intersect(deg_ctrl_csv[,1], deg_lphc_csv[,1])

# Map ENSEMBL IDs to ENTREZ IDs
common_genes_entrez <- mapIds(org.Mm.eg.db, keys = common_genes, keytype = "ENSEMBL", column = "ENTREZID")

# Remove any NA values (genes that couldn't be mapped)
common_genes_entrez <- common_genes_entrez[!is.na(common_genes_entrez)]

# GO enrichment analysis
ego_common <- enrichGO(gene = common_genes_entrez, 
                       OrgDb = org.Mm.eg.db, 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

# View the top enriched GO terms
head(ego_common)

# KEGG enrichment analysis
kegg_common <- enrichKEGG(gene = common_genes_entrez, 
                          organism = 'mmu', 
                          pvalueCutoff = 0.05)

# View the top enriched KEGG pathways
head(kegg_common)

# Dot plot for ego_common
plot <- dotplot(ego_common, showCategory=10) + ggtitle("GO Enrichment - Common")
ggsave(paste0("Plots/GO_Common.png"), plot = plot, dpi = 300, width = 6, height = 5)

# Dot plot for kegg_common
plot <- dotplot(kegg_common, showCategory=10) + ggtitle("KEGG Common")
ggsave(paste0("Plots/KEGG_Common.png"), plot = plot, dpi = 300, width = 4.75, height = 6)

# Visualize the most significantly enriched pathway for control
pathview(gene.data = setNames(kegg_common$GeneRatio, kegg_common$Description)[1], pathway.id = kegg_common$ID[1], species = "mmu")


# Cell cycle: This pathway is essential for cell division and controls the progression of cells through DNA synthesis and mitosis. 
# The enrichment of this pathway indicates a significant role of the common DEGs in cell cycle regulation.
# 
# Progesterone-mediated oocyte maturation: This pathway is essential for the maturation of oocytes, preparing them for fertilization.
# 
# Protein digestion and absorption: This pathway is vital for the breakdown of dietary proteins into amino acids and small 
# peptides for absorption.
# 
# Renin-angiotensin system: This system is a hormone system that regulates blood pressure and fluid balance.
# It has a variety of effects, including vasoconstriction and the stimulation of aldosterone production from the adrenal glands.
# 
# ECM-receptor interaction: The extracellular matrix (ECM) interacts with cells through various receptors,
# influencing cell behavior, proliferation, and differentiation.
# Oocyte meiosis: This pathway is related to the maturation of oocytes and their progression through the 
# meiotic cell cycle.

## heatmap common genes
counts_data <- read.csv("bigomics/counts.csv")


sample_names <- colnames(counts_data)[-1] # exclude the gene ID column
conditions <- c(rep("Ctrl", 5), rep("Ctrl_IGF1", 3), rep("LPHC_IGF1", 5), rep("LPHC", 5))
coldata <- data.frame(row.names=sample_names, condition=conditions)
# Renaming the levels
coldata$condition <- factor(coldata$condition, levels = c("Ctrl", "Ctrl_IGF1", "LPHC_IGF1", "LPHC"))
levels(coldata$condition) <- c("Ctrl", "Ctrl_IGF1", "LPHC_IGF1", "LPHC")

y <- DGEList(counts=counts_data)
y <- calcNormFactors(y)
log2cpm <- cpm(y, log=TRUE, prior.count=3)
rownames(log2cpm) <- counts_data$X

common_matrix <- log2cpm[rownames(log2cpm) %in% common_genes,]
common_matrix <- as.data.frame(common_matrix)
common_matrix$symbol <- mapIds(org.Mm.eg.db, keys = rownames(common_matrix), keytype = "ENSEMBL", column = "SYMBOL")
common_matrix <- na.omit(common_matrix)
rownames(common_matrix) <- common_matrix$symbol
common_matrix <- common_matrix[, !(names(common_matrix) == "symbol")]

png(filename = "Plots/heatmap_common_genes.png", width = 5, height = 7, units = "in", res = 600)
pheatmap(common_matrix, 
         color = viridis(100), 
         annotation_col = coldata, 
         fontsize_row = 2,  # Adjust the font size to prevent overlap
         scale = "row", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean")
dev.off()


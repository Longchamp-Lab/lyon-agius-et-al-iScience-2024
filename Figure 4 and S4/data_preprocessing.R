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
                              design = ~ condition)

# Normalize the data
dds <- DESeq(dds)

## PCA ####
# Variance stabilizing transformation
vsd <- vst(dds)

# Compute PCA on VST transformed data
pca_res <- prcomp(t(assay(vsd)))

# Plot PCA with ellipses
p <- fviz_pca_ind(pca_res, 
                  geom.ind = "point", 
                  geom.ind.sup = "arrow", 
                  col.ind = coldata$condition, 
                  palette = "jco",
                  addEllipses = T,
                  ellipse.type = "confidence",
                  legend.title = "Groups",
                  repel = TRUE) +
  theme_minimal() +
  labs(title="PCA of VST-transformed data") +
  theme(legend.position="right", legend.title=element_blank(),
        legend.text=element_text(size=12), title=element_text(size=14),
        axis.title=element_text(size=12))

print(p)
ggsave(paste0("Plots/PCA_vstdata.png"), plot = p, dpi = 300, width = 6, height = 6)

# plot tsne map
tsne_res <- Rtsne(t(assay(vsd)), perplexity=5)

# Plot t-SNE
p_tsne <- ggplot(data.frame(tsne_res$Y), aes(x=X1, y=X2, color=coldata$condition)) + 
  geom_point(size=2) +
  scale_color_manual(values=c("#D55E00", "#009E73", "#0072B2", "#F0E442")) +
  labs(title="t-SNE of VST-transformed data") +
  theme_minimal() +
  theme(legend.position="right", legend.title=element_blank(),
        legend.text=element_text(size=12), title=element_text(size=14),
        axis.title=element_text(size=12))

print(p_tsne)
ggsave(paste0("Plots/tSNE_vstdata.png"), plot = p_tsne, dpi = 300, width = 6, height = 6)

## Try to create ellipse
library(ellipse)
# Compute the ellipses for each group
ellipses <- data.frame(X1=numeric(), X2=numeric(), group=character())
for(cond in unique(ellipses_data$group)) {
  el_data <- with(subset(ellipses_data, group == cond), 
                  ellipse(cor(X1, X2), centre=c(mean(X1), mean(X2)), level=0.68))
  
  # Rename columns from x and y to X1 and X2
  colnames(el_data) <- c("X1", "X2")
  
  # Convert el_data to data frame and add the "group" column
  el_data <- data.frame(el_data)
  el_data$group <- cond
  
  # Print column names for debugging
  cat("Column names for el_data:", names(el_data), "\n")
  cat("Column names for ellipses:", names(ellipses), "\n")
  
  # Ensure that the column names match before binding
  if (identical(names(ellipses), names(el_data))) {
    ellipses <- rbind(ellipses, el_data)
  } else {
    stop(paste("Column names do not match for group:", cond))
  }
}

# Plot t-SNE with ellipses
p_tsne <- ggplot(ellipses_data, aes(x=X1, y=X2, color=group)) + 
  #geom_polygon(data=ellipses, aes(fill=group), alpha=0.3, level = 0.95) + # Add the ellipses
  geom_point(size=2) +
  scale_color_manual(values=c("#D55E00", "#009E73", "darkmagenta", "darkgrey")) +
  labs(title="t-SNE of VST-transformed data", x="t-SNE X", y="t-SNE Y") +
  theme_minimal() +
  theme(legend.position="right", legend.title=element_blank(),
        legend.text=element_text(size=12), title=element_text(size=14),
        axis.title=element_text(size=12))

print(p_tsne)
ggsave(paste0("Plots/tSNE_vstdata_elp.png"), plot = p_tsne, dpi = 300, width = 6, height = 6)

## DEGS ####
# Differential expression for Standard Diet with IGF1 vs. Standard Diet without IGF1
res_ctrl <- results(dds, contrast=c("condition", "Ctrl_IGF1", "Ctrl"))
rownames(res_ctrl) <- counts_data$X

# Differential expression for LPHC Diet with IGF1 vs. LPHC Diet without IGF1
res_lphc <- results(dds, contrast=c("condition", "LPHC_IGF1", "LPHC"))
rownames(res_lphc) <- counts_data$X

# Filter DEGs for Ctrl_IGF1 vs. Ctrl
deg_ctrl <- res_ctrl[which(res_ctrl$padj < 0.05 & abs(res_ctrl$log2FoldChange) > 1), ]
# Export deg_ctrl to CSV
write.csv(deg_ctrl, file="deg_ctrl.csv", row.names=TRUE)

# Filter DEGs for LPHC_IGF1 vs. LPHC
deg_lphc <- res_lphc[which(res_lphc$padj < 0.05 & abs(res_lphc$log2FoldChange) > 1), ]
# Export deg_lphc to CSV
write.csv(deg_lphc, file="deg_lphc.csv", row.names=TRUE)

# Checking the number of DEGs for both comparisons
n_deg_ctrl <- nrow(deg_ctrl)
n_deg_lphc <- nrow(deg_lphc)

n_deg_ctrl
n_deg_lphc


## volcano plots
source("volcano_function.R")
# Plotting the volcano plot for Ctrl vs Ctrl_IGF1
p1 <- plot_volcano(res_ctrl, "Ctrl vs Ctrl_IGF1")
print(p1)
ggsave(paste0("Plots/volcano_Ctrl_Ctrl_IGF1.png"), plot = p1, dpi = 300, width = 6, height = 6)

# Plotting the volcano plot for LPHC vs LPHC_IGF1
p2 <- plot_volcano(res_lphc, "LPHC vs LPHC_IGF1")
print(p2)
ggsave(paste0("Plots/volcano_LPHC_LPHC_IGF1.png"), plot = p2, dpi = 300, width = 6, height = 6)

## Pathway analysis ####
# Now, map the gene IDs to Entrez IDs again
deg_ctrl_entrez <- mapIds(org.Mm.eg.db, keys=row.names(deg_ctrl), keytype="ENSEMBL", column="ENTREZID")
deg_lphc_entrez <- mapIds(org.Mm.eg.db, keys=row.names(deg_lphc), keytype="ENSEMBL", column="ENTREZID")

# Remove any NA values (genes that couldn't be mapped)
deg_ctrl_entrez <- deg_ctrl_entrez[!is.na(deg_ctrl_entrez)]
deg_lphc_entrez <- deg_lphc_entrez[!is.na(deg_lphc_entrez)]

# GO and KEGG ####
# GO enrichment analysis
ego_ctrl <- enrichGO(gene = deg_ctrl_entrez, 
                     OrgDb = org.Mm.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

ego_lphc <- enrichGO(gene = deg_lphc_entrez, 
                     OrgDb = org.Mm.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

# View the top enriched GO terms
head(ego_ctrl)
head(ego_lphc)
write.csv(ego_ctrl, file="ego_ctrl.csv", row.names=TRUE)
write.csv(ego_lphc, file="ego_lphc.csv", row.names=TRUE)

# KEGG enrichment analysis
kegg_ctrl <- enrichKEGG(gene = deg_ctrl_entrez, 
                        organism = 'mmu', # Mouse
                        pvalueCutoff = 0.05)

kegg_lphc <- enrichKEGG(gene = deg_lphc_entrez, 
                        organism = 'mmu', # Mouse
                        pvalueCutoff = 0.05)

# View the top enriched KEGG pathways
head(kegg_ctrl)
head(kegg_lphc)
write.csv(kegg_ctrl, file="kegg_ctrl.csv", row.names=TRUE)
write.csv(kegg_lphc, file="kegg_lphc.csv", row.names=TRUE)


#Visualization
# Dot plot for ego_ctrl
plot <- dotplot(ego_ctrl, showCategory=10) + ggtitle("GO Enrichment - Control")
ggsave(paste0("Plots/GO_Control.png"), plot = plot, dpi = 300, width = 6, height = 5)

# Dot plot for ego_lphc
plot <- dotplot(ego_lphc, showCategory=10) + ggtitle("GO Enrichment - LPHC")
ggsave(paste0("Plots/GO_LPHC.png"), plot = plot, dpi = 300, width = 6, height = 5)

# Dot plot for kegg_ctrl
plot <- dotplot(kegg_ctrl, showCategory=10) + ggtitle("KEGG Enrichment - Control")
ggsave(paste0("Plots/KEGG_Control.png"), plot = plot, dpi = 300, width = 6, height = 5)

# Dot plot for kegg_lphc
plot <- dotplot(kegg_lphc, showCategory=10) + ggtitle("KEGG Enrichment - LPHC")
ggsave(paste0("Plots/KEGG_LPHC.png"), plot = plot, dpi = 300, width = 6, height = 5)

# Visualize the most significantly enriched pathway for control
pathview(gene.data = setNames(kegg_ctrl$GeneRatio, kegg_ctrl$Description)[1], pathway.id = kegg_ctrl$ID[1], species = "mmu")

# Visualize the most significantly enriched pathway for LPHC
pathview(gene.data = setNames(kegg_lphc$GeneRatio, kegg_lphc$Description)[1], pathway.id = kegg_lphc$ID[1], species = "mmu")

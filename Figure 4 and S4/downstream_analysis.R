source("functions.R")
library(tidyverse)
library(ggplot2)
library(org.Mm.eg.db)
library(clusterProfiler)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(stats)
library(stringr)
library(DESeq2)
library(reshape2)

set.seed(8564)

#Load data from csv
data.dm <- read.csv("WGCNA files/geneInfo_from_darkmagenta_merged.csv",header = T, row.names = 1)
data.dg <- read.csv("WGCNA files/geneInfo_from_darkgrey_merged.csv",header = T, row.names = 1)
data.sb <- read.csv("WGCNA files/geneInfo_from_steelblue_merged.csv",header = T, row.names = 1)

deg_ctrl_csv <- read.csv("deg_ctrl.csv")
deg_lphc_csv <- read.csv("deg_lphc.csv")
common_genes <- intersect(deg_ctrl_csv[,1], deg_lphc_csv[,1])

# Check if genes in module are DEGs
geneDEG.dm <- data.dm[data.dm$ensembl.id %in% common_genes,]
geneDEG.dg <- data.dg[data.dg$ensembl.id %in% common_genes,]
geneDEG.sb <- data.sb[data.sb$ensembl.id %in% c(deg_ctrl_csv$X, deg_lphc_csv$X),]

# universe
expr_matrix <- read.csv("bigomics/counts.csv", row.names = 1)
keep_genes <- rowSums(expr_matrix > 10) > (0.5 * ncol(expr_matrix))
expr_matrix_filtered <- expr_matrix[keep_genes, ]
universe <- mapIds(org.Mm.eg.db, keys = rownames(expr_matrix_filtered), keytype = "ENSEMBL", column = "ENTREZID")
universe <- universe[!is.na(universe)]

# Map ENSEMBL IDs to ENTREZ IDs
data <- list(geneDEG.dm,geneDEG.dg,geneDEG.sb) # run again with each modules
module <- c("darkmagenta","darkgrey","steelblue")

for (i in 1:length(module)) {
  tmp <- data[[i]]
  
  geneDEG_entrez <- mapIds(org.Mm.eg.db, keys = tmp$ensembl.id, keytype = "ENSEMBL", column = "ENTREZID")
  
  # Remove any NA values (genes that couldn't be mapped)
  geneDEG_entrez <- geneDEG_entrez[!is.na(geneDEG_entrez)]
  
  # GO enrichment analysis
  ego_geneDEG <- enrichGO(gene = geneDEG_entrez, 
                         universe = universe, 
                         OrgDb = org.Mm.eg.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)
  
  # View the top enriched GO terms
  #head(ego_geneDEG)
  
  # KEGG enrichment analysis
  kegg_geneDEG <- enrichKEGG(gene = geneDEG_entrez,
                            universe = universe, 
                            organism = 'mmu', 
                            pvalueCutoff = 0.05)

  #head(kegg_geneDEG)
  
  # Dot plot for ego_common
  plot <- dotplot(ego_geneDEG, showCategory=10) + ggtitle(paste0("GO Enrichment - ", module[i]))
  ggsave(paste0("Plots/WGCNA/GO_deg_", module[i], ".png"), plot = plot, dpi = 300, width = 6, height = 5)
  write.csv(ego_geneDEG, file=paste0("WGCNA files/GO_deg_", module[i], ".csv"), row.names=TRUE)
  
  # Dot plot for kegg_common
  plot <- dotplot(kegg_geneDEG, showCategory=10) + ggtitle(paste0("KEGG Enrichment - ", module[i]))
  ggsave(paste0("Plots/WGCNA/KEGG_deg_", module[i], ".png"), plot = plot, dpi = 300, width = 6, height = 5)
  write.csv(kegg_geneDEG, file=paste0("WGCNA files/KEGG_deg_", module[i], ".csv"), row.names=TRUE)
}


## Combined graph
data <-geneDEG.sb # do the same with geneDEG.dg, geneDEG.sb

geneDEG_entrez <- mapIds(org.Mm.eg.db, keys = data$ensembl.id, keytype = "ENSEMBL", column = "ENTREZID")

# Remove any NA values (genes that couldn't be mapped)
geneDEG_entrez <- geneDEG_entrez[!is.na(geneDEG_entrez)]

# KEGG enrichment analysis
kegg_geneDEG.sb <- enrichKEGG(gene = geneDEG_entrez, # do the same with dg and sb
                           universe = universe, 
                           organism = 'mmu', 
                           pvalueCutoff = 0.05)

# Function to remove the unwanted text and update the description
clean_description <- function(description) {
  gsub("\\s*-\\s*Mus musculus \\(house mouse\\)", "", description)
}
# Function to extract necessary information and assign module name
extract_info <- function(res, module_name) {
  data <- as.data.frame(res)
  data$module <- module_name
  data$Description <- sapply(data$Description, clean_description)
  data$Description <- str_wrap(data$Description, width = 30)  # Wrap long pathway descriptions
  data.frame(
    Pathway = data$Description,
    pvalue = -log10(data$p.adjust),
    module = data$module
  )
}

# Convert the enrichResult objects into data frames with module names
df_dm <- extract_info(kegg_geneDEG.dm, "darkmagenta")
df_dg <- extract_info(kegg_geneDEG.dg, "darkgrey")
df_sb <- extract_info(kegg_geneDEG.sb, "steelblue")

# Combine the data frames
plot_data <- rbind(df_dm, df_dg, df_sb)

# Sort the data based on pvalue to get the top pathways

# Apply the function to the Description column
plot_data$Pathway <- sapply(plot_data$Pathway, clean_description)

plot_data <- plot_data[order(plot_data$pvalue, decreasing = TRUE), ]

# Plot the data using ggplot2
custom_colors <- c(darkmagenta = "darkmagenta", darkgrey = "darkgrey", steelblue = "steelblue")

p <- ggplot(plot_data, aes(x=reorder(Pathway, -pvalue), y=pvalue, fill=module)) + 
  geom_bar(stat="identity", width=0.5) +  # Adjust the width of the bars
  coord_flip() + 
  labs(title="Enriched Pathways", 
       x="Pathways", 
       y="-log10(Adjusted p-value)") +
  scale_fill_manual(values = custom_colors, name = "Module") +
  theme_minimal() +
  theme(text = element_text(size = 14),  # Increase global text size
        axis.text.y = element_text(size = 12))  # Increase axis text size

print(p)
ggsave(paste0("Plots/WGCNA/KEGG_common_deg_sign_modules.png"), plot = p, dpi = 300, width = 6.5, height = 5)




#################################
## z-score for heatmap DEG gene in modules
expr_matrix <- read.csv("bigomics/counts.csv", row.names = 1)

# Creating the metadata
sample_names <- colnames(expr_matrix)  # exclude the gene ID column
conditions <- c(rep("Ctrl", 5), rep("Ctrl_IGF1", 3), rep("LPHC_IGF1", 5), rep("LPHC", 5))
coldata <- data.frame(row.names=sample_names, condition=conditions)

# Renaming the levels
coldata$condition <- factor(coldata$condition, levels = c("Ctrl", "Ctrl_IGF1", "LPHC_IGF1", "LPHC"))
levels(coldata$condition) <- c("Ctrl", "Ctrl_IGF1", "LPHC_IGF1", "LPHC")

# Filter out lowly expressed genes
keep_genes <- rowSums(expr_matrix > 10) > (0.5 * ncol(expr_matrix))
expr_matrix_filtered <- expr_matrix[keep_genes, ]

# Normalization using DESeq2 (as already done)
dds <- DESeqDataSetFromMatrix(countData = expr_matrix_filtered, colData = coldata, design = ~1)
dds <- estimateSizeFactors(dds)
vst_transform <- vst(dds, blind=FALSE)  # Directly get VST transformation
datExpr <- as.data.frame(assay(vst_transform))
datExpr$gene <- rownames(datExpr)

modGenes <- c(geneDEG.dm$ensembl.id,geneDEG.dg$ensembl.id,geneDEG.sb$ensembl.id)

exp_modGenes <- datExpr %>% dplyr::filter(gene %in% modGenes)

exp_modGenes$Symbol <- mapIds(org.Mm.eg.db, rownames(exp_modGenes), keytype = "ENSEMBL", column = "SYMBOL")
exp_modGenes <- na.omit(exp_modGenes)
rownames(exp_modGenes) <- exp_modGenes$Symbol
exp_modGenes <- subset(exp_modGenes, select = - c(gene,Symbol))

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
mean_modGenes <- apply(exp_modGenes, 1, function(x)by(x, all.group, mean))
mean_exp_modGenes_zscore <- t(apply(mean_modGenes, 1, cal_z_score))

# Perform percentile transformation
percentile_matrix <- percentile_transformation(exp_modGenes)
mean_percentile_matrix <- apply(percentile_matrix, 1, function(x)by(x, all.group, mean))
mean_exp_modGenes_zscore <- mean_percentile_matrix

# heatmap of gene expression from the selected module
module <- na.omit(data.frame(c(geneDEG.dm$moduleColor,geneDEG.dg$moduleColor,geneDEG.sb$moduleColor),
                             c(geneDEG.dm$Symbol,geneDEG.dg$Symbol,geneDEG.sb$Symbol)))

row.names(module) <- module[,2]
colnames(module)[1] <- "module"
module <- subset(module, select = - 2)

my_colour = list("module" = 
                   c("darkmagenta" = "darkmagenta", 
                     "darkgrey" = "darkgrey",
                     "steelblue" = "steelblue"))

pdf("Plots/WGCNA/Heatmap_deg_module_p_transform.pdf",onefile=T,width=13,height=6)
mybreaks <- c(
  #seq(min(mean_kid_modGenes), -0.01, length.out=50), #_zscore
  #seq(0, max(mean_kid_modGenes),length.out=50)
  seq(0, 1,length.out=100)
) 
color <- viridis(99) #colorRampPalette(c("blue","white","red"))(100)
pheatmap(mean_exp_modGenes_zscore,#_zscore,
         cluster_rows = F,
         cluster_cols = T,
         #show_colnames = TRUE,
         #show_rownames = TRUE,
         main = "Percentile-score of mean gene expression",
         color = color,
         #breaks = mybreaks,
         #cutree_col = 4, 
         annotation_col = module, 
         annotation_colors = my_colour, 
         annotation_legend = TRUE,
         border_color = NA,
         cellheight = 15,
         cellwidth = 12,
         fontsize_col = 12
         #filename = "heatmap.pdf"
)
dev.off()


##################################
# Individual gene plots

library(ggpubr)

expr_matrix <- read.csv("bigomics/counts.csv", row.names = 1)

# Creating the metadata
sample_names <- colnames(expr_matrix)  # exclude the gene ID column
conditions <- c(rep("Ctrl", 5), rep("Ctrl_IGF1", 3), rep("LPHC_IGF1", 5), rep("LPHC", 5))
coldata <- data.frame(row.names=sample_names, condition=conditions)

# Renaming the levels
coldata$condition <- factor(coldata$condition, levels = c("Ctrl", "Ctrl_IGF1", "LPHC_IGF1", "LPHC"))
levels(coldata$condition) <- c("Ctrl", "Ctrl_IGF1", "LPHC_IGF1", "LPHC")

# Filter out lowly expressed genes
keep_genes <- rowSums(expr_matrix > 10) > (0.5 * ncol(expr_matrix))
expr_matrix_filtered <- expr_matrix[keep_genes, ]

# Normalization using DESeq2 (as already done)
dds <- DESeqDataSetFromMatrix(countData = expr_matrix_filtered, colData = coldata, design = ~1)
dds <- estimateSizeFactors(dds)
vst_transform <- vst(dds, blind=FALSE)  # Directly get VST transformation
datExpr <- as.data.frame(assay(vst_transform))

geneList <- c("Top2a","Col3a1","Flvcr1","G6pc","Irs2","Hsd11b1")

# Convert ensembl ID to gene names
datExpr$Symbol <- mapIds(org.Mm.eg.db, rownames(datExpr), keytype = "ENSEMBL", column = "SYMBOL")
datExpr <- na.omit(datExpr)

datExpr$Symbol <- make.unique(as.character(datExpr$Symbol))
rownames(datExpr) <- datExpr$Symbol
datExpr <- subset(datExpr, select = - c(Symbol))

gene2graph <- datExpr[rownames(datExpr) %in% geneList,]
gene2graph$gene <- rownames(gene2graph)

datMelted <- melt(gene2graph, id.vars = c("gene"))

# Add the conditions/groups based on the 'variable' column
datMelted$condition <- coldata$condition[datMelted$variable]


##Stats
# Split the condition into diet and treatment
datMelted$diet <- ifelse(grepl("Ctrl", datMelted$condition), "Ctrl", "LPHC")
datMelted$treatment <- ifelse(grepl("IGF1", datMelted$condition), "IGF1", "Vehicle")

# Function to run two-way ANOVA and Tukey's post-hoc test for each gene
anova_test <- function(data){
  fit <- aov(value ~ diet * treatment, data = data)
  pval_interaction <- summary(fit)[[1]]["diet:treatment", "Pr(>F)"]
  
  # Tukey's post-hoc test
  posthoc <- TukeyHSD(fit, which = "diet:treatment")
  
  # Extract p-values from Tukey's HSD result
  if (!is.null(posthoc[['diet:treatment']])) {
    pval_tukey <- posthoc[['diet:treatment']][, "p adj"]
  } else {
    pval_tukey <- NULL
  }
  
  return(list(pval_interaction = pval_interaction, pval_tukey = pval_tukey))
}


# Extract p-values for interaction
pvals_interaction <- sapply(results, function(res) res$pval_interaction)

# Extract Tukey's post-hoc p-values (This will be a list of p-values for each pairwise comparison)
pvals_tukey <- lapply(results, function(res) res$pval_tukey)

# For the sake of simplicity, we'll just display the interaction term p-value on the plot. 
# You can explore the `pvals_tukey` list separately to examine the pairwise comparisons.

# Plot
# Compute the maximum value for each gene to position the p-values
max_values <- aggregate(value ~ gene, data=datMelted, max)

# Merge the max values and p-values
annotation_data <- merge(max_values, data.frame(gene = names(pvals_interaction), pvalue = unname(pvals_interaction)), by="gene")

p <- ggplot(datMelted, aes(x=condition, y=value)) + 
  geom_boxplot() +
  
  # Y-axis label
  ylab("Transcript expression \nlevel (VST normalized)") +
  
  # X-axis label (will also be used for the legend)
  xlab('') +
  
  # Use facet_wrap to display the gene name at the top, all in one line
  facet_wrap(~gene, scales="free_y", ncol=length(unique(datMelted$gene))) +
  
  # Add statistical comparisons (annotate with the p-value of the interaction term)
  geom_text(data=annotation_data, 
            aes(x=2.5, y=value + 0.5,  # Adjust y position here
                label=sprintf("p=%.3f", pvalue)), hjust=1.1, vjust=0.5, size=3) +  # Adjust font size if needed
  
  # Theme adjustments
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10))

print(p)
ggsave(paste0("Plots/Common_deg_sign_from_modules.png"), plot = p, dpi = 300, width = 12, height = 3)


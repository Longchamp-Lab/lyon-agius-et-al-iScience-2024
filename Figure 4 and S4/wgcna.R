# Load necessary libraries
source("functions.R")
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(WGCNA)
library(flashClust)
library(org.Mm.eg.db)
library(visNetwork)

set.seed(8564)


# Load the expression matrix
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

# Normalization using DESeq2 (as you've already done)
dds <- DESeqDataSetFromMatrix(countData = expr_matrix_filtered, colData = coldata, design = ~1)
dds <- estimateSizeFactors(dds)
vst_transform <- vst(dds, blind=FALSE)  # Directly get VST transformation
datExpr <- t(assay(vst_transform))

gsg=goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK # if TRUE ok to continue


# Sample Clustering
# Compute the dissimilarity based on correlation
dissimilarity <- 1 - cor(t(datExpr))
sampleTree <- flashClust(as.dist(dissimilarity), method="average")

pdf(file = "Plots/WGCNA/1-n-sampleClustering.pdf", width = 12, height = 9);
par(cex = 1.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

# Calculate sample distance and cluster the samples
sampleTree = hclust(dist(datExpr), method = "average");
# plot sample tree
pdf(file = "Plots/WGCNA/1b-n-sampleClustering.pdf", width = 12, height = 9);
par(cex = 1.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

#filter outliers samples
#keep_samples <- !rownames(datExpr) %in% c("Ctrl3", "Ctrl5", "LPHC.IGF1_2")
#datExpr_filtered <- datExpr[keep_samples, ]

# Compute the dissimilarity based on correlation
dissimilarity <- 1 - cor(t(datExpr_filtered))
sampleTree <- flashClust(as.dist(dissimilarity), method="average")

pdf(file = "Plots/WGCNA/1-n-sampleClustering_filtered.pdf", width = 12, height = 9);
par(cex = 1.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

#===============================================================================
#
#  Choose soft threshold parameter
#
#===============================================================================

# Choose a set of soft threshold parameters
powers = c(c(1:20), seq(from = 22, to=50, by=2))
sft = pickSoftThreshold(datExpr, 
                        powerVector = powers,
                        networkType = "signed hybrid",
                        verbose = 5)
## from this plot, we would choose a power of 11 because it's the lowest power for 
## which the scale free topology index reaches the plateau above 0.80

## allowWGCNAThreads(nThreads = 22)
# sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, 
#                         #blockSize = 16700, 
#                         networkType="signed",
#                         corFnc = "bicor", 
#                         corOptions = list(use = 'p', maxPOutliers = 0.1))

# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file = "Plots/WGCNA/2-n-sft.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


#===============================================================================
#
#  Turn data expression into topological overlap matrix
#
#===============================================================================

# Turn data expression into topological overlap matrix
power=sft$powerEstimate #power=5

#TOM = TOMsimilarityFromExpr(datExpr, power = power)
#dissTOM = 1-TOM 

#adjacency=adjacency(datExpr, power=softPower, type="signed hybrid")
adjacency = adjacency(datExpr, power = power, #nThreads = 22,
                      type = "signed hybrid");
#corFnc = "bicor", corOptions = "use = 'p', maxPOutliers = 0.1");

#save(sft, adjacency, file = "adjacency.signedhybrid.softpwer.RData")

## Topological Overlap Matrix (TOM)
# Turn adjacency into topological overlap, i.e. translate the adjacency into 
# topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType = "signed", verbose = 5);
dissTOM = 1-TOM;
colnames(TOM) = rownames(TOM) = colnames(datExpr) # from PKLab Harvard 
save(TOM, file = "./TOMsimilarity.signed.RData")

#===============================================================================
#
#  Construct modules
#
#===============================================================================

##load TOM map
load("WGCNA files/TOMsimilarity.signed.RData")

# Module identification using dynamic tree cut
## generate a clustered gene tree
geneTree <- flashClust(as.dist(dissTOM), method="average")

pdf(file = "Plots/WGCNA/3-gene_cluster_alt.pdf", width = 12, height = 9);
plot(geneTree, xlab = "", sub = "", main = "Gene Clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# This sets the minimum number of genes to cluster into a module
minModuleSize <- 30 
# generating modules and assigning them colors
# Module identification using dynamic tree cut: 
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, 
                             deepSplit = 2, pamRespectsDendro = FALSE, 
                             minClusterSize = minModuleSize)
table(dynamicMods)

## Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

## Plot the dendrogram and colors underneath
pdf(file = "Plots/WGCNA/4-module_tree_alt.pdf", width = 8, height = 6);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

#===============================================================================
#
#  Merge modules
#
#===============================================================================

## Merging of modules whose expression profiles are very similar
# Calculate eigengenes
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs) #or bicor
# Cluster module eigengenes
METree <- flashClust(as.dist(MEDiss), method= "average")


#plots tree showing how the eigengenes cluster together
MEDissThres <- 0.3
plot(METree, main= "Clustering of module eigengenes", xlab = "", sub = "") %>%
  abline(h = MEDissThres, col = "red")

## Note: the tutorial chooses a height cut of 0.25, corresonding to correlation
## of 0.75, to merge
# set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
#MEDissThres <- 0.25;

pdf(file = "Plots/WGCNA/5-Module_Dendro.pdf", wi = 8, he = 6)
par(mar = c(.75, 2.75, 0.6, 1) + 0.1)            # The default is ‘c(5, 4, 4, 2) + 0.1’ c(bottom, left, top, right)’.
plot(METree, main= "Clustering of module eigengenes", xlab = "", sub = "") %>%
  abline(h = MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, 
                           corFnc = "bicor", verbose = 5)
# The merged module colors
mergedColors <- merge$colors
length(unique(mergedColors))
table(mergedColors)
# Eigengenes of the new merged modules:
mergedMEs <- merge$newMEs
length(mergedMEs)

# formated
pdf(file = "Plots/WGCNA/5-merged_geneDendro_format.pdf", wi = 6, he = 4)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    cex.rowText = 1.3)
dev.off()

write.table(merge$oldMEs,file="WGCNA files/oldMEs.txt");
write.table(merge$newMEs,file="WGCNA files/newMEs.txt");

## In subsequent analysis, the merged module colors in mergedColors will be used
# Rename to moduleColors
moduleColors <- mergedColors
# Construct numerical labels corresonding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs
dim(MEs)

#=====================================================================================
#
#  Correlation between gene modules and traits
#
#=====================================================================================

## Correlate eigengenes with external traits and look for the most significant associations
# Define numbers of genes and samples
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
# Recalculate MEs with color labels
MEsO <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEsO)
# Read data as traits
traits = read.table("traits.csv", header = T, sep = ",")
rownames(traits) = traits[, 1]
traits = traits[, -c(1)]
traits = traits[, -c(1)]

# sample names should be consistent in eigen genes and traits !!!!
traits = traits[match(rownames(MEs), rownames(traits)), ]
table(rownames(MEs) == rownames(traits))
moduleTraitCor <- bicor(MEs, traits, use = "p") #bicor, cor
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

write.table(moduleTraitCor,file="WGCNA files/moduleTrait_correlation.txt")
write.table(moduleTraitPvalue,file="WGCNA files/moduleTrait_pValue.txt")

## Color code each association by correlation value
# sizeGrWindow(13, 10)
# Will display correlations and their p-values
pdf(file = "Plots/WGCNA/mod-trait.relations.pdf", width = 12, height = 9)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(5, 10.1, 3, 1.5) + 0.1);            # The default is ‘c(5, 4, 4, 2) + 0.1’.
# Display the correlation values within a heatmap plot
labeledHeatmap.multiPage(Matrix = moduleTraitCor,
                         xLabels = names(traits),
                         yLabels = names(MEs),
                         ySymbols = names(MEs),
                         colorLabels = FALSE,
                         colors = blueWhiteRed(50),
                         textMatrix = textMatrix,
                         setStdMargins = FALSE,
                         cex.text = 0.7,
                         textAdj = c(0.5, 0.5),
                         zlim = c(-1,1),
                         maxColsPerPage = 15,
                         maxRowsPerPage = 25,
                         main = paste("Module-trait relationships"))
dev.off()

#=====================================================================================
#
#  Plot heatmap of module-traits relationship
#
#=====================================================================================

# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(datExpr, traits$hIGF1 , use = 'p')
gene.signf.corr.pvals <- data.frame(corPvalueStudent(gene.signf.corr, nSamples))

gene.signf.corr.pvals$symbol <- mapIds(org.Mm.eg.db, rownames(gene.signf.corr.pvals),
                                       keytype = "ENSEMBL", column = "SYMBOL")

write.csv(gene.signf.corr.pvals %>% 
            as.data.frame() %>% 
            arrange(corPvalueStudent.gene.signf.corr..nSamples.), "WGCNA files/gene.signf.corr.pval.csv")

#=====================================================================================
#
#   Intramodular analysis: identifying genes with high geneModuleMembership & geneTraitSignificance
#
#=====================================================================================

# names (colors) of the modules
modNames = substring(names(MEs), 3)

MET = orderMEs(cbind(MEs, traits))

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");


geneTraitSignificance = as.data.frame(cor(datExpr, traits, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(traits), sep="");
names(GSPvalue) = paste("p.GS.", names(traits), sep="");

# Plot the dendrogram
pdf(file = "Plots/WGCNA/6-Eigengene_dendrogram.pdf", width = 6, height = 4.5)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
pdf(file = "Plots/WGCNA/6-Eigengene adjacency heatmap.pdf", width = 6, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MET,  "Eigengene adjacency heatmap", marHeatmap = c(10,10,1,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

#======================================

modules = "steelblue" #up: darkmagenta down: #darkgrey
# Rename to moduleColors
moduleColors = mergedColors
column = match(modules, modNames);
moduleGenes = moduleColors==modules;


pdf(file = "Plots/WGCNA/7-Module membership vs gene significance violet.pdf", width = 7, height = 7)
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 4]),
                   lmFnc = lm,
                   #abline = T,
                   abline.color = "red",
                   xlab = paste("Module Membership in", modules, "module"),
                   ylab = "Gene significance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "blue")
dev.off()

#=====================================================================================
#
#  Extract gene list from module
#
#=====================================================================================

names(merge$newMEs)
modules = c(substring(names(merge$newMEs)[5], 3)); #13 3 7
genes = colnames(datExpr)
inModule = is.finite(match(mergedColors,modules))
modGenes = data.frame(ensembl.id = genes[inModule],
                      Symbol = mapIds(org.Mm.eg.db, genes[inModule],
                                      keytype = "ENSEMBL", column = "SYMBOL"),
                      moduleColor = mergedColors[inModule],
                      geneTraitSignificance[genes[inModule],],
                      GSPvalue[genes[inModule],])

# Order modules by their significance for igf1
modOrder = order(-abs(cor(MEs, traits$hIGF1, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(modGenes)
  modGenes = data.frame(modGenes, geneModuleMembership[genes[inModule], modOrder[mod]], 
                        MMPvalue[genes[inModule], modOrder[mod]]);
  names(modGenes) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                      paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the modGenes variable first by module color, then by geneTraitSignificance
geneOrder = order(modGenes$moduleColor, -abs(modGenes$GS.hIGF1));
geneInfo = modGenes[geneOrder, ]

write.csv(geneInfo, paste("WGCNA files/geneInfo_from_",paste(modules,collapse = "_"), "_merged.csv", sep = ""))

#=====================================================================================
#
#  Visualizing Gene Network
#
#=====================================================================================

# Extract genes from the "darkmagenta" module
modules = "darkmagenta"
genes = colnames(datExpr)
inModule = is.finite(match(mergedColors, modules))

modGenes = data.frame(
  ensembl.id = genes[inModule],
  Symbol = mapIds(org.Mm.eg.db, genes[inModule],
                  keytype = "ENSEMBL", column = "SYMBOL"),
  moduleColor = mergedColors[inModule]
)

# Subset TOM matrix for genes in the module
TOM_module = TOM[inModule, inModule]
rownames(TOM_module) <- genes[inModule]
colnames(TOM_module) <- genes[inModule]

# Convert the TOM into an Edge List:
threshold = 0.02
adjacency_thresholded = ifelse(TOM_module > threshold, TOM_module, 0)

edges <- which(adjacency_thresholded != 0, arr.ind = TRUE)
edge_list <- data.frame(
  from = rownames(adjacency_thresholded)[edges[, 1]],
  to = rownames(adjacency_thresholded)[edges[, 2]],
  value = adjacency_thresholded[edges]
)

# Create node data frame using the Symbols as IDs
nodes <- na.omit(data.frame(id = modGenes$ensembl.id, 
                    group = as.factor(modGenes$moduleColor)))

library(igraph)
g <- graph_from_data_frame(edge_list, directed=FALSE, vertices=nodes)
plot(g)


edge_list_no_self_loops <- edge_list[edge_list$from != edge_list$to, ]
visNetwork(nodes, edge_list, width = "100%") %>%
  visLayout(randomSeed = 42, improvedLayout = TRUE)

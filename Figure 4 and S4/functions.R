if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("DESeq2")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

if (!requireNamespace("factoextra", quietly = TRUE)) {
  install.packages("factoextra")
}

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("clusterProfiler")
}

if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("org.Mm.eg.db")
}

if (!requireNamespace("pathview", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("pathview")
}

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("WGCNA")
}

if (!requireNamespace("flashClust", quietly = TRUE)) {
  install.packages("flashClust")
}

if (!requireNamespace("visNetwork", quietly = TRUE)) {
  install.packages("visNetwork")
}

if (!requireNamespace("Rtsne", quietly = TRUE)) {
  install.packages("Rtsne")
}

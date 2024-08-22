library(ggplot2)
library(org.Mm.eg.db)
library(ggrepel)
library(gridExtra)

plot_volcano <- function(res, title="", lfc =1, p.adj=0.05) {
  # Create a data frame from the DESeq2 results
  df <- as.data.frame(res)
  
  # Add gene symbols to the data frame
  df$Symbol <- mapIds(org.Mm.eg.db, rownames(df), 
                      keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
  
  # Define thresholds
  lfc_threshold <- lfc
  pval_threshold <- p.adj
  
  # Define color based on significance and fold change direction
  df$color <- ifelse(df$padj < pval_threshold & df$log2FoldChange > lfc_threshold, "Upregulated",
                     ifelse(df$padj < pval_threshold & df$log2FoldChange < -lfc_threshold, "Downregulated", "Neutral"))
  
  # Plot
  p <- ggplot(df, aes(x=log2FoldChange, y=-log10(padj), label=Symbol)) +
    
    # Add colored points based on direction of fold change
    geom_point(aes(color=color), na.rm=TRUE, alpha=0.75, size=2) +
    
    # Highlight top genes (here top by adjusted p-value)
    geom_text_repel(data=subset(df, padj < pval_threshold & abs(log2FoldChange) > lfc_threshold & rank(padj) <= 10),
                    nudge_y=1, box.padding=0.35, point.padding=0.5, segment.color='transparent',
                    size=5) +  # Adjust font size here
    
    # Add vertical lines at log2FoldChange thresholds and horizontal lines at p-value threshold
    geom_vline(xintercept = 0, linetype="solid", color = "black") +  # Centered vertical line
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype="dotted") +
    geom_hline(yintercept = -log10(pval_threshold), linetype="dotted") +
    
    # Set custom colors
    scale_color_manual(values=c("Upregulated"="red", "Downregulated"="blue", "Neutral"="lightgray")) +
    
    # Set titles and labels
    labs(title=title, x="log2(Fold Change)", y="-log10(p-adj)") +
    
    # Theme adjustments
    theme_bw() +
    theme(legend.position="none",
          panel.border=element_blank(),
          panel.grid.major=element_line(color="gray", size=0.5),
          axis.text.x=element_text(angle=45, hjust=1, size=16),  # Adjusted size for axis text (numbers)
          axis.text.y=element_text(size=16),  # Adjusted size for y axis text (numbers)
          axis.title.x=element_text(size=20),  # Adjusted size for x axis label
          axis.title.y=element_text(size=20),
          plot.title = element_text(size=20))  # Adjusted size for y axis label
  
  
  return(p)
}
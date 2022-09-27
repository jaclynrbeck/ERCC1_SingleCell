# Helper functions used by "Step02_QualityControl". Some of these functions
# may be unused after final code cleanup. 

# Author: Jaclyn Beck
# Copied from https://github.com/jaclynrbeck/TCellsAD2022, modified to use
# human genes

library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)


# Create a Seurat object from 10X data. Adds percent mitochondrial genes and
# percent ribosomal genes as metadata. We immediately filter out any cells
# that have < 200 distinct genes counted, and each gene must have at least 10 
# cells that express it
makeSeurat <- function(file_dir, gene_column = 2, min_cells = 10, 
                       min_features = 200) {
  data <- Read10X(file_dir, gene.column=gene_column)
  scRNA <- CreateSeuratObject(data, 
                              project="ERCC1", 
                              min.cells = min_cells, 
                              min.features = min_features, 
                              names.field = 2, 
                              names.delim = "_")
  
  rm(data)
  
  # We are counting mitochondrial ribosomal proteins here. These patterns will 
  # recognize genes starting with "MT-" or "MRPS" or "MRPL"
  scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
  scRNA[["percent.mt"]] <- scRNA[["percent.mt"]] + PercentageFeatureSet(scRNA, pattern = "^MRP[S|L]")
  
  # This pattern matches genes starting with "RPL" or "RPS"
  scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[S|L]")
  
  scRNA[["complexity"]] = scRNA[["nCount_RNA"]] / scRNA[["nFeature_RNA"]]
  
  scRNA
}


##### Unused functions #####

# Plots a cumulative line graph of UMI vs number of cells
# data: matrix read by Read10X
# combine_types = TRUE: all genotypes are merged into one line
# combine_types = FALSE: each genotype gets its own line
plotUMICurve <- function(data, combine_types = FALSE) {
  barcodes <- data.frame(Barcode = colnames(data),
                         Sample = str_replace(colnames(data), ".*_", ""))
  
  umis <- data.frame(Cells=c(), UMIs=c(), Sample=c())
  
  if (combine_types == TRUE) {
    cs <- colSums(data > 0)
    umis <- data.frame(Cells=1:length(cs),
                       UMIs=sort(cs, decreasing=TRUE), 
                       Sample="Combined")
  }
  else {
    for (s in levels(factor(barcodes$Sample))) {
      samp_bcs <- subset(barcodes, Sample == s)
      cs <- colSums(data[, samp_bcs$Barcode] > 0)
      scs <- data.frame(Cells=1:length(cs), 
                        UMIs=sort(cs, decreasing=TRUE), 
                        Sample=s)
      
      umis <- rbind(umis, scs)
    }
  }
  
  plt = ggplot(umis, aes(x=Cells, y=UMIs, color=Sample)) + 
    geom_line() +
    scale_x_log10(n.breaks=12) + scale_y_log10(n.breaks=12) +
    ylab("UMIs") + theme_bw() +
    scale_color_brewer(palette="Set1")
  plt
}


# Finds the confidence interval for nCount_RNA using its log-normal distribution.
# 2 sds gives a ~95% confidence interval, 3 gives ~99.7
getConfidenceInterval <- function(scRNA, sds = 2) {
  log.rna <- log10(scRNA$nCount_RNA)
  ci <- sd(log.rna)*sds
  high.cutoff <- 10^(mean(log.rna) + ci)
  low.cutoff <- 10^(mean(log.rna) - ci)
  
  list("high" = high.cutoff, "low" = low.cutoff)
}


# Plots an NxN grid of PC scatterplots. Quick way to plot multiple PC dimensions
# against each other.
plotPcaGrid <- function(scRNA, dims = 1:4, group = "orig.ident") {
  plots <- list()
  
  for (d1 in dims) {
    for (d2 in dims) {
      plt <- DimPlot(scRNA, dims=c(d1,d2), reduction = "pca", group.by = group, 
                     shuffle = TRUE, raster = TRUE) + 
        theme(legend.position = "none") + 
        labs(title = NULL)
      plots[[(d1-dims[1])*length(dims) + (d2-dims[1]+1)]] = plt
    }
  }
  
  plot_grid(plotlist = plots, nrow = length(dims), labels = c("A", "B", "C"))
}
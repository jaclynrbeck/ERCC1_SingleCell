# Plotting functions commonly used for multiple figures
# Author: Jaclyn Beck
# Copied from https://github.com/jaclynrbeck/TCellsAD2022, some functions
# are unused. Slightly modified for this data set.

library(ggplot2)
library(dplyr)
library(viridis)
library(reshape2)
library(stringr)
library(ggprism)

# Creates a dot plot of the top 5 most upregulated markers for each cluster
# scRNA: Seurat object
# markers_file: full file path to an RDS file that contains the output of FindAllMarkers
# excluded genes: a list of sets of genes to exclude from display
# FDR: p value cutoff for significance
dotplotTop5Markers <- function( scRNA, markers_file, excluded_genes, FDR = 0.01 ) {
  all.markers <- readRDS(markers_file)
  sig.markers <- subset(all.markers, p_val_adj <= FDR & pct.1 >= 0.1)
  sig.markers <- subset(sig.markers, !(gene %in% unlist(excluded_genes)))
  
  markers.to.plot <- sig.markers %>% group_by(cluster) %>%
    top_n(5, wt = avg_log2FC)
  
  plt <- doViridisDotPlot(scRNA, unique(markers.to.plot$gene))
}


# Creates a dotplot with specific viridis colors and a few other tweaks
# scRNA: Seurat object
# feats: vector of gene names to include in the plot
doViridisDotPlot <- function( scRNA, feats ) {
  plt <- DotPlot(scRNA, features = feats, dot.scale = 3) + 
    RotatedAxis() + coord_flip() + 
    scale_color_viridis(option = "A", direction=-1, begin=0.6) +
    theme(axis.text=element_text(size=8), axis.title = element_blank(),
          plot.background = element_rect(fill = "white"))
}


# Creates a bar graph of population distributions, with the 4 genotypes
# side by side for each cluster.
# scRNA: Seurat object
# geno.colors: vector of colors to use for genotype bars. Should be in the same
#              order as levels(scRNA$genotype)
# sig.df: a data frame containing details about which comparisons are 
#         significant and where to draw significance bars and stars. Can be
#         null if no significance. See "Figure5.R" for an example format, but
#         briefly it should have the following fields: 
#           cluster: cluster name, repeated to fill out the column
#           group1: the first genotype in each comparison
#           group2: the second genotype in each comparison, paired with group1
#           p.adj: the p value for each comparison
#           xmin: x coordinate of the left side of each significance bar
#           xmax: x coordinate of the right side of each significance bar
#           y.position: y coordinate of each significance bar
populationBarGraph <- function( scRNA, geno.colors, sig.df = NULL) {
  pop <- table(Idents(scRNA), scRNA$orig.ident)
  cell.counts <- colSums(pop)
  cell.counts <- do.call(rbind, lapply(1:nrow(pop), function(x) cell.counts))
  pct <- melt(pop / cell.counts)
  
  pop <- melt(pop)
  pop <- cbind(pop, pct$value*100)
  colnames(pop) <- c("Cluster", "Sample", "Count", "Percent")
  pop$Genotype <- str_replace(pop$Sample, "-[123]", "")
  
  # Just for ordering of clusters in the graph
  means = aggregate(pop, list(Clust = pop$Cluster), mean)
  pop$Cluster <- factor(pop$Cluster, 
                        levels = means$Clust[order(means$Percent, decreasing = TRUE)])
  pop$Genotype <- factor(pop$Genotype, levels = c("WT", "KO"))
  
  # Used to make error bars
  data <- pop %>% group_by(Cluster, Genotype) %>% 
    summarize(Pct = mean(Percent), SD = sd(Percent), Max = max(Percent), Min = min(Percent))
  data$Max[is.na(data$SD)] <- NA
  data$Min[is.na(data$SD)] <- NA
  
  # Part 1 of plot: bar graph with separate bars for each genotype, grouped
  # by cluster
  plt <- ggplot(data, aes(x = Cluster, y = Pct)) +
    geom_bar(aes(fill = Genotype), stat = "identity", color = "black", 
             position = "dodge", width = 0.6,
             size = 0.2) + 
    theme_classic() +
    xlab("Cluster") + ylab("Percent of Cell Population") + 
    scale_fill_discrete(type = geno.colors) +
    scale_x_discrete(guide = guide_axis(angle = 45)) + 
    theme(axis.text.x = element_text(size=8), axis.title.x = element_blank()) + 
    geom_errorbar(aes(group = Genotype, ymin = Min, ymax = Max), 
                  position = position_dodge(width = 0.6), width = 0.25, size = 0.2)
  
  # Part 2: Add significance bars
  if (!is.null(sig.df)) {
    plt <- plt + add_pvalue(sig.df, label = "p.adj.signif", label.size = 2.5,
                            xmin = "xmin", xmax = "xmax", 
                            tip.length = 0.01, bracket.size = 0.3, 
                            bracket.shorten = 0.04)
  }
  
  plt 
}

populationBarGraphZoom <- function( scRNA, remove.idents = c(), geno.colors, sig.df = NULL) {
  pop <- table(Idents(scRNA), scRNA$orig.ident)
  cell.counts <- colSums(pop)
  cell.counts <- do.call(rbind, lapply(1:nrow(pop), function(x) cell.counts))
  pct <- melt(pop / cell.counts)
  
  pop <- melt(pop)
  pop <- cbind(pop, pct$value*100)
  colnames(pop) <- c("Cluster", "Sample", "Count", "Percent")
  pop$Genotype <- str_replace(pop$Sample, "-[123]", "")
  
  # Just for ordering of clusters in the graph
  means = aggregate(pop, list(Clust = pop$Cluster), mean)
  pop$Cluster <- factor(pop$Cluster, 
                        levels = means$Clust[order(means$Percent, decreasing = TRUE)])
  pop$Genotype <- factor(pop$Genotype, levels = c("WT", "KO"))
  
  # Used to make error bars
  data <- pop %>% group_by(Cluster, Genotype) %>% 
    summarize(Pct = mean(Percent), SD = sd(Percent), Max = max(Percent), Min = min(Percent))
  data$Max[is.na(data$SD)] <- NA
  data$Min[is.na(data$SD)] <- NA
  
  data <- subset(data, !(Cluster %in% remove.idents))
  
  # Part 1 of plot: bar graph with separate bars for each genotype, grouped
  # by cluster
  plt <- ggplot(data, aes(x = Cluster, y = Pct)) +
    geom_bar(aes(fill = Genotype), stat = "identity", color = "black", 
             position = "dodge", width = 0.6,
             size = 0.2) + 
    theme_classic() +
    xlab("Cluster") + ylab("Percent of Cell Population") + 
    scale_fill_discrete(type = geno.colors) +
    scale_x_discrete(guide = guide_axis(angle = 45)) + 
    theme(axis.text.x = element_text(size=8), axis.title.x = element_blank()) + 
    geom_errorbar(aes(group = Genotype, ymin = Min, ymax = Max), 
                  position = position_dodge(width = 0.6), width = 0.25, size = 0.2)
}


# Downsamples each genotype to have the same number of cells. Returns a new
# Seurat object with the downsampled data.
downsample <- function( scRNA ) {
  counts <- table(scRNA$genotype)
  downsamp <- min(counts)
  cells <- c()
  
  set.seed(10)
  for (N in names(counts)) {
    subs <- filter(scRNA@meta.data, genotype == N)
    cells <- c(cells, sample(rownames(subs), downsamp))
  }
  
  scRNA.down <- subset(scRNA, cells = cells)
  scRNA.down
}


# Draws UMAPs split by genotype, colored by cluster, where the cells have been
# downsampled to be equal across genotypes
umapSplitDownsampled <- function( scRNA, ident.name = "clusters", clust.colors ) {
  scRNA.down <- downsample(scRNA)
  Idents(scRNA.down) <- scRNA.down@meta.data[,ident.name]
  
  plt <- DimPlot(scRNA.down, reduction = "umap", split.by = "genotype", 
                 ncol = 2, pt.size = 0.01) + 
    NoLegend() + scale_color_discrete(type = clust.colors)
  
  rm(scRNA.down)
  
  plt
}


# Draws UMAPs of gene expression with 5 columns per row, using the viridis
# color scheme and a few visual tweaks
# scRNA: Seurat object
# genes: vector of genes to make UMAPs for
# coord.fixed: whether to plot so the x and y axes of the UMAP have the same
#              scale, or not.
umapGeneExprPlots <- function( scRNA, genes, coord.fixed = TRUE ) {
  plt <- FeaturePlot(scRNA, features = genes,
                     order = TRUE, ncol = 5, pt.size = 0.01, 
                     coord.fixed = coord.fixed) &
    scale_color_viridis(begin = 0.1, option = "D") & 
    theme(axis.title = element_blank(), axis.text = element_blank(),
          title = element_text(size = 10), legend.position = "none",
          axis.line = element_blank(), axis.ticks = element_blank())
}


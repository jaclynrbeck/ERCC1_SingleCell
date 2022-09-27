# Reads 10X data from the output of cellranger aggr and turns it into a Seurat
# object that has been normalized. 

# Author: Jaclyn Beck

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
source(file.path("functions", "QualityControl_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))
source("Filenames.R")


##### Create the Seurat object #####

scRNA <- makeSeurat(dir_filtered_counts)

# Samples WT-3 and KO-3 were run in a different batch from the other samples
scRNA$batch <- "Batch 1"
scRNA$batch[scRNA$orig.ident %in% c("WT-3", "KO-3")] <- "Batch 2"
gc()


##### Exploratory plots ##### 

VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                            "percent.rb", "complexity"), ncol=3, pt.size = 0)

plt <- list()
for (L in levels(scRNA$orig.ident)) {
  
  plt[[L]] <- FeatureScatter(scRNA, feature1 = "nFeature_RNA", 
                             feature2 = "nCount_RNA", 
                             cells = colnames(scRNA)[scRNA$orig.ident == L],
                             shuffle=TRUE, pt.size = 0.1) +
    scale_y_log10(n.breaks=16) + scale_x_log10(n.breaks=16) + theme_bw()
}
(plt[[1]] | plt[[2]] | plt[[3]]) / (plt[[4]] | plt[[5]] | plt[[6]])

for (L in levels(scRNA$orig.ident)) {
  plt[[L]] <- FeatureScatter(scRNA, feature1 = "complexity", 
                             feature2 = "nCount_RNA", 
                             cells = colnames(scRNA)[scRNA$orig.ident == L],
                             shuffle=TRUE, pt.size = 0.1) +
    scale_y_log10(n.breaks=16) + scale_x_continuous(n.breaks=16) + theme_bw()
}
(plt[[1]] | plt[[2]] | plt[[3]]) / (plt[[4]] | plt[[5]] | plt[[6]])

for (L in levels(scRNA$orig.ident)) {
  plt[[L]] <- FeatureScatter(scRNA, feature1 = "nFeature_RNA", 
                             feature2 = "percent.mt",
                             cells = colnames(scRNA)[scRNA$orig.ident == L],
                             shuffle=TRUE, pt.size = 0.1) +
    scale_y_log10(n.breaks=16) + scale_x_log10(n.breaks=16) + theme_bw()
}
(plt[[1]] | plt[[2]] | plt[[3]]) / (plt[[4]] | plt[[5]] | plt[[6]])

for (L in levels(scRNA$orig.ident)) {
  plt[[L]] <- FeatureScatter(scRNA, feature1 = "nFeature_RNA", 
                             feature2 = "percent.rb", 
                             cells = colnames(scRNA)[scRNA$orig.ident == L],
                             shuffle=TRUE, pt.size = 0.1) +
    scale_y_continuous(n.breaks=16) + scale_x_log10(n.breaks=16) + theme_bw() +
    geom_hline(yintercept = 1.2)
}
(plt[[1]] | plt[[2]] | plt[[3]]) / (plt[[4]] | plt[[5]] | plt[[6]])

for (L in levels(scRNA$orig.ident)) {
  plt[[L]] <- FeatureScatter(scRNA, feature1 = "percent.rb", 
                             cells = colnames(scRNA)[scRNA$orig.ident == L],
                             feature2 = "percent.mt", 
                             shuffle=TRUE, pt.size = 0.1) +
    scale_y_log10(n.breaks=16) + scale_x_continuous(n.breaks=16) + theme_bw() +
    geom_vline(xintercept = 1.2)
}
(plt[[1]] | plt[[2]] | plt[[3]]) / (plt[[4]] | plt[[5]] | plt[[6]])


##### Remove low-quality cells #####

scRNA <- subset(scRNA, subset = nFeature_RNA >= 1200 & 
                  percent.mt <= 10 &
                  percent.rb >= 1.2 &
                  percent.rb <= 20 & 
                  complexity <= 4.0)

# Batch 1 had higher counts in general and needs a higher nFeature_RNA threshold
good <- colnames(scRNA)[scRNA$batch == "Batch 2" | 
                          scRNA$nFeature_RNA >= 2000]
scRNA <- subset(scRNA, cells = good)


##### Finalize un-normalized Seurat object #####

table(scRNA$orig.ident)

# Remove genes expressed in < 10 cells after cell filtering
scRNA <- removeLowExpressedGenes(scRNA, 10)
scRNA <- DietSeurat(scRNA, scale.data = FALSE)
scRNA$genotype <- str_replace(scRNA$orig.ident, "-[123]", "")

# At this point scRNA contains all cells that passed QC
saveRDS(scRNA, file = file_seurat_unnorm)
gc()

# Done

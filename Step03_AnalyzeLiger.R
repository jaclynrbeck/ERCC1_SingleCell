# There are significant batch differences between the two runs. After testing
# different methods of batch correction for clustering (Seurat integration, 
# Harmony, LIGER), I found LIGER to give the best correction and clustering.
#
# This script runs the LIGER pipeline on the unnormalized Seurat object to 
# integrate the samples, then creates a new Seurat object with the proper 
# metadata and the LIGER cell encodings in a dimension reduction object called 
# "inmf". 
#
# This script also runs liger's differential gene expression function
# on the clusters found by liger (which are the same as those found by Seurat
# later), using the Wilcoxon test. This was informational for finding bad
# clusters and removing them. 
#
# Further analysis of "good" clusters is done using the Seurat object and
# several figures are produced, corresponding to figures in my thesis. 
#
# Author: Jaclyn Beck

library(Seurat)
library(rliger)
library(stringr)
library(dplyr)
library(writexl)
library(ggplot2)
library(Matrix.utils)
library(pheatmap)
library(multcomp)
library(msigdbr)
library(fgsea)
library(emmeans)

source("Filenames.R")
source(file.path("functions", "Analysis_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "DPA_HelperFunctions.R"))
source(file.path("functions", "Figure_PlottingFunctions.R"))

scRNA <- readRDS(file_seurat_unnorm)

# Set up LIGER object
lig <- seuratToLiger(scRNA, combined.seurat = TRUE, meta.var = "orig.ident")
rm(scRNA)
gc()

lig <- normalize(lig)
lig <- selectGenes(lig)
lig <- scaleNotCenter(lig)

# k = 40 seems to give the best biologically-relevant shape / clusters
lig <- optimizeALS(lig, k = 40)
lig <- quantile_norm(lig)

# Using same UMAP parameters that Seurat uses, for consistency
lig <- runUMAP(lig, distance = "cosine", n_neighbors = 30, min_dist = 0.3)
lig <- louvainCluster(lig, resolution = 0.1)

plotByDatasetAndCluster(lig, axis.labels = c("UMAP1","UMAP2"))
calcAlignment(lig)

# Save intermediate step
saveRDS(lig, file_liger_norm)
gc()

#gene_loadings <- plotGeneLoadings(lig, do.spec.plot = FALSE, return.plots = TRUE)
#gene_loadings[[1]]

res <- runWilcoxon(lig, compare.method = "clusters")
res <- subset(res, padj <= 0.01)
res$Ensembl.Id <- geneNameToEnsembl(res$feature)

sheets <- list()

for (G in unique(res$group)) {
  res.sub <- subset(res, group == G)
  sheets[[paste0("Cluster ", G)]] <- res.sub
}

write_xlsx(sheets, path = file.path(dir_output, "diffgenes_liger_bycluster.xlsx"))


# Convert to Seurat object for the rest of the analysis
scRNA <- ligerToSeurat(lig, by.dataset = FALSE)
rm(lig)
gc()

# Add back missing metadata
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
scRNA[["percent.mt"]] <- scRNA[["percent.mt"]] + PercentageFeatureSet(scRNA, pattern = "^MRP[S|L]")

# This pattern matches genes starting with "RPL" or "RPS"
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[S|L]")

scRNA[["complexity"]] = scRNA[["nCount_RNA"]] / scRNA[["nFeature_RNA"]]

scRNA$genotype <- factor(str_replace(scRNA$orig.ident, "-[123]", ""), 
                         levels = c("WT", "KO"))

scRNA$batch <- "Batch 1"
scRNA$batch[scRNA$orig.ident %in% c("WT-3", "KO-3")] <- "Batch 2"

scRNA <- NormalizeData(scRNA)
scRNA <- RunUMAP(scRNA, reduction = "inmf", dims = 1:40)
Idents(scRNA) <- scRNA$orig.ident
DimPlot(scRNA, shuffle = TRUE)

scRNA <- FindNeighbors(scRNA, reduction = "inmf", dims = 1:40)
scRNA <- FindClusters(scRNA, resolution = 0.1)
DimPlot(scRNA, shuffle = TRUE, label = TRUE)
DimPlot(scRNA, shuffle = TRUE, label = FALSE, split.by = "orig.ident", ncol = 3)

# Remove bad clusters: 
#   7 has high expression of hemoglobin genes (HBA1, HBB, HBD)
#   8 appears to be macrophages (RNASE1, SELENOP, CD163, MRC1)
#   9 appears to be low quality cells? Lower than average percent.rb, higher 
#     average nCount_RNA, nFeature_RNA, and complexity
scRNA <- subset(scRNA, seurat_clusters %in% c(0:6))

Idents(scRNA) <- scRNA$seurat_clusters
scRNA <- RenameIdents(scRNA, list("0" = "Homeostatic",
                                  "1" = "Proliferating",
                                  "2" = "IFN 1 / IL-1\u03b2",
                                  "3" = "DAM",
                                  "4" = "IFN 2",
                                  "5" = "MHC II",
                                  "6" = "Stressed"))

scRNA$clusters <- Idents(scRNA)

DimPlot(scRNA, shuffle = TRUE, label = TRUE)

# Save final Seurat object
saveRDS(scRNA, file_seurat_analyzed)


# Differences between KO and WT cells
Idents(scRNA) <- scRNA$genotype

markers <- FindMarkers(scRNA, ident.1 = "KO", ident.2 = "WT", test.use = "MAST", 
                       latent.vars = c("nCount_RNA", "percent.rb", "percent.mt"),
                       logfc.threshold = 0.1, min.pct = 0.01)

sig.markers <- subset(markers, p_val_adj <= 0.01) %>% arrange(desc(avg_log2FC))
sig.markers$Ensembl.Id <- geneNameToEnsembl(rownames(sig.markers))

write.csv(sig.markers, file_markers_genotypes)
gc()

##### Heatmap of top 100 changed genes between genotypes #####
# Figure 2.9A of my thesis

markers <- sig.markers

excluded <- grep("^MT-|^RPL|^RPS", rownames(scRNA), value = TRUE)
markers <- markers[!(rownames(markers) %in% excluded),]

top100 <- top_n(markers, wt = abs(avg_log2FC), n = 100)
genes <- rownames(top100)

assay <- GetAssayData(scRNA, slot = "counts", assay = "RNA")
meta <- scRNA@meta.data

# Remove genes expressed in fewer than 10 cells
ok <- rowSums(assay >= 1) >= 10
assay <- assay[ok,]

bulk <- aggregate.Matrix(t(assay), groupings = meta$orig.ident, fun = "sum")

norm <- DESeq2::varianceStabilizingTransformation(as.matrix(t(bulk)), blind = TRUE)

plt <- pheatmap(norm[genes,], scale = "row", fontsize_row = 5, 
                treeheight_row = 0, treeheight_col = 0, cluster_cols = T)

ggsave(file.path(dir_figures, "fig2.9a.png"),
       plot = plt, width = 2.25, height = 8, units = "in", dpi = "print")


##### Batch differences #####

Idents(scRNA) <- scRNA$batch
markers <- FindMarkers(scRNA, ident.1 = "Batch 2", ident.2 = "Batch 1", 
                       test.use = "MAST", logfc.threshold = 0.1, min.pct = 0.01,
                       latent.vars = c("nCount_RNA", "percent.rb", "percent.mt"))

sig.markers <- subset(markers, p_val_adj <= 0.01) %>% arrange(desc(avg_log2FC))
sig.markers$Ensembl.Id <- geneNameToEnsembl(rownames(sig.markers))

write.csv(sig.markers, file_markers_batch)
gc()


##### Differences between clusters #####

Idents(scRNA) <- scRNA$clusters
markers <- FindAllMarkers(scRNA, test.use = "MAST", 
                          latent.vars = c("nCount_RNA", "percent.rb", "percent.mt"))
saveRDS(markers, file_markers_all)

sig.markers <- subset(markers, p_val_adj <= 0.01)

writeDifferentialGenes(sig.markers, file_markers_clusters)


excluded_genes <- grep("^MT-|^RPL|^RPS|XIST", rownames(scRNA), value = TRUE)

# Dotplot of top 5 markers per cluster -- color scheme and order of genes is 
# different than in the thesis figure, because I switched to using the function 
# I wrote for the T cell project, but data is the same. 
Idents(scRNA) <- scRNA$clusters

plt <- dotplotTop5Markers(scRNA, file_markers_all, excluded_genes, 
                          FDR = 0.01)
plt

ggsave(file.path(dir_figures, "fig2.8c.png"),
       plot = plt, width = 5.25, height = 5.25, units = "in", dpi = "print")


##### DPA #####

Idents(scRNA) <- scRNA$clusters

pop <- table(Idents(scRNA), scRNA$orig.ident)
cell.counts <- colSums(pop)
cell.counts <- do.call(rbind, lapply(1:nrow(pop), function(x) cell.counts))
pct <- melt(pop / cell.counts)

pop <- melt(pop)
pop <- cbind(pop, pct$value*100)
colnames(pop) <- c("Cluster", "Sample", "Count", "Percent")
pop$Genotype <- str_replace(pop$Sample, "-[123]", "")

res <- aov(Percent ~ Cluster * Genotype, data = pop)
summary(res) # This will show no sig. interaction between genotype and cluster

# .. but we'll run pair-wise comparisons for fun
mult <- glht(res, linfct=emm(pairwise~Genotype|Cluster, adjust="BH"))
summary(mult) # Nothing significant here


# Population bargraph - colors aren't exactly like in the thesis but close.

Idents(scRNA) <- scRNA$clusters
plt <- populationBarGraph(scRNA, c("#1A9865", "#E74557"), NULL)
plt

ggsave(file.path(dir_figures, "fig2.8b_all.png"),
       plot = plt, width = 5, height = 3.2, units = "in", dpi = "print")

plt <- populationBarGraphZoom(scRNA, remove.idents = c("Homeostatic"), 
                              geno.colors = c("#1A9865", "#E74557"), NULL)

plt

ggsave(file.path(dir_figures, "fig2.8b_zoom.png"),
       plot = plt, width = 4, height = 3.2, units = "in", dpi = "print")


##### GSEA #####

# This can be run on any dataset in MSigDB, but I found that Reactome terms were
# the most informative

Idents(scRNA) <- scRNA$genotype
DefaultAssay(scRNA) <- "RNA"

# Need to re-run this with wider parameters than before, for GSEA
markers <- FindMarkers(scRNA, ident.1 = "KO", ident.2 = "WT", 
                       logfc.threshold = 0.01, min.pct = 0, test.use = "MAST",
                       latent.vars = c("nCount_RNA", "percent.rb", "percent.mt"))
saveRDS(markers, file = file_gsea_markers)
gc()

# Only use genes expressed in at least 10 cells for GSEA. Gives 26300 genes
tmp <- GetAssayData(scRNA, slot = "counts")
tmp <- tmp > 0
ok <- rowSums(tmp)

good_genes <- names(ok[ok >= 10])
rm(tmp, ok)

as.data.frame(msigdbr_collections())
reac <- msigdbr(species = "human", category = "C2", subcategory = "CP:REACTOME")
go_bp <- msigdbr(species = "human", category = "C5", subcategory = "GO:BP")
imm <- msigdbr(species = "human", category = "C7", subcategory = "IMMUNESIGDB")
hall <- msigdbr(species = "human", category = "H")

gene_set <- reac #go_bp
gene_set_list <- split(x = gene_set$ensembl_gene, f = gene_set$gs_name)
set_name <- "Reactome" #"GO:BP"
replace_string <- "REACTOME" # "HALLMARK", "REACTOME", or "GOBP". ImmuneSigDB doesn't have anything to remove

markers <- readRDS(file_gsea_markers)
good_genes <- intersect(good_genes, rownames(markers))

# 9416 genes left with a LFC threshold of 0.01
ranks <- markers[good_genes, "avg_log2FC"]
names(ranks) <- geneNameToEnsembl(good_genes)
ranks <- ranks[order(ranks)]

set.seed(112233)

res <- fgsea(pathways = gene_set_list, stats = ranks, maxSize = 500)
res <- subset(res, !is.na(padj))
res <- subset(res, padj <= 0.05)

collapsed <- collapsePathways(res, gene_set_list, ranks)

res$leadingEdgeGeneSymbols <- lapply(res$leadingEdge, ensemblToGeneName)
res$leadingEdgeGeneSymbols <- sapply(res$leadingEdgeGeneSymbols, str_c, collapse = ", ")
res$leadingEdge <- sapply(res$leadingEdge, str_c, collapse = ", ")

gs_source <- subset(gene_set, gs_name %in% res$pathway)
gs_source <- unique(gs_source[,c("gs_name", "gs_exact_source")])

res <- merge(res, gs_source, by.x = "pathway", by.y = "gs_name")

res <- res[order(res$NES, decreasing = TRUE),]

res$IsMainPathway <- res$pathway %in% collapsed$mainPathways

res$pathway <- str_replace(res$pathway, replace_string, "") %>% 
  str_replace_all("_", " ") %>% str_trim() %>% str_to_title()

#res.display <- res[res$IsMainPathway]
res.display <- res
res.display$pathway <- str_trunc(res.display$pathway, width = 50, side = "right")


plt <- ggplot(res.display, aes(reorder(pathway, NES), NES, fill = NES)) +
  geom_col() +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=paste0(set_name, " pathway enrichment")) + 
  theme_classic() + theme(axis.text.y = element_text(size = 6))
plt

ggsave(file.path(dir_figures, "fig2.9b.png"), plot = plt,
       width = 5.5, height = 8, units = "in", dpi = "print")

write_xlsx(res, path = file.path(dir_go, "gsea_KOvWT_reactome.xlsx"))

# Done






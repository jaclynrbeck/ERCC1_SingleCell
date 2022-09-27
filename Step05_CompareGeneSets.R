# This script takes DGE sets from multiple papers / databases on aging and
# compares those genes to this data set. At the end this script generates 
# simple heatmaps to show overlap between the data sets. In these figures,
# the color-coding is as follows: 
#   Red = "significantly increased in aged cells / ERCC1 KO cells"
#   Light red = "increased but not significant"
#   Gray = "gene is missing from data set"
#   Light blue = "decreased but not significant"
#   Blue = "significantly decreased in aged cells / ERCC1 KO cells"

library(dplyr)
library(readxl)
library(writexl)
library(stringr)
library(Matrix.utils)
library(Seurat)

source(file.path("functions", "Analysis_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))
source("Filenames.R")

diff.genes <- read.csv(file_markers_genotypes)
colnames(diff.genes)[1] <- "Gene"

sheets <- list()
gene.sets <- list()

##### Lopes gene set (human aging) #####

lopes <- read_xlsx(file.path(dir_external, "Lopes_Aging_DEG.xlsx"), 
                   sheet = "Table14_age-related_genes", skip = 2)

lopes$ensembl_id <- str_replace(lopes$ensembl_id, "\\..*", "")

lopes$Direction <- sign(lopes$logFC)
lopes$Direction[lopes$adj.P.Val > 0.05] <- 0
lopes$Ensembl.Id <- lopes$ensembl_id
lopes$Gene.name <- lopes$symbol

gene.sets[["Human Mg Aging (Lopes 2022)"]] <- lopes

lopes <- subset(lopes, adj.P.Val <= 0.05)

lu <- subset(lopes, logFC > 0)
ld <- subset(lopes, logFC < 0)

genes <- unique(c(lu$symbol, ensemblToGeneName(lu$ensembl_id)))
matches <- intersect(genes, diff.genes$Gene)
sheets[["Lopes up"]] <- subset(diff.genes, Gene %in% matches)

genes <- unique(c(ld$symbol, ensemblToGeneName(ld$ensembl_id)))
matches <- intersect(genes, diff.genes$Gene)
sheets[["Lopes down"]] <- subset(diff.genes, Gene %in% matches)


##### Galatro gene set (human aging) #####

galatro <- list("Down" = read_xlsx(file.path(dir_external, "Galatro_Aging_DEG.xlsx"), sheet = "age_down"),
                "Up" = read_xlsx(file.path(dir_external, "Galatro_Aging_DEG.xlsx"), sheet = "age_up"))

names(galatro[["Down"]])[1] <- "ensembl_gene_id" # Correct a typo

galatro.df <- rbind(galatro[["Up"]], galatro[["Down"]])
galatro.df$Direction <- sign(galatro.df$logFC)
galatro.df$Direction[galatro.df$P.Value > 0.05] <- 0 
galatro.df$Ensembl.Id <- galatro.df$ensembl_gene_id
galatro.df$Gene.name <- galatro.df$external_gene_name

gene.sets[["Human Mg Aging (Galatro 2017)"]] <- galatro.df

gu <- galatro[["Up"]]

# Account for synonyms
genes <- unique(c(gu$external_gene_name, ensemblToGeneName(gu$ensembl_gene_id)))

matches <- intersect(genes, diff.genes$Gene)
sheets[["Galatro up"]] <- subset(diff.genes, Gene %in% matches)

gd <- galatro[["Down"]]
genes <- unique(c(gd$external_gene_name, ensemblToGeneName(gd$ensembl_gene_id)))

matches <- intersect(genes, diff.genes$Gene)
sheets[["Galatro down"]] <- subset(diff.genes, Gene %in% matches)


##### Olah gene set (middle age vs old age) #####

olah <- read_xlsx(file.path(dir_external, "Olah_MidvsOldAge_DEG.xlsx"), sheet = 1)

olah$Direction <- sign(olah$logFC)
olah$Direction[olah$adj.P.Val > 0.05] <- 0
olah$Ensembl.Id <- geneNameToEnsembl(olah$ID)
olah$Gene.name <- olah$ID

gene.sets[["Human Mg Aging (Olah 2018)"]] <- olah

olah <- subset(olah, adj.P.Val <= 0.05) %>% arrange(desc(logFC))

ou <- subset(olah, logFC > 0)
od <- subset(olah, logFC < 0)

matches <- intersect(ou$ID, diff.genes$Gene)
sheets[["Olah up"]] <- subset(diff.genes, Gene %in% matches)

matches <- intersect(od$ID, diff.genes$Gene)
sheets[["Olah down"]] <- subset(diff.genes, Gene %in% matches)


##### CellAge Senescence genes #####

cellage <- read_xlsx(file_cellage_edited)
cellage$Direction <- 0
cellage$Direction[cellage$senescence_effect == "Induces"] <- 1
cellage$Direction[cellage$senescence_effect == "Inhibits"] <- -1
cellage$Ensembl.Id <- geneNameToEnsembl(cellage$gene_name)
cellage$Gene.name <- cellage$gene_name

gene.sets[["Senescence Genes (CellAge)"]] <- cellage

cu <- subset(cellage, senescence_effect == "Induces")
cd <- subset(cellage, senescence_effect == "Inhibits")

matches <- intersect(cu$gene_name, diff.genes$Gene)
sheets[["CellAge Induces"]] <- subset(diff.genes, Gene %in% matches)

matches <- intersect(cd$gene_name, diff.genes$Gene)
sheets[["CellAge Inhibits"]] <- subset(diff.genes, Gene %in% matches)


##### GenAge brain aging-related genes #####

genage <- list("Up" = read_xlsx(file_genage_edited, sheet = "Genes_overexpressed"),
               "Down" = read_xlsx(file_genage_edited, sheet = "Genes_underexpressed"))

genage.df <- rbind(genage[["Up"]], genage[["Down"]])
genage.df$Direction <- sign(genage.df$Mouse_brain)
valid <- !(is.na(genage.df$Human_brain))
genage.df$Direction[valid] <- sign(genage.df$Human_brain[valid]) # Use human data where possible
genage.df$Ensembl.Id <- geneNameToEnsembl(genage.df$Symbol)
genage.df$Gene.name <- genage.df$Symbol

gene.sets[["Brain Aging (GenAge Meta-Analysis)"]] <- genage.df

matches <- intersect(genage[["Up"]]$Symbol, diff.genes$Gene)
sheets[["GenAge Up"]] <- subset(diff.genes, Gene %in% matches)

matches <- intersect(genage[["Down"]]$Symbol, diff.genes$Gene)
sheets[["GenAge Down"]] <- subset(diff.genes, Gene %in% matches)


##### HAM genes #####

ham <- list("Ham up" = read_xlsx(file.path(dir_external, "Srinivasan_HAM_Genes.xlsx"), sheet = "HAM Up"),
            "Ham down" = read_xlsx(file.path(dir_external, "Srinivasan_HAM_Genes.xlsx"), sheet = "HAM Down"))

ham[["Ham up"]]$Direction <- 1
ham[["Ham down"]]$Direction <- -1

ham.df <- rbind(ham[["Ham up"]], ham[["Ham down"]])
ham.df$Ensembl.Id <- geneNameToEnsembl(ham.df$Gene)
ham.df$Gene.name <- ham.df$Gene

gene.sets[["HAM Genes (Srinivasan 2020)"]] <- ham.df

matches <- intersect(ham[["Ham up"]]$Gene, diff.genes$Gene)
sheets[["HAM up"]] <- subset(diff.genes, Gene %in% matches)

matches <- intersect(ham[["Ham down"]]$Gene, diff.genes$Gene)
sheets[["HAM down"]] <- subset(diff.genes, Gene %in% matches)


##### DAM genes #####

dam <- list("Dam up" = read_xlsx(file.path(dir_external, "KerenShaul_DAM_Homologs.xlsx"), sheet = "DAM Up"),
            "Dam down" = read_xlsx(file.path(dir_external, "KerenShaul_DAM_Homologs.xlsx"), sheet = "DAM Down"))

dam[["Dam up"]]$Direction <- 1
dam[["Dam down"]]$Direction <- -1
dam.df <- rbind(dam[["Dam up"]], dam[["Dam down"]])
dam.df$Ensembl.Id <- dam.df$Human.gene.stable.ID
dam.df$Gene.name <- dam.df$Human.gene.name

gene.sets[["DAM genes (Keren-Shaul 2017)"]] <- dam.df

genes <- unique(c(dam[["Dam up"]]$Human.gene.name, ensemblToGeneName(dam[["Dam up"]]$Human.gene.stable.ID)))

matches <- intersect(genes, diff.genes$Gene)
sheets[["DAM up"]] <- subset(diff.genes, Gene %in% matches)

genes <- unique(c(dam[["Dam down"]]$Human.gene.name, ensemblToGeneName(dam[["Dam down"]]$Human.gene.stable.ID)))

matches <- intersect(genes, diff.genes$Gene)
sheets[["DAM down"]] <- subset(diff.genes, Gene %in% matches)


##### Holtman set (aging and Ercc1 mice) #####

holtman <- lapply(c("Aged Up", "Aged Down", "Ercc1 Up", "Ercc1 Down"), function (S) {
  read_xlsx(file.path(dir_external, "Holtman_Aging_Homologs.xlsx"), sheet = S)
})
names(holtman) <- c("Aged Up", "Aged Down", "Ercc1 Up", "Ercc1 Down")

holtman[["Aged Up"]]$Direction <- 1
holtman[["Ercc1 Up"]]$Direction <- 1
holtman[["Aged Down"]]$Direction <- -1
holtman[["Ercc1 Down"]]$Direction <- -1

for (N in names(holtman)) {
  holtman[[N]]$Ensembl.Id <- holtman[[N]]$Human.gene.stable.ID
  holtman[[N]]$Gene.name <- holtman[[N]]$Human.gene.name
}

ha <- rbind(holtman[["Aged Up"]], holtman[["Aged Down"]])
he <- rbind(holtman[["Ercc1 Up"]], holtman[["Ercc1 Down"]])

gene.sets[["Mouse Mg Aging (Holtman 2015)"]] <- ha
gene.sets[["Mouse Mg Ercc1 \u0394/- (Holtman 2015)"]] <- he

genes <- unique(c(holtman[["Aged Up"]]$Human.gene.name, ensemblToGeneName(holtman[["Aged Up"]]$Human.gene.stable.ID)))
matches <- intersect(genes, diff.genes$Gene)
sheets[["Holtman aged up"]] <- subset(diff.genes, Gene %in% matches)

genes <- unique(c(holtman[["Aged Down"]]$Human.gene.name, ensemblToGeneName(holtman[["Aged Down"]]$Human.gene.stable.ID)))
matches <- intersect(genes, diff.genes$Gene)
sheets[["Holtman aged down"]] <- subset(diff.genes, Gene %in% matches)

genes <- unique(c(holtman[["Ercc1 Up"]]$Human.gene.name, ensemblToGeneName(holtman[["Ercc1 Up"]]$Human.gene.stable.ID)))
matches <- intersect(genes, diff.genes$Gene)
sheets[["Holtman Ercc1 up"]] <- subset(diff.genes, Gene %in% matches)

genes <- unique(c(holtman[["Ercc1 Down"]]$Human.gene.name, ensemblToGeneName(holtman[["Ercc1 Down"]]$Human.gene.stable.ID)))
matches <- intersect(genes, diff.genes$Gene)
sheets[["Holtman Ercc1 down"]] <- subset(diff.genes, Gene %in% matches)

# Save it all to a file
write_xlsx(sheets, file.path(dir_output, "Geneset_Comparisons.xlsx"))


##### Combine gene sets into one matrix #####

scRNA <- readRDS(file_seurat_analyzed)

assay <- GetAssayData(scRNA, slot = "counts", assay = "RNA")
meta <- scRNA@meta.data

# Remove genes expressed in fewer than 10 cells
ok <- rowSums(assay >= 1) >= 10
assay <- assay[ok,]

bulk <- aggregate.Matrix(t(assay), groupings = meta$orig.ident, fun = "sum")
bulk <- edgeR::cpm(t(bulk))
bulk <- aggregate.Matrix(t(bulk), groupings = c("KO", "KO", "KO", "WT", "WT", "WT"), fun = "mean")
direction <- sign(bulk["KO",] - bulk["WT",])

sig <- names(direction) %in% diff.genes$Gene
direction[sig] <- direction[sig] * 2

combined <- data.frame(ERCC1.KO.Mg = direction, row.names = geneNameToEnsembl(names(direction)))
colnames(combined)[1] <- "ERCC1 KO Mg"

for (N in names(gene.sets)) {
  gs <- gene.sets[[N]]
  
  # Account for multiple mappings between Ensembl ID and gene name
  bad <- !(gs$Ensembl.Id %in% rownames(combined))
  gs$Ensembl.Id[bad] <- geneNameToEnsembl(gs$Gene.name[bad])
  
  ok <- which(gs$Ensembl.Id %in% rownames(combined))
  
  gs <- gs[ok,]
  
  combined[, N] <- 0
  combined[gs$Ensembl.Id, N] <- gs$Direction * 2
}

# Only keep genes that are significant in at least one external gene set 
# (in addition to ERCC1 KO data)
keep <- rowSums(combined != 0) > 1
combined <- combined[keep,]
rownames(combined) <- ensemblToGeneName(rownames(combined))

library(pheatmap)

# Lopes + Galatro (Human aging)
tmp <- combined[, c(1, 2, 3)]
keep <- rowSums(tmp != 0) > 1
tmp <- tmp[keep,]

colors <- c("#0571B0", "#92C5De", "#F7F7F7", "#F4A582", "#CA0020")
plt <- pheatmap(t(tmp), treeheight_col = 0, treeheight_row = 0, 
                show_colnames = FALSE, cluster_rows = FALSE,
                color = colors, fontsize = 6, border_color = NA,
                filename = file.path(dir_figures, "human_aging_heatmap_allgenes.png"),
                width = 5, height = 1)
dev.off()
plt

tmp <- tmp[rownames(tmp) %in% diff.genes$Gene,]
plt <- pheatmap(t(tmp), treeheight_col = 0, treeheight_row = 0,
                show_colnames = TRUE, cluster_rows = FALSE, border_color = NA,
                color = colors, fontsize = 6, fontsize_col = 2,
                filename = file.path(dir_figures, "human_aging_heatmap_significant.png"),
                width = 7, height = 1)
dev.off()
plt

# Senescence 
tmp <- combined[, c(1, 5)]
keep <- rowSums(tmp != 0) > 1
tmp <- tmp[keep,]

colors <- c("#0571B0", "#92C5De", "#F7F7F7", "#F4A582", "#CA0020")
plt <- pheatmap(t(tmp), treeheight_col = 0, treeheight_row = 0, 
                show_colnames = FALSE, cluster_rows = FALSE,
                color = colors, fontsize = 6, border_color = NA,
                filename = file.path(dir_figures, "senescence_heatmap_allgenes.png"),
                width = 5, height = 0.6)
dev.off()
plt

tmp <- tmp[rownames(tmp) %in% diff.genes$Gene,]
plt <- pheatmap(t(tmp), treeheight_col = 0, treeheight_row = 0,
                show_colnames = TRUE, cluster_rows = FALSE, border_color = NA,
                color = colors, fontsize = 6, fontsize_col = 4,
                filename = file.path(dir_figures, "senescence_heatmap_significant.png"),
                width = 4, height = 1)
dev.off()
plt

# GenAge 
tmp <- combined[, c(1, 6)]
keep <- rowSums(tmp != 0) > 1
tmp <- tmp[keep,]

colors <- c("#0571B0", "#92C5De", "#F7F7F7", "#F4A582", "#CA0020")
plt <- pheatmap(t(tmp), treeheight_col = 0, treeheight_row = 0, 
                show_colnames = FALSE, cluster_rows = FALSE,
                color = colors, fontsize = 6, border_color = NA,
                filename = file.path(dir_figures, "genAge_heatmap_allgenes.png"),
                width = 5, height = 0.6)
dev.off()
plt

tmp <- tmp[rownames(tmp) %in% diff.genes$Gene,]
plt <- pheatmap(t(tmp), treeheight_col = 0, treeheight_row = 0,
                show_colnames = TRUE, cluster_rows = FALSE, border_color = NA,
                color = colors, fontsize = 6, fontsize_col = 4,
                filename = file.path(dir_figures, "genAge_heatmap_significant.png"),
                width = 3, height = 1)
dev.off()
plt

# DAM genes
tmp <- combined[, c(1, 8)]
keep <- rowSums(tmp != 0) > 1
tmp <- tmp[keep,]

colors <- c("#0571B0", "#92C5De", "#F7F7F7", "#F4A582", "#CA0020")
plt <- pheatmap(t(tmp), treeheight_col = 0, treeheight_row = 0, 
                show_colnames = FALSE, cluster_rows = FALSE,
                color = colors, fontsize = 6, border_color = NA,
                filename = file.path(dir_figures, "dam_heatmap_allgenes.png"),
                width = 5, height = 1)
dev.off()
plt

tmp <- tmp[rownames(tmp) %in% diff.genes$Gene,]
plt <- pheatmap(t(tmp), treeheight_col = 0, treeheight_row = 0,
                show_colnames = TRUE, cluster_rows = FALSE, border_color = NA,
                color = colors, fontsize = 6, fontsize_col = 3,
                filename = file.path(dir_figures, "dam_heatmap_significant.png"),
                width = 5, height = 1)
dev.off()
plt

# Holtman
tmp <- combined[, c(1, 9, 10)]
keep <- rowSums(tmp != 0) > 1
tmp <- tmp[keep,]

colors <- c("#0571B0", "#92C5De", "#F7F7F7", "#F4A582", "#CA0020")
plt <- pheatmap(t(tmp), treeheight_col = 0, treeheight_row = 0, 
                show_colnames = FALSE, cluster_rows = FALSE,
                color = colors, fontsize = 6, border_color = NA,
                filename = file.path(dir_figures, "mouse_aging_heatmap_allgenes.png"),
                width = 7, height = 1)
dev.off()
plt

tmp <- tmp[rownames(tmp) %in% diff.genes$Gene,]
plt <- pheatmap(t(tmp), treeheight_col = 0, treeheight_row = 0,
                show_colnames = TRUE, cluster_rows = FALSE, border_color = NA,
                color = colors, fontsize = 6, fontsize_col = 3,
                filename = file.path(dir_figures, "mouse_aging_heatmap_significant.png"),
                width = 7.8, height = 1)
dev.off()
plt

# Mouse + human aging
tmp <- combined[, c(1, 2, 9)]
keep <- rowSums(tmp != 0) > 1
tmp <- tmp[keep,]

colors <- c("#0571B0", "#92C5De", "#F7F7F7", "#F4A582", "#CA0020")
plt <- pheatmap(t(tmp), treeheight_col = 0, treeheight_row = 0, 
                show_colnames = FALSE, cluster_rows = FALSE,
                color = colors, fontsize = 6, border_color = NA,
                filename = file.path(dir_figures, "human_mouse_aging_heatmap_allgenes.png"),
                width = 7, height = 1)
dev.off()
plt

tmp <- tmp[rownames(tmp) %in% diff.genes$Gene,]
plt <- pheatmap(t(tmp), treeheight_col = 0, treeheight_row = 0,
                show_colnames = TRUE, cluster_rows = FALSE, border_color = NA,
                color = colors, fontsize = 6, fontsize_col = 3,
                filename = file.path(dir_figures, "human_mouse_aging_heatmap_significant.png"),
                width = 10, height = 1)
dev.off()
plt

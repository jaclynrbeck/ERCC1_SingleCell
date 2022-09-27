# Helper functions that are used in the analysis step

# Author: Jaclyn Beck
# Copied from https://github.com/jaclynrbeck/TCellsAD2022, modified for human
# genes and clonotype functions removed. 

library(readxl)
library(writexl)
library(biomaRt)


##### File I/O #####

# Writes output of FindAllMarkers to an Excel file. Automatically replaces
# gene names that have been altered with "--1" at the end, with the actual
# gene name. Adds the Ensembl ID of the gene as a column.
# sig.markers should be a data frame as output by FindAllMarkers.
# out.file should be a string with the full file path of the output file.
writeDifferentialGenes <- function( sig.markers, out.file ) {
  sheets <- c()
  
  sig.markers$Ensembl.ID <- geneNameToEnsembl(sig.markers$gene)
  sig.markers$Gene <- str_replace(sig.markers$gene, "--.*", "")
  
  for (i in levels(sig.markers$cluster)) {
    sub <- filter(sig.markers, cluster==i)
    sub <- sub[order(sub$avg_log2FC, decreasing=TRUE),]
    
    # cut out "p_val", "cluster", and old "gene" columns
    sub <- sub[, c("Ensembl.ID", "Gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
    
    # Some clusters have a "/" in their names (i.e. Naive/CM), which is an
    # illegal character for sheet names. Replace it with "+".
    sheets[[str_replace(i, "/", "+")]] <- sub
  }
  
  write_xlsx(sheets, path=out.file)
}


# Calculates differential genes between genotypes, running comparisons on 
# <genotype> vs all other genotypes, and <genotype> vs each individual genotype.
# Automatically replaces gene names that have been altered with "--1" at the 
# end, with the actual gene name. Adds the Ensembl ID of the gene as a column.
# scRNA: Seurat object
# genotypes: list or vector of genotype names, corresponding to what is in
#            scRNA$genotype. i.e. c("WT", "5XFAD", "PS19", "PS-5X")
# FDR: threshold for significance, i.e. 0.01
# out.file: string with full file path to output file
writeGenotypeDifferentialGenes <- function( scRNA, genotypes, FDR, out.file ) {
  sheets <- c()
  
  for (G1 in genotypes) {
    # G1 vs All
    markers <- FindMarkers(scRNA, ident.1 = G1, group.by = "genotype",
                           test.use = "MAST", min.pct = 0.05,
                           latent.vars = c("nCount_RNA", "StressScoreVariable"))
    
    sig.markers <- filter(markers, p_val_adj <= FDR & pct.1 >= 0.05)
    sig.markers$Ensembl.ID <- geneNameToEnsembl(rownames(sig.markers))
    sig.markers$Gene <- str_replace(rownames(sig.markers), "--.*", "")
    
    sub <- sig.markers[order(sig.markers$avg_log2FC, decreasing=TRUE),]
    
    # cut out "p_val"
    sub <- sub[, c("Ensembl.ID", "Gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
    
    sheets[[paste(G1, "vs All")]] <- sub
    
    # G1 vs G2 for each G2 in <genotypes>
    for (G2 in genotypes[genotypes != G1]) {
      markers <- FindMarkers(scRNA, ident.1 = G2, ident.2 = G1, min.pct = 0.05,
                             group.by = "genotype", test.use = "MAST",
                             latent.vars = c("nCount_RNA", "StressScoreVariable"))
      
      sig.markers <- filter(markers, p_val_adj <= FDR & pct.1 >= 0.05)
      sig.markers$Ensembl.ID <- geneNameToEnsembl(rownames(sig.markers))
      sig.markers$Gene <- str_replace(rownames(sig.markers), "--.*", "")
      
      sub <- sig.markers[order(sig.markers$avg_log2FC, decreasing=TRUE),]
      
      # cut out "p_val"
      sub <- sub[, c("Ensembl.ID", "Gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
      
      sheets[[paste(G2, "vs", G1)]] <- sub
    }
    
    write_xlsx(sheets, path=out.file)
  }
}


# Calculates differential genes between genotypes, on a per-cluster basis. 
# Only runs <genotype> vs All comparison. Automatically accounts for gene names
# with "--1" at the end. 
# scRNA: Seurat object
# genotypes: list or vector of genotype names, corresponding to what is in
#            scRNA$genotype. i.e. c("WT", "5XFAD", "PS19", "PS-5X")
# FDR: threshold for significance, i.e. 0.01
# out.file: string with full file path to output file
writeClusterVGenotypeDiffGenes <- function( scRNA, genotypes, FDR, out.file ) {
  sheets <- c()
  clust <- unique(Idents(scRNA))
  
  for (C in clust) {
    for (G in genotypes) {
      markers <- FindMarkers(scRNA, ident.1 = G, subset.ident = C, 
                             group.by = "genotype", 
                             test.use = "MAST",
                             latent.vars = c("nCount_RNA", "StressScoreVariable"))
      
      sig.markers <- filter(markers, p_val_adj <= FDR & pct.1 >= 0.1)
      sig.markers$Ensembl.ID <- geneNameToEnsembl(rownames(sig.markers))
      sig.markers$Gene <- str_replace(rownames(sig.markers), "--.*", "")
      
      sub <- sig.markers[order(sig.markers$avg_log2FC, decreasing=TRUE),]
      
      # cut out "p_val"
      sub <- sub[, c("Ensembl.ID", "Gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
      
      sheets[[paste(G, "vs All -", C)]] <- sub
    }
  }
  
  write_xlsx(sheets, path=out.file)
}


##### Things that print to the console #####

# Prints the top 10 markers for each cluster
# sig.markers: data frame output by FindAllMarkers
# pos.only: TRUE = display only positively-changed genes. 
#           FALSE = display both positively and negatively changed genes
printTop10Markers <- function ( sig.markers, pos.only = TRUE ) {
  if (pos.only) {
    sorted <- sig.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  }
  else {
    sorted <- sig.markers %>% group_by(cluster) %>% top_n(n = 10, wt = abs(avg_log2FC))
  }
  visual_sorted <- data.frame(Cluster=unique(sorted$cluster),
                              Genes=" ")
  genes <- sapply(visual_sorted$Cluster, 
                  function(x,y) 
                    list(subset(y, cluster==x)$gene),
                  sorted)
  visual_sorted$Genes <- genes
  
  tmp <- data.frame(visual_sorted$Genes[1])
  for (r in 2:nrow(visual_sorted)) {
    tmp <- cbind(tmp, visual_sorted$Genes[r])
  }
  colnames(tmp) <- paste("Cluster", 1:nrow(visual_sorted)-1)
  print(tmp)
  tmp
}


# Prints the percentage of each genotype in each cluster, and the percentage
# of each cluster in each genotype. 
printClusterDistributions <- function(scRNA) {
  summ_ident_table = table(scRNA$orig.ident, Idents(scRNA))
  
  # Distribution across clusters of each group separately (sum of columns = 100)
  clust_per_group = round(t(summ_ident_table / rowSums(summ_ident_table))*100, 2)
  print("Dist. of each cluster within each genotype (columns sum to 100)")
  print(clust_per_group)
  
  # Distribution of groups in each cluster (sum of rows = 100)
  print("")
  print("Dist. of each genotype within each cluster (rows sum to 100)")
  print(round((clust_per_group / rowSums(clust_per_group))*100, 2))
}


# Gets the 10 highest expressed genes (by scaled count) in each cluster.
# Uses the "integrated" assay. 
getHighestExpressedGenes <- function( scRNA, ngenes = 10 ) {
  DefaultAssay(scRNA) <- "integrated"
  highest.expressed = data.frame()
  for (i in levels(scRNA)) {
    sub <- subset(scRNA, idents=i)
    data <- GetAssayData(object = sub, slot = "scale.data")
    top.genes <- sort(rowMeans(data), decreasing=TRUE)[1:ngenes]
    if (length(highest.expressed) == 0) {
      highest.expressed <- data.frame("0"=names(top.genes))
    }
    else {
      highest.expressed <- cbind(highest.expressed, names(top.genes))
    }
  }
  
  colnames(highest.expressed) <- paste("Cluster", levels(scRNA), sep=" ")
  DefaultAssay(scRNA) <- "RNA"
  highest.expressed
}


##### BioMart Queries #####

# Queries BioMart to get mouse homologs for human genes
# human.genes: vector of gene names
getHumanMouseHomologs <- function(human.genes) {
  # Temporarily using archive because new ensembl update crashes biomart
  ens.human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
                       host = "https://dec2021.archive.ensembl.org/")
  ens.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", 
                       host = "https://dec2021.archive.ensembl.org/")
  
  homologs <- getLDS(attributes = c('hgnc_symbol'),
                     filters = 'hgnc_symbol', 
                     values = unique(human.genes), 
                     mart = ens.human,
                     attributesL = c('mgi_symbol'), 
                     martL = ens.mouse)
}


getOfficialMouseGenes <- function(mouse.genes) {
  ens.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", 
                       host = "https://dec2021.archive.ensembl.org/")
  
  # Get official symbol names, accounting for synonyms
  tmp1 <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                filters = "external_gene_name",
                values = mouse.genes,
                mart = ens.mouse)
  tmp1$external_synonym <- ""
  
  unknown <- setdiff(mouse.genes, tmp1$external_gene_name)
  tmp2 <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "external_synonym"),
                filters = "external_synonym",
                values = unknown,
                mart = ens.mouse)
  
  mouse.genes.official <- rbind(tmp1, tmp2) %>% unique()
}


getMouseHumanHomologs <- function(mouse.genes) {
  # Temporarily using archive because new ensembl update crashes biomart
  ens.human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
                       host = "https://dec2021.archive.ensembl.org/")
  ens.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", 
                       host = "https://dec2021.archive.ensembl.org/")
  
  homologs <- getLDS(attributes = c('external_gene_name', 'ensembl_gene_id'),
                     filters = c('external_gene_name'), 
                     values = unique(mouse.genes), 
                     mart = ens.mouse,
                     attributesL = c('external_gene_name', 'ensembl_gene_id'), 
                     martL = ens.human)
  colnames(homologs) <- c("Mouse.gene.name", "Mouse.gene.stable.ID", 
                          "Human.gene.name", "Human.gene.stable.ID")
  return(homologs)
}


##### GO analysis functions #####

# Runs GProfiler to get statisticaly significant GO terms.
# markers: data frame of significant markers as output by FindMarkers
# sources: list of sources compatible with the "sources" argument of gost.
#          See "?gost" for possible values. 
runGOAnalysis <- function( markers, sources = c("GO", "REAC") ) {
  m.up <- subset(markers, avg_log2FC > 0)
  m.down <- subset(markers, avg_log2FC < 0)
  
  # Order by -log10(p value) : most significant first
  m.up <- m.up[order(-log10(m.up$p_val_adj), decreasing = TRUE),]
  m.down <- m.down[order(-log10(m.down$p_val_adj), decreasing = TRUE),]
  
  go.res.up <- list()
  go.res.down <- list()
  
  if (nrow(m.up) > 0) {
    go.res.up <- gost(query = m.up$Ensembl.ID, 
                      organism = "mmusculus", ordered_query = TRUE, 
                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                      measure_underrepresentation = FALSE, evcodes = TRUE, 
                      user_threshold = 0.05, correction_method = "g_SCS", 
                      domain_scope = "annotated", custom_bg = NULL, 
                      numeric_ns = "", sources = sources, 
                      as_short_link = FALSE)
  }
  
  if (nrow(m.down) > 0) {
    go.res.down <- gost(query = m.down$Ensembl.ID, 
                        organism = "mmusculus", ordered_query = TRUE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                        measure_underrepresentation = FALSE, evcodes = TRUE, 
                        user_threshold = 0.05, correction_method = "g_SCS", 
                        domain_scope = "annotated", custom_bg = NULL, 
                        numeric_ns = "", sources = sources, 
                        as_short_link = FALSE)
  }
  
  if (!is.null(go.res.up)) {
    go.res.up$result$Phenotype <- "+1"
  }
  if (!is.null(go.res.down)) {
    go.res.down$result$Phenotype <- "-1"
  }
  
  return(list("Up" = go.res.up, "Down" = go.res.down))
}


# Converts gost output to a Gem file for use with Cytoscape
# go.res.list: a list with data frames output by gost. Items in list should
#              be named "Up" and "Down"
# outfile: string with the full file path to the output gem file
# use.ensembl: TRUE = use Ensembl IDs in the gem file. (recommended)
#              FALSE = use gene names in the gem file.
goResultsToGem <- function( go.res.list, outfile, use.ensembl = TRUE ) {
  go.res.up <- go.res.list[["Up"]]
  go.res.down <- go.res.list[["Down"]]
  
  go.res <- NULL
  
  if (!is.null(go.res.up) & !is.null(go.res.down)) {
    go.res <- rbind(go.res.up$result, go.res.down$result)
  }
  else if (!is.null(go.res.up)) {
    go.res <- go.res.up
  }
  else if (!is.null(go.res.down)) {
    go.res <- go.res.down
  }
  
  if (!is.null(go.res)) {
    go.res <- subset(go.res, term_size <= 1000)
    gem <- go.res[,c("term_id", "term_name", "p_value", "intersection", "Phenotype")]  
    colnames(gem) <- c("GO.ID", "Description", "p.Val", "Genes", "Phenotype")
    gem$FDR <- gem$p.Val
    gem <- gem[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
    
    if (!use.ensembl) {
      tmp <- sapply(str_split(gem$Genes, pattern = ","), ensemblToGeneName)
      gem$Genes <- sapply(tmp, str_c, collapse = ",")
    }
    
    write.table(gem, file = outfile, sep = "\t", quote = F, row.names = F)
  }
}

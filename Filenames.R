# A list of filenames used in multiple places in the pipeline, to avoid mistakes
# or using the wrong file. 

# Main data directory
dir_data <- "data"


# Raw data from CellRanger
dir_cellranger <- file.path(dir_data, "CellRanger")
dir_filtered_counts <- file.path(dir_cellranger, "count", "filtered_feature_bc_matrix")


# Gene symbols / ensembl IDs files
file_gene_symbols <- file.path(dir_data, "gene_symbols.rds")


# Un-normalized and normalized Seurat objects, and analyzed objects
dir_seurat <- file.path(dir_data, "Seurat")
file_seurat_unnorm <- file.path(dir_seurat, "seurat_unnorm_2022-09-26.rds")
file_seurat_norm <- file.path(dir_seurat, "seurat_liger_2022-09-26.rds")
file_seurat_analyzed <- file.path(dir_seurat, "seurat_analyzed_2022-09-26.rds")


# Cluster / Genotype markers
dir_output <- file.path(dir_data, "Output")
dir_go <- file.path(dir_output, "GO_Analysis")
dir_dpa <- file.path(dir_output, "DPA")
file_markers_all <- file.path(dir_output, "all_markers_2022-09-26.rds")
file_markers_clusters <- file.path(dir_output, 'DiffGenes_byCluster_2022-09-26.xlsx')
file_markers_genotypes <- file.path(dir_output, 'DiffGenes_byGenotype_2022-09-26.xlsx')
file_markers_clusters_vs_genotype <- file.path(dir_output, 'DiffGenes_byGenotypePerCluster_2022-09-26.xlsx')


# Figures
dir_figures <- file.path("figures")


# Reference gene sets
dir_external <- file.path(dir_data, "external")
dir_gene_sets_unedited <- file.path(dir_external, "originals")
file_holtman_original <- file.path(dir_gene_sets_unedited, "Holtman_Aging_DEG.xlsx")
file_holtman_homologs <- file.path(dir_external, "Holtman_Aging_Homologs.xlsx")
file_ks_original <- file.path(dir_gene_sets_unedited, "KerenShaul_DAM2vsHomeostatic_DEG.xlsx")
file_ks_homologs <- file.path(dir_external, "KerenShaul_DAM_Homologs.xlsx")
file_cellage_original <- file.path(dir_gene_sets_unedited, "CellAge", "CellAgeDatabase.csv")
file_cellage_edited <- file.path(dir_external, "CellAgeDatabase_Clean.xlsx")
file_genage_original <- file.path(dir_gene_sets_unedited, "CellAge", "GenAge_Meta-analysis_tables.xls")
file_genage_edited <- file.path(dir_external, "GenAge_Meta-analysis_Clean.xlsx")


# Ensure all directories we plan to read/write to exist. 
dirs <- grep("dir_", ls(), value = TRUE)
for (D in dirs) {
  dir_eval <- eval(parse(text = D))
  if (!file.exists(dir_eval)) {
    dir.create(dir_eval, recursive = TRUE)
    print(paste0("Created directory ", dir_eval))
  }
}

rm(D)



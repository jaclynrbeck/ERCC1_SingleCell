# Output from CellRanger aggr doesn't carry genotype over into its barcodes. 
# It only increases the number at the end of the barcode after each sample it 
# adds to the aggregated data.
#
# This script re-writes the "barcodes.tsv.gz" file with the genotype embedded in
# the barcodes for easy reading into Seurat. 
#
# Lastly, gene names in the "features.tsv.gz" file aren't necessarily unique,
# which leads Seurat to ignore duplicate names. This is fixed in this script as
# well by adding "--1" to duplicated gene names and re-writing the file. 

# Author: Jaclyn Beck
# Copied from https://github.com/jaclynrbeck/TCellsAD2022, clonotype pieces
# removed. 

library(stringr)

source("Filenames.R")

file_barcodes <- file.path(dir_filtered_counts, "barcodes.tsv.gz")
file_features <- file.path(dir_filtered_counts, "features.tsv.gz")

# Copy the original barcodes and features files as a backup
file.copy(from = file_barcodes, 
          to = file.path(dir_filtered_counts, "barcodes.original.tsv.gz"))

file.copy(from = file_features, 
          to = file.path(dir_filtered_counts, "features.original.tsv.gz"))

# Read in barcodes and annotation
barcodes <- read.table(file = gzfile(file_barcodes))
barcodes <- barcodes$V1

anno <- read.csv(file = file.path(dir_cellranger, "aggregation.csv"))
samples <- anno$sample_id

# Replace the number at the end of each barcode with the sample ID. Also change 
# the "-" separator to "_" so the sample names are read correctly by Seurat
for (N in 1:length(samples)) {
  barcodes <- str_replace(barcodes, paste0("-", N), paste0("_", samples[N]))
}

# Write a new 'barcodes.tsv.gz' file for Read10X
write.table(barcodes, 
            file=gzfile(file_barcodes), sep=' ', row.names=FALSE, 
            col.names=FALSE, quote=FALSE)


# Fix the features file
features <- read.table(file = gzfile(file_features))
features$V2 <- make.unique(features$V2, sep = "--")

write.table(features, 
            file=gzfile(file_features), sep='\t', row.names=FALSE, 
            col.names=FALSE, quote=FALSE)

features <- features[,1:2]
colnames(features) <- c("Ensembl.Id", "Gene.Symbol")
saveRDS(features, file_gene_symbols)

# Clear data
rm(list=ls())
gc()

# Done. 

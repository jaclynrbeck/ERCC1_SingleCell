# Steps 04 and 05 are for comparing differential genes from ERCC1 cells to
# known data sets.
# Some supplementary tables from other papers need editing for easier reading
# into R or contain mouse genes that need to be converted to human genes. 

library(dplyr)
library(readxl)
library(writexl)
library(stringr)

source(file.path("functions", "Analysis_HelperFunctions.R"))
source("Filenames.R")


##### Holtman data set (Aged and Ercc1 delta mice) #####

h_orig <- read_xlsx(file_holtman_original, sheet = "DE_ALL_Datasets")
colnames(h_orig)[1] <- "Gene"

# Throw out everything but the aged set and the Ercc1 set
aged <- subset(h_orig, AGED_FDR_p <= 0.05) %>% arrange(desc(AGED_logFC))
ercc1 <- subset(h_orig, ERCC1_FDR_p <= 0.05) %>% arrange(desc(ERCC1_logFC))

# Get rid of capitalization, these are mouse genes and not human genes
aged$Gene <- str_to_title(aged$Gene)
ercc1$Gene <- str_to_title(ercc1$Gene)

aged.genes.up <- getOfficialMouseGenes(aged$Gene[aged$AGED_logFC > 0])
aged.genes.down <- getOfficialMouseGenes(aged$Gene[aged$AGED_logFC < 0])
ercc1.genes.up <- getOfficialMouseGenes(ercc1$Gene[ercc1$ERCC1_logFC > 0])
ercc1.genes.down <- getOfficialMouseGenes(ercc1$Gene[ercc1$ERCC1_logFC < 0])

agedUp.hom <- getMouseHumanHomologs(aged.genes.up$external_gene_name)
agedDown.hom <- getMouseHumanHomologs(aged.genes.down$external_gene_name)
ercc1Up.hom <- getMouseHumanHomologs(ercc1.genes.up$external_gene_name)
ercc1Down.hom <- getMouseHumanHomologs(ercc1.genes.down$external_gene_name)

homologs <- list("Aged Up" = agedUp.hom, 
                 "Aged Down" = agedDown.hom,
                 "Ercc1 Up" = ercc1Up.hom,
                 "Ercc1 Down" = ercc1Down.hom)

write_xlsx(homologs, path = file_holtman_homologs)

# A few homologous genes aren't coming up, so these will be manually input 
missing <- setdiff(aged.genes.up$external_gene_name, agedUp.hom$Mouse.gene.name)
write.table(missing, file = file.path(dir_external, "missing_holtmanAgedUp.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
missing <- setdiff(aged.genes.down$external_gene_name, agedDown.hom$Mouse.gene.name)
write.table(missing, file = file.path(dir_external, "missing_holtmanAgedDown.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
missing <- setdiff(ercc1.genes.up$external_gene_name, ercc1Up.hom$Mouse.gene.name)
write.table(missing, file = file.path(dir_external, "missing_holtmanErcc1Up.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
missing <- setdiff(ercc1.genes.down$external_gene_name, ercc1Down.hom$Mouse.gene.name)
write.table(missing, file = file.path(dir_external, "missing_holtmanErcc1Down.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)


##### Keren-Shaul data set (DAM genes in mice) #####

dam_orig <- read_xlsx(file_ks_original, sheet = "Diff.expression_mic3tomic1")
colnames(dam_orig)[1] <- "Gene"
dam_orig <- arrange(dam_orig, desc("up/down"))

genes.up <- getOfficialMouseGenes(dam_orig$Gene[dam_orig$`up/down` == 1])
genes.down <- getOfficialMouseGenes(dam_orig$Gene[dam_orig$`up/down` == -1])

damUp.hom <- getMouseHumanHomologs(genes.up$external_gene_name)
damDown.hom <- getMouseHumanHomologs(genes.down$external_gene_name)

homologs <- list("DAM Up" = damUp.hom,
                 "DAM Down" = damDown.hom)

write_xlsx(homologs, path = file_ks_homologs)

# A few homologous genes aren't coming up, so these will be manually input 
missing <- setdiff(genes.up$external_gene_name, damUp.hom$Mouse.gene.name)
write.table(missing, file = file.path(dir_external, "missing_damUp.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
missing <- setdiff(genes.down$external_gene_name, damDown.hom$Mouse.gene.name)
write.table(missing, file = file.path(dir_external, "missing_damDown.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)


##### CellAge senescence genes (human, cleaning up format) #####

age <- read.table(file_cellage_original, sep = ";", header = TRUE)
age <- arrange(age, senescence_effect)
write_xlsx(age, file_cellage_edited)


##### GenAge meta-analysis aging genes (human, cleaning up) #####

genage <- list("Up" = read_xls(file_genage_original, sheet = "Genes_overexpressed", skip = 14),
               "Down" = read_xls(file_genage_original, sheet = "Genes_underexpressed", skip = 9))

for (N in c("Up", "Down")) {
  remove <- grep("^Mouse|^Rat", genage[[N]]$Symbol, value = T)
  genage[[N]] <- subset(genage[[N]], !is.na(Symbol) & !(Symbol %in% remove))
  genage[[N]]$Human_brain <- as.numeric(str_replace(genage[[N]]$Human_brain, "\\*", ""))
  genage[[N]]$Mouse_brain <- as.numeric(str_replace(genage[[N]]$Mouse_brain, "\\*", ""))
}

# Only take genes where the human brain data agrees with direction. If no 
# human data available, use mouse brain data
genage[["Up"]] <- subset(genage[["Up"]], Human_brain > 0 | 
                           (is.na(Human_brain) & Mouse_brain > 0))
genage[["Down"]] <- subset(genage[["Down"]], Human_brain < 0 | 
                             (is.na(Human_brain) & Mouse_brain < 0))

names(genage) <- c("Genes_overexpressed", "Genes_underexpressed")

write_xlsx(genage, file_genage_edited)


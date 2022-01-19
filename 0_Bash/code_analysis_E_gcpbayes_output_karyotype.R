# 		 Script for creation of a Karyotype Plot from the GCPBayes output
# ================================================================================
# Summary: 
# creation of a Karyotype Plot from the GCPBayes output 
# Genes with potential pleiotropic effects are considered and showed on a Karyotype Plot based on their position
# NOTE: the script uses the HGNC symbol (for each input gene symbol) extracted from Ensembl database, 
# So, it would NOT be able to show the genes that their HGNC symbol did not find.
# ================================================================================
# Written by: Yazdan Asgari
# Edited Date: 1/2022
# https://cesp.inserm.fr/en/equipe/exposome-and-heredity
# ================================================================================
# for calculation the running time of the program
start_time <- Sys.time()

# ================================================================================
# libraries used in the code
# please install the following packages if needed
#BiocManager::install("biomaRt")
#BiocManager::install("regioneR")
#BiocManager::install("karyoploteR")
library(biomaRt)
library(regioneR)
library(karyoploteR)

# ================================================================================
# 								DEFINITION SECTION
# ================================================================================
# directory in which input data exist
path_input_data = "~/BCAC_OCAC/"

# the file name of input data
gene_list <- "gene_list.txt"

# ================================================================================
# 								RUNNING SECTION
# ================================================================================
# reading the input data
geneSymbols <- read.table(file=paste0(path_input_data, gene_list), header = FALSE)
geneSymbols <- as.character(unlist(geneSymbols))

# getting information from biomaRt
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# finding position of the genes
genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                         filters = 'hgnc_symbol', values =geneSymbols, mart = ensembl))

seqlevelsStyle(genes) <- "UCSC"

head(genes)

# All chromosomes (1-22) in one plot
kp <- plotKaryotype("hg19")
kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, text.orientation = "vertical",
              r1=0.5, cex=0.6, adjust.label.position = TRUE, ignore.chromosome.ends = TRUE ,
              label.color = "blue")

# Showing chromosomes in two plots (1-10 and 11-22)
kp <- plotKaryotype("hg19", plot.type=2, chromosomes = c("chr1","chr2", "chr3", "chr4","chr5", "chr6", 
                                                         "chr7","chr8", "chr9", "chr10"))
kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, text.orientation = "vertical",
              r1=0.5, cex=0.6, adjust.label.position = TRUE, ignore.chromosome.ends = TRUE ,
              label.color = "blue")

kp <- plotKaryotype("hg19", plot.type=2, chromosomes = c("chr11","chr12", "chr13", "chr14","chr15", "chr16", 
                                                         "chr17","chr18", "chr19", "chr20", "chr21", "chr22"))
kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, text.orientation = "vertical",
              r1=0.5, cex=0.6, adjust.label.position = TRUE, ignore.chromosome.ends = TRUE ,
              label.color = "blue")

# ================================================================================
# write the details information of the computer and softwares related to this code 
# please keep this info since it might be needed for any scientist in the future!
sessionInfo()
# calculation the running time of the program
end_time <- Sys.time()
end_time - start_time
# ================================================================================

# 		                Script for creation of Annotation files 
# ================================================================================
# Summary: 
# creation of the coding-gene annotation file for Human organism based on the file downloaded from GENCODE webpage 
# Then, adding a gene column (based on the created annotation file) to a GWAS summary statistics data (for one of the traits)
# For preparation of the input file (annot_gencode_v38lift37_modified_gene_class.txt), please see our GitHub page:
# https://github.com/CESP-ExpHer/Gene_Annotation/tree/main/1_hg19/2_Creation_annot_file 
# ================================================================================
# Written by: Yazdan Asgari
# Initial Creation Date: 07/2021
# Edited Date: 1/2022
# ================================================================================
# for calculation the running time of the program
start_time <- Sys.time()

# ================================================================================
# 						DEFINITION SECTION
# should be changed by a user
# ================================================================================
path_input_data = "~/BCAC_OCAC/"
input_data = "annot_gencode_v38lift37_modified_gene_class.txt"

path_output_data_coding = "~/BCAC_OCAC/"
output_data_coding = "annot_gencode_v38lift37_modified_gene_class_coding.txt"

path_output_data_plink = "~/BCAC_OCAC/"
output_data_plink = "annot_gencode_v38lift37_modified_gene_class_coding_chr_num_plink_input.txt"

# ================================================================================
# libraries used in the code
# please install the following packages if needed
#BiocManager::install("dplyr")
library(dplyr)

# ================================================================================
# 						RUNNING SECTIONS
# ================================================================================
# extraction of "protein_coding" genes from the annotation file
annot_all <- read.table(file = paste0(path_input_data, input_data), header = TRUE)

coding <- filter(annot_all, annot_all$gene_type == "protein_coding")

write.table(coding, file = paste0(path_output_data_coding, output_data_coding), 
            quote = FALSE, col.names = TRUE, sep = "\t")

# changing all chromosomes characters to numbers
# add a chromosome columns which are all numbers
x <- substring(coding[, 1], 4, nchar(coding[, 1]))
coding$chr <- x

coding$chr <- gsub("X", "23", coding$chr)
coding$chr <- gsub("x", "23", coding$chr)
coding$chr <- gsub("Y", "24", coding$chr)
coding$chr <- gsub("y", "24", coding$chr)
coding$chr <- gsub("M", "25", coding$chr)
coding$chr <- gsub("m", "25", coding$chr)

coding$chr <- as.numeric(coding$chr)

# creation of the annotation file in the Plink input format
Annot_plink <- coding[, c("chr", "start", "end", "gene_name")]

write.table(Annot_plink, file = paste0(path_output_data_plink, output_data_plink), sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
			
# ================================================================================
# write the details information of the computer and softwares related to this code 
# please keep this info since it might be needed for any scientist in the future!
sessionInfo()
# calculation the running time of the program
end_time <- Sys.time()
end_time - start_time
# ================================================================================
# 		 Script for creation of an overview table for the GCPBayes output
# ================================================================================
# Summary: 
# creation of a Table from the GCPBayes output 
# Genes with potential pleiotropic effects are considered 
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
#BiocManager::install("readxl")
#BiocManager::install("dplyr")
library(readxl)
library(dplyr)

# ================================================================================
# 								DEFINITION SECTION
# ================================================================================
# directory in which input data exist
path_input_data = "~/BCAC_OCAC/"

# the file name of input data
# It could be extracted from the "...._pleiotropy.txt" file
# The file of genes with potential pleiotropic effects MUST contains the following column names
# chr, gene_length, snp_number
gene_list <- "gene_list.txt"

# the file name of output data
output_table <- "gcpbayes_pleiotropy_summary_table.txt"

# ================================================================================
# 								RUNNING SECTION
# ================================================================================
# reading the input data
GCPBayes_pleio <- read.table(file=paste0(path_input_data, gene_list), header = TRUE)

# selection of desired columns
GCPBayes_pleio_selection <- GCPBayes_pleio %>% select(chr, gene_length, snp_number)

# calculation of the summary
summary_pleio <- GCPBayes_pleio_selection %>%
  group_by(chr) %>%
  summarise(
    Count_by_chr = length(chr),
    Min_Gene_Length = format(min(gene_length), big.mark="," , scientific=FALSE),
    Max_Gene_Length = format(max(gene_length), big.mark=",", scientific=FALSE),
    Min_SNP_number = min(snp_number),
    Max_SNP_number = max(snp_number)
  ) %>%
  arrange(chr)

# saving the data as a Table
write.table(summary_pleio, file=paste0(path_input_data, output_table), col.names = TRUE, 
            row.names = FALSE, quote = FALSE)

# ================================================================================
# write the details information of the computer and softwares related to this code 
# please keep this info since it might be needed for any scientist in the future!
sessionInfo()
# calculation the running time of the program
end_time <- Sys.time()
end_time - start_time
# ================================================================================

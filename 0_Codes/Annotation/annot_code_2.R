# 	Script for changing format of GWAS sum stat file into a PLINK input file 
# ================================================================================
# Summary: 
# Creation of GWAS input file for PLINK based on BCAC_2020 GWAS reformatted file 
# (the file created after the standardization step) 
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
input_data = "BCAC_2020_onco_ALL_reformatted.txt"

path_output_data_plink = "~/BCAC_OCAC/"
output_data_plink = "BCAC_2020_onco_ALL_reformatted.assoc"

# ================================================================================
# libraries used in the code
# please install the following packages if needed
#BiocManager::install("vroom")
library(vroom)

# ================================================================================
# 						RUNNING SECTIONS
# ================================================================================
gwas <- vroom(file=paste0(path_input_data, input_data))
gwas <- as.data.frame(gwas)
gwas <- gwas[, c(2, 1, 3, 8)]
colnames(gwas) <- c("CHR", "SNP", "BP", "P")

vroom_write(gwas, file=paste0(path_output_data_plink, output_data_plink))

# ================================================================================
# write the details information of the computer and softwares related to this code 
# please keep this info since it might be needed for any scientist in the future!
sessionInfo()
# calculation the running time of the program
end_time <- Sys.time()
end_time - start_time
# ================================================================================
# 			reformatting the PLINK Annotation output 
# ================================================================================
# Summary: 
# Reading PLINK Annotation output "plink.annot" 
# Splitting "P-value" and "ANNOT" into two dolumns 
# In the PLINK output, all columns are TAB separated except "P-value" and "ANNOT" columns
# So this step would split them to two columns (TAB separated)
# ================================================================================
# Written by: Yazdan Asgari
# Initial Creation Date: 07/2021
# Edited Date: 12/2021
# ================================================================================
# for calculation the running time of the program
start_time <- Sys.time()

# ================================================================================
# 						DEFINITION SECTION
# should be changed by a user
# ================================================================================
path_input_data = "~/BCAC_OCAC/"
input_data = "plink.annot"

path_output_data = "~/BCAC_OCAC/"
output_data = "Annot_BCAC_2020_onco_ALL_reformatted_coding.txt"

# ================================================================================
# libraries used in the code
# please install the following packages if needed
#BiocManager::install("vroom")
#BiocManager::install("splitstackshape")
library(vroom)
library(splitstackshape)

# ================================================================================
# 						 RUNNING SECTIONS
# ================================================================================
gwas <- vroom(file=paste0(path_input_data, input_data))
gwas <- as.data.frame(gwas)
colnames(gwas) <- c("CHR", "SNP", "BP", "P_ANNOT")

gwas <- cSplit(gwas, "P_ANNOT", " ")
colnames(gwas) <- c("CHR","SNP","BP","P", "ANNOT")

vroom_write(gwas, file=paste0(path_output_data, output_data))

# ================================================================================
# write the details information of the computer and softwares related to this code 
# please keep this info since it might be needed for any scientist in the future!
sessionInfo()
# calculation the running time of the program
end_time <- Sys.time()
end_time - start_time
# ================================================================================
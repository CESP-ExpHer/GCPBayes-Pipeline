# =====================================================================
#   *** Separation of list of groups for GCPBayes analysis according to a selected threshold ***
# =====================================================================
# Summary:
# Separating GCPBayes inputs into subparts corresponding to
# - inputs for Genes with length < threshold (without LD clumping)
# - inputs for Genes with length > threshold (with LD clumping)
# ================================================================================
# Initially Written by: PE Sugier
# Modified by: Yazdan Asgari
# Initial Creation Date: 10/2022
# Edited Date: 03/2023
# ================================================================================

# for calculation the running time of the program
start_time <- Sys.time()

# ================================================================================
# libraries used in the code
library(tictoc)

# ================================================================================
# 						DEFINITION SECTION
# ================================================================================
# PATH in which a GCPBayes input file exists
path_input_data_read <- path_input_data

# ================================================================================
# 						            RUNNING SECTIONS
# ================================================================================
# GCPBayes input file name
# Without clumping
load(paste0(path_input_data_read,"D1_Matrices_output_pipeline_", output_step4wc, ".Rdata"))
load(paste0(path_input_data_read,"D1_Matrices_extra_info_output_pipeline_", output_step4wc, ".Rdata"))
Matrix_in_noclump <- gcpbayes_input_data_final
Matrix_in_noclump_extra_info <- gcpbayes_extra_input_data_final


# PATH for writing output files
path_output_data = path_output_data


# Parameters used in this script
# A threshold for SNP numbers to split Genes based on that 
SNP_number_threshold <- group_clump_threshold
# A threshold for considering a Gene with a potential pleiotropic effect 
theta_threshold <- theta_exploration

out_noclump <- list()
out_noclump_extra <- list()

out_noclump <- Matrix_in_noclump
out_noclump_extra <- Matrix_in_noclump_extra_info

# Save .RData genes
gcpbayes_input_data_final <- out_noclump
gcpbayes_extra_input_data_final <- out_noclump_extra
save(gcpbayes_input_data_final, file=paste0(path_output_data, "D3_Matrices_output_pipeline_", output_step5wc, ".Rdata"))
save(gcpbayes_extra_input_data_final, file=paste0(path_output_data, "D3_Matrices_extra_info_output_pipeline_", output_step5wc,".Rdata"))


# ================================================================================
# write the details information of the computer and softwares related to this code 
# please keep this info since it might be needed for any scientist in the future!
sessionInfo()
# calculation the running time of the program
end_time <- Sys.time()
end_time - start_time
# ================================================================================


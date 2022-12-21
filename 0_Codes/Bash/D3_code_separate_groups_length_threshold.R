#!/usr/bin/env Rscript

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
# Edited Date: 12/2022
# ================================================================================

# for calculation the running time of the program
start_time <- Sys.time()

# ================================================================================
# libraries used in the code
# please install the following packages if needed
library(tictoc)
library(optparse)

# ================================================================================
# Read of command inputs
# ================================================================================

option_list = list(
  make_option(c("--path"), type="character", default=NULL, 
              help="working directory", metavar="character"),
  make_option(c("--file_noclump"), type="character", default=NULL, 
              help="input file name for not clumped genes", metavar="character"),
  make_option(c("--file_clump"), type="character", default=NULL, 
              help="input file name for clumped genes", metavar="character"),
  make_option(c("--t_SNP_number"), type="numeric", default=NULL, 
              help="threshold used for minimum MAF", metavar="number"),
  make_option(c("--out_noclump"), type="character", default="gcpbayes_greater_out.txt", 
              help="output file name for not clumped genes [default= %default]", metavar="character"),
  make_option(c("--out_clump"), type="character", default="gcpbayes_greater_out.txt", 
              help="output file name for clumped genes [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# ================================================================================
# 						DEFINITION SECTION
# ================================================================================
# PATH in which a GCPBayes input file exists
path_input_data_read <- opt$path

# ================================================================================
# 						            RUNNING SECTIONS
# ================================================================================
# GCPBayes input file name
# Without clumping
load(paste0(path_input_data_read,"D1_Matrices_output_pipeline_",opt$file_noclump,".Rdata"))
load(paste0(path_input_data_read,"D1_Matrices_extra_info_output_pipeline_",opt$file_noclump,".Rdata"))
Matrix_in_noclump <- gcpbayes_input_data_final
Matrix_in_noclump_extra_info <- gcpbayes_extra_input_data_final
# With clumping
load(paste0(path_input_data_read,"D2_Matrices_output_pipeline_",opt$file_clump,".Rdata"))
load(paste0(path_input_data_read,"D2_Matrices_extra_info_output_pipeline_",opt$file_clump,".Rdata"))
Matrix_in_clump <- gcpbayes_input_data_final
Matrix_in_clump_extra_info <- gcpbayes_extra_input_data_final

# PATH for writing output files
path_output_data = opt$path


# Parameters used in this script
# A threshold for SNP numbers to split Genes based on that 
SNP_number_threshold <- opt$t_SNP_number
# A threshold for considering a Gene with a potential pleiotropic effect 
theta_threshold <- opt$t_theta


# To get the list of genes with size <= SNP_number_threshold  --> list of short genes
# To get the list of genes with size > SNP_number_threshold  --> list of long genes
# Genes are sorted by size
out_noclump <- list()
out_noclump_extra <- list()
genes_clump <- list()
out_clump <- list()
out_clump_extra <- list()
n <- length(Matrix_in_noclump_extra_info)

if(Matrix_in_noclump_extra_info[[n]]$snp_number<=SNP_number_threshold){
    out_noclump <- Matrix_in_noclump
    out_noclump_extra <- Matrix_in_noclump_extra_info
} else{
    k <- 1
    for (i in 1:length(Matrix_in_noclump_extra_info)){
      if(Matrix_in_noclump_extra_info[[i]]$snp_number<=SNP_number_threshold){
        out_noclump[[i]] <- Matrix_in_noclump[[i]]
        out_noclump_extra[[i]] <- Matrix_in_noclump_extra_info[[i]]
      } else{
        genes_clump[[k]] <- Matrix_in_noclump_extra_info[[i]]$gene
        k <- k + 1
      }
    }
    m <- 1
    for (i in 1:(k-1)) {
      for (j in 1:length(Matrix_in_clump_extra_info)) {
        if (genes_clump[[i]] == Matrix_in_clump_extra_info[[j]]$gene) { # Some genes from "noclump" list could not be retrieved in "clump" list
          out_clump[[m]] <- Matrix_in_clump[[j]]
          out_clump_extra[[m]] <- Matrix_in_clump_extra_info[[j]]
          m <- m + 1
        }
      }
    }
}

#########################################
# Saving outputs

# Save .RData short genes
gcpbayes_input_data_final <- out_noclump
gcpbayes_extra_input_data_final <- out_noclump_extra
save(gcpbayes_input_data_final, file=paste0(path_output_data, "D3_Matrices_output_pipeline_", opt$out_noclump, ".Rdata"))
save(gcpbayes_extra_input_data_final, file=paste0(path_output_data, "D3_Matrices_extra_info_output_pipeline_", opt$out_noclump,".Rdata"))

# Save .RData long genes
gcpbayes_input_data_final <- out_clump
gcpbayes_extra_input_data_final <- out_clump_extra
save(gcpbayes_input_data_final, file=paste0(path_output_data, "D3_Matrices_output_pipeline_", opt$out_clump, ".Rdata"))
save(gcpbayes_extra_input_data_final, file=paste0(path_output_data, "D3_Matrices_extra_info_output_pipeline_", opt$out_clump,".Rdata"))


# ================================================================================
# write the details information of the computer and softwares related to this code 
# please keep this info since it might be needed for any scientist in the future!
sessionInfo()
# calculation the running time of the program
end_time <- Sys.time()
end_time - start_time
# ================================================================================


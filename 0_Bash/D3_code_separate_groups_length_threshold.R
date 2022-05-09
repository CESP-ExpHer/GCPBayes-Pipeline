#!/usr/bin/env Rscript
library("optparse")

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


#   ***Running GCPBayes for a pair of traits***
# =====================================================================
# Summary:
# To get correct list of groups for further analysis according to wanted threshold (length of groups)
# Separating GCPBayes inputs into subparts corresponding to
# - inputs for Genes with length < threshold (without LD clumping)
# - inputs for Genes with length > threshold (with LD clumping)
# ================================================================================
# Initially Written by: PE Sugier
# Modified by 
# Edited Date: 01/2022
# ================================================================================
# for calculation the running time of the program
start_time <- Sys.time()

# ================================================================================
# 						DEFINITION SECTION
# should be changed by a user
# ================================================================================
# PATH in which a GCPBayes input file exists
path_input_data_read <- opt$path

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

# ================================================================================
# libraries used in the code
# please install the following packages if needed
# BiocManager::install("arm")
# library(devtools)
# devtools::install_github("https://github.com/cran/bglm")
# devtools::install_github("https://github.com/nyiuab/BhGLM.git")
# devtools::install_github("https://github.com/tbaghfalaki/GCPBayes")
library(tictoc)
library(GCPBayes)
library(BhGLM)

# ================================================================================
# 						RUNNING SECTION
# ================================================================================

# To get the list of genes with size <= SNP_number_threshold  --> list of short genes
# To get the list of genes with size > SNP_number_threshold  --> list of long genes
# Genes are sorted by size
out_noclump <- list()
out_noclump_extra <- list()
for (i in 1:length(Matrix_in_noclump_extra_info)){
	if(Matrix_in_noclump_extra_info[[i]]$snp_number<=SNP_number_threshold){
		out_noclump[[i]] <- Matrix_in_noclump[[i]]
		out_noclump_extra[[i]] <- Matrix_in_noclump_extra_info[[i]]
	}
}

# To get the name of the last gene to drop
n_last_short <- length(out_noclump_extra)
name_last_short <- out_noclump_extra[[n_last_short]]$gene

# To get the indice of the gene in the clumped data
j=1
while (Matrix_in_clump_extra_info[[j]]$gene!=name_last_short){
	j=j+1
}
# To drop all short genes from the clumped data (already considered in the short data) of the list of the long genes
out_clump <- Matrix_in_clump[-seq(1,j)]
out_clump_extra <- Matrix_in_clump_extra_info[-seq(1,j)]


#########################################
# Replacing sigma matrix if required





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


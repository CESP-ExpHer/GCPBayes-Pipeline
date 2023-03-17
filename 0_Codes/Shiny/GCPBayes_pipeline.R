# ================================================================================
# 					      *** GCPBAYES PIPELINE ***
# ================================================================================
# Summary:
# ================================================================================
# Written first by: Yazdan Asgari
# Modified by: Yazdan Asgari
# Initial Creation Date: 02/2023
# Edited Date: 03/2023
# ================================================================================
# function for getting the path of current script
library(tidyverse)
getCurrentFileLocation <-  function()
{
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}

# ==============================================================
setwd(getCurrentFileLocation())

source("GCPBayes_pipeline_parameters.R")


# Run script C1 (Finding Common SNPs)
source("C1_code_find_common_snps_one_pair.R")

# Run script C2 (Running PLACO)
if (toclump == TRUE){
  # Checking existence of PLACO output file
  check_PLACO_file <- paste0(path_input_data, output_step2, ".txt")
  if (file.exists(check_PLACO_file)) {
    message("The file ", output_step2, " which is the output of PLACO step, exists. Therefore, this step does not need to be run.")
  }else {
    source("C2_code_run_PLACO_decor_one_pair.R")
  }
}

# Run script C3 (Running LD Clumping)
if (toclump == TRUE){
  # Checking existence of LD Clumping output file
  check_LD_file <- paste0(path_input_data, output_step3, ".txt")
  if (!file.exists(check_LD_file)) {
    # Checking existence of Path for GWAS reference file
    if (!file.exists(ref_path_b_files)) {
      message("The path ", ref_path_b_files, " does not exist. During running LD Clumping, GCPBayes Pipeline looks for bfiles in this Path.")
      stop("Path for GWAS reference files not found !")
    }else {
      source("C3_code_ldclumping_local.R")
    }
  }  else {
    message("The file ", output_step3, " which is the output of LD Clumping step, exists. Therefore, this step does not need to be run.")
  }
}

# Run script D1 (Main Pipeline Steps)
source("D1_code_pipeline_annot_coding_withoutldclumping_extra_info.R")

# Run script D2 (Main Pipeline Steps if LD Clumping done before)
if (toclump == TRUE){
  # Checking existence of LD file
  check_LD_file <- paste0(path_input_data, output_step3, ".txt")
  if (!file.exists(check_LD_file)) {
    message("The file ", check_LD_file, " which is required does not exist. Please first do the LD Clumping Step.")
    stop("File not found !")
  } else {
    source("D2_code_pipeline_annot_coding_ldclumping_extra_info.R")
  }
}

# Run script D3 (Preparation of GCPBayes Input Files)
if (toclump == TRUE){
  # Run script D3
  source("D3_code_separate_groups_length_threshold.R")
} else {
  # Run script D3_no_clump
  source("D3_code_separate_groups_length_threshold_noclump.R")
}

# Run script E1 (Running GCPBayes)
source("E1_code_gcpbayes_less_extra_info.R")

# Run script E1 (Running GCPBayes if LD Clumping done before)
if (toclump == TRUE){
  source("E1_code_gcpbayes_clump_less_extra_info.R")
}



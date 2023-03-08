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

source("GCPBayes_pipeline_parameters_test.R")

# Run script C1 (Finding Common SNPs)
source("R_C1_code_find_common_snps_one_pair_test.R")

# Run script D1 (Main Pipeline Steps)
source("R_D1_code_pipeline_annot_coding_withoutldclumping_extra_info_test.R")

# Run script D3_no_clump
source("R_D3_code_separate_groups_length_threshold_noclump_test.R")

# Run script E1 (Running GCPBayes)
source("R_E1_code_gcpbayes_less_extra_info_test.R")




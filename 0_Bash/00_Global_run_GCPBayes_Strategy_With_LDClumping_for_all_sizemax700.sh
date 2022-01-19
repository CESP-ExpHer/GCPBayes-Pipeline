#===============================================================================
#           USAGE:  ./00_Global_run_GCPBayes_Strategy_With_LDClumping_for_all_sizemax700.sh
#
#     DESCRIPTION:  global codes to run the GCPBayes pipeline (encapsulate all the different procedures in one code, with a global specification of the parameters for the user)
#
#===============================================================================
#          AUTHOR:  Pierre-Emmanuel Sugier, Yazdan Asgari
#         VERSION:  Pierre-Emmanuel Sugier
#         CREATED:  2021-11-23
#        REVISION:  2022-01-18
# WARNINGS N BUGS:  
#                   
#    NOTES N TODO:  
#                   Run HS
#                   Manhattan Plots (global figure)
#                   Figure per group with theta > theta_exploration
#                   Figure for PLACO?
#                   Parameter chr_start chr_end : add in scripts + maybe it is better to consider a list (more flexible)
#                   All strategy based on covariance matrix estimated from a genetic LD matrix from reference, to build IN the codes step 4 (gcpaybes_corr_ref to implement)
#                   
# ================================================================================
#      PARAMETERS:  
#               
#                  ---to specify by user
#     
#                       work_dir                             - working directory (to put files of each step)
#                       output_dir                           - output directory (to put all final outputs)
#                       
#                       theta_exploration=0.1              - threshold for theta values obtained from DS to run the group with HS
#                       theta_plot=0.5                     - theta value threshold (with DS or HS) to plot the group
#                       gcpaybes_corr_ref="FALSE"          - to perform GCPBayes by using a reference matrix (values: "TRUE" or "FALSE")
#
#                       group_clump_threshold=500          - length threshold used to consider groups (more longer genes are only performed with previous LD clumping if toclump="TRUE")
#                       group_absolute_threshold=1200      - length threshold used to consider groups (analysis will not be performed for more longer genes / If toclump=TRUE", this threshold is evaluated after clumping)
#                       placo_pval_threshold=0.05          - threshold used for decorrelating the Z-scores (in PLACO)
#                       toclump="TRUE"                     - If perform LD Clumping or not (value: TRUE or FALSE)
#                       clump_threshold_r2   =0.8          - LD threshold used for the LD clumping based on r²
#                       clump_threshold_kb=10000           - Maximum distance in based pair used to consider clumping
#                       clump_threshold_p=0.99             - Threshold based on p for ld clumping
#                       
#                       short_data1                          - short name to call dataset 1
#                       short_data2                          - short name to call dataset 2
#
#                       Input_file_trait1                    - input file name for dataset 1
#                       path_inputfile_trait1                - path for dataset 1
#                       Input_file_trait2                    - input file name for dataset 2
#                       path_inputfile_trait2                - path for dataset 2
#                       
#                       chr_start=1                        - first chromosome considered for analysis
#                       chr_end=22                         - last chromosome considered for analysis
#                       
#                       info_threshold=0.9                 - threshold used for minimum imputation quality of SNPs (GWAS QC)
#                       MAF_threshold=0.05                 - threshold used for minimum MAF (GWAS QC)
#
#                       
#                   ---Hard coded inside scripts that are called by the pipeline ---    
#
#                       -> Maximum distance in based pair used to consider clumping (used for example kb=10000)
#                       -> Threshold based on p (used for example p=0.99)
#                         
#                       
# ================================================================================					
#           INPUTS:  
#                   - GWAS data reformatted for each dataset			
#					   An example of the first lines of an input file (reformatted) is:
#                   - reference bfiles for ldclumping
#                   - reference file for annotations
#                   - 
# ================================================================================
#          DIFFERENT STEPS:
#                    
#                    Step 1 : To extract common SNPs between dataset
#                    Step 2 : To run PLACO on each dataset
#                    Step 3 : To perform LD Clumping (locally)
#                    Step 4 : To prepare the GCPBayes inputs in the right format
#                    Step 5 : To get correct list of groups for further analysis according to wanted threshold (length of groups)
#                    Step 6 : To run GCPBayes (DS)
#                    Step 7 : To run GCPBayes (HS) on groups with theta > theta_exploration
#                    Step 8 : To plot figures
#                    
# ================================================================================                   
#          OUTPUTS:  
#                   Results tables for DS
#                   Results tables for HS (to add)
#                   Results for PLACO for each phenotype
#                   Manhattan Plot (to add)
#                   Figures per genes (to add)
#                   Figures for PLACO ?
# ================================================================================                 
#    REQUIREMENTS:  R
# 	 REQUIRED LIBRARIES
# libraries used in the code
# please install the following packages if needed
#     BiocManager::install("devtools")
#     install.packages('tictoc')
#     BiocManager::install("arm")
#     library(devtools)
#     devtools::install_github("https://github.com/cran/bglm")
#     devtools::install_github("https://github.com/nyiuab/BhGLM.git")
#     devtools::install_github("https://github.com/tbaghfalaki/GCPBayes")
#     BiocManager::install("vroom")
#     BiocManager::install("dplyr")
#     BiocManager::install("data.table")
#     BiocManager::install("tidyr")
#     BiocManager::install("tidyverse")
#     install.packages("optparse")
# ================================================================================


##########################################################################################################################################
##       USER SPECIFICATIONS       ##
#####################################

# working directory
work_dir="/home1/6_AMLAP/GCPBayes/Pipeline/00.pipeline_global/work/"
# output directory
output_dir="/home1/6_AMLAP/GCPBayes/Pipeline/00.pipeline_global/work/"

# short names
short_data1="BCAC"
short_data2="OCAC"

# Input files - datasets
Input_file_trait1="BCAC_2020_onco_ALL_reformatted.txt"
path_inputfile_trait1="/home1/6_AMLAP/GCPBayes/Pipeline/00.pipeline_global/inputs/"
Input_file_trait2="OCAC_BCAC_2020_onco_ALL_reformatted.txt"
path_inputfile_trait2="/home1/6_AMLAP/GCPBayes/Pipeline/00.pipeline_global/inputs/"

# chr considered for analysis
#chr_start=1
#chr_end=22

# Input files - annotation
path_gwas_annot=${path_inputfile_trait1} #directory in which the annotated data for the first GWAS exist
file_gwas_annot="Annot_BCAC_2020_onco_ALL_reformatted_coding" #the file name of the annotated data for the first GWAS
path_annot="/home1/6_AMLAP/GCPBayes/Pipeline/00_data_ref/"
file_annot="annot_gencode_v38lift37_modified_gene_class.txt"

# Input files - reference data
ref_path_b_files="/home1/6_AMLAP/GCPBayes/Pipeline/00_data_ref/1000K_ref_data/EUR/EUR" # directory path for the reference bfiles

# QC parameters
info_threshold=0.9
MAF_threshold=0.05

# Parameters
theta_exploration=0.5
theta_plot=0.5
group_clump_threshold=0 # More longer groups will be "ld-clumped" before analysis (if clumping is FALSE, then they will not be performed)
group_absolute_threshold=700 # More longer groups (after LD clumping) will not be performed. We recommand to not put this threshold upper than 1500
placo_pval_threshold=0.05
toclump="TRUE" # TRUE or FALSE
clump_threshold_r2=0.8 # LD threshold used for the LD clumping based on r²
clump_threshold_kb=10000 # Maximum distance in based pair used to consider clumping
clump_threshold_p=0.99 # Threshold based on p for ld clumping



##########################################################################################################################################
##       HARD CODED      ##
###########################

output_step1="step1_output"
output_step2="step2_PLACO_output"
output_step3="step3_ldclumping_output"
output_step4c=${short_data1}_${short_data2}_coding_clumping_${clump_threshold_r2}
output_step4wc=${short_data1}_${short_data2}_coding_withoutclumping
output_step5c=${short_data1}_${short_data2}_genes_longer_than_${group_clump_threshold}
output_step5wc=${short_data1}_${short_data2}_genes_shorter_than_${group_clump_threshold}
output_step6wc=output_GCPBayes_${short_data1}_${short_data2}_without_ldclumping
output_step6c=output_GCPBayes_${short_data1}_${short_data2}_with_ldclumping_${clump_threshold_r2}

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##       START RUN       ##
###########################

# STEP 1 - To extract common SNPs between dataset
#####################################################

Rscript C1_code_find_common_snps_one_pair.R \
											--path1 ${path_inputfile_trait1} \
											--file1 ${Input_file_trait1}.txt \
											--path2 ${path_inputfile_trait2} \
											--file2 ${Input_file_trait2}.txt \
											--pathout ${work_dir} \
											--out1 ${output_step1}_${short_data1} \
											--out2 ${output_step1}_${short_data2}_common_${short_data1} \
											--info ${info_threshold} \
											--maf ${MAF_threshold}

# Step 2 - To run PLACO on each dataset
#####################################################

if [Clumping="TRUE"]; then
  if [ -e ${work_dir}/${output_step2}_${short_data1}_${short_data2}.txt ] # this step does not run if a corresponding output file already exist
  then
    echo "A corresponding PLACO output is present in the work directory and is going to be used."
  else
    Rscript C2_code_run_PLACO_decor_one_pair.R \
        --path ${work_dir} \
        --file1 ${output_step1}_${short_data1} \
        --file2 ${output_step1}_${short_data2}_common_${short_data1} \
        --pathout ${work_dir} \
        --out ${output_step2}_${short_data1}_${short_data2} \
        --pval ${placo_pval_threshold}
  fi
fi

# Step 3 - To perform LD Clumping (locally)
#####################################################

if [Clumping="TRUE"]; then
  Rscript C3_code_ldclumping_local.R \
        --path ${work_dir} \
        --file ${output_step2}_${short_data1}_${short_data2} \
        --pathout ${work_dir} \
        --out ${output_step3}_${clump_threshold_r2}_${short_data1}_${short_data2} \
        --ref ${ref_path_b_files} \
        --r2 ${clump_threshold_r2} \
        --kb ${clump_threshold_kb} \
        --p ${clump_threshold_p}
fi

# Step 4 - To prepare the GCPBayes inputs in the right format
#####################################################

Rscript D1_code_pipeline_annot_coding_withoutldclumping_extra_info.R \
        --path ${work_dir} \
        --path1 ${path_inputfile_trait1} \
        --path2 ${path_inputfile_trait2} \
        --file1 ${Input_file_trait1} \
        --file2 ${Input_file_trait2} \
        --name1 ${short_data1} \
        --name2 ${short_data2} \
        --annot_gwas_path ${path_gwas_annot} \
        --annot_gwas_file ${file_gwas_annot} \
        --annot_path ${path_annot} \
        --annot_file ${file_annot} \
        --info ${info_threshold} \
        --maf ${MAF_threshold} \
        --out ${output_step4wc}

if [Clumping="TRUE"]; then
  Rscript D2_code_pipeline_annot_coding_ldclumping_extra_info.R \
        --path ${work_dir} \
        --path1 ${path_inputfile_trait1} \
        --path2 ${path_inputfile_trait2} \
        --file1 ${Input_file_trait1} \
        --file2 ${Input_file_trait2} \
        --name1 ${short_data1} \
        --name2 ${short_data2} \
        --annot_gwas_path ${path_gwas_annot} \
        --annot_gwas_file ${file_gwas_annot} \
        --annot_path ${path_annot} \
        --annot_file ${file_annot} \
        --ldclumpout ${output_step3}_${clump_threshold_r2}_${short_data1}_${short_data2} \
        --info ${info_threshold} \
        --maf ${MAF_threshold} \
        --out ${output_step4c}
fi

# Step 5 - To get correct list of groups for further analysis according to wanted threshold (length of groups)
#####################################################

Rscript D3_code_separate_groups_length_threshold.R \
        --path ${work_dir} \
        --file_noclump ${output_step4wc} \
        --file_clump ${output_step4c} \
        --out_noclump ${output_step5wc} \
        --out_clump ${output_step5c} \
        --t_SNP_number ${group_clump_threshold}
        

# Step 6 - To run GCPBayes (DS)
#####################################################

# A - To run GCPBayes (DS) on groups without ld clumping (length <= group_clump_threshold)
Rscript E1_code_gcpbayes_less_extra_info.R \
        --path ${work_dir} \
        --file ${output_step5wc} \
        --t_SNP_number ${group_clump_threshold} \
        --t_theta ${theta_exploration} \
        --pathout ${output_dir} \
        --out ${output_step6wc}
# WARNINGS : theta_exploration is used (not theta_plot)

# B - To run GCPBayes (DS) on groups after ld clumping (group_clump_threshold < length <= group_absolute_threshold)
if [Clumping="TRUE"]; then
Rscript E1_code_gcpbayes_less_extra_info.R \
        --path ${work_dir} \
        --file ${output_step5c} \
        --t_SNP_number ${group_absolute_threshold} \
        --t_theta ${theta_exploration} \
        --pathout ${output_dir} \
        --out ${output_step6c}
# WARNINGS : theta_exploration is used (not theta_plot)
fi

# Step 7 - To run GCPBayes (HS) on groups with theta > theta_exploration
#####################################################

# Step 8 - To plot figures
#####################################################





































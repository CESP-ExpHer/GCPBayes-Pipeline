# ================================================================================
# 			This is a file for definition of Parameters using for GCPBayes Pipeline
# --------------------------------------------------------------------------------
#     Please read comments or NOTES available for each parameter before changing
# --------------------------------------------------------------------------------
#                         AUTHOR:  Yazdan Asgari
#                         REVISION:  2023-03-01
# ================================================================================

# ===================================
##       PATH SPECIFICATIONS       ##
# ===================================
# working directory
work_dir <- "C:/CESP/AMLAP/4_Manuscripts/GCPBayes-Pipeline/0_code_r_all_in_one/test_dataset/"

# ===================================
##    INPUT GWAS FILES NAMES       ##
# ===================================
# files names for GWAS input files
# NOTE: All names MUST be in quotation and devided by comma
input <- c("gwas_BCAC_chr5.txt", "gwas_OCAC_chr5.txt")
# A short name for each GWAS input file
# NOTE: All names MUST be in quotation and devided by comma
input_shortname <- c("BCAC", "OCAC")

# ===================================
##    GWAS HEADERS NAMES         ##
# ===================================
# GWAS #1
g1_rsid <- "snp"           #(SNP rs id)
g1_chr  <- "chr"   #(chromosome number)
g1_pos  <- "bp_hg19"          #(base pair position)
g1_EA   <- "Effect_A"           #(Effect Allele)
g1_nE_A <- "Baseline_A"           #(non-Effect Allele)
g1_beta <- "beta"         #(beta value)
g1_se   <- "se"           #(standard error)
g1_pval <- "pval"         #(P-value)
g1_info <- "info"         #(imputation quality value)
g1_EAF  <- "EAF"          #(Effect Allele Frequency)
g1_MAF  <- "MAF"          #(Minor Allele Frequency)

# GWAS #2
g2_rsid <- "snp"          #(SNP rs id)
g2_chr  <- "chr"  #(chromosome number)
g2_pos  <- "bp_hg19"         #(base pair position)
g2_EA   <- "Effect_A"          #(Effect Allele)
g2_nE_A <- "Baseline_A"          #(non-Effect Allele)
g2_beta <- "beta"            #(beta value)
g2_se   <- "se"           #(standard error)
g2_pval <- "pval"            #(P-value)
g2_info <- "info"         #(imputation quality value)
g2_EAF  <- "EAF"          #(Effect Allele Frequency)
g2_MAF  <- "MAF"          #(Minor Allele Frequency)

# ===================================
##       ANNOTATION FILES NAMES    ##
# ===================================
# GENCODE Annotation file name
file_annot <- "annot_gencode_v38lift37_modified_gene_class.txt"
# the file name of the annotated data for the first GWAS
file_gwas_annot <- "Annot_BCAC_2020_onco_ALL_reformatted_coding.txt"

# ===================================
##          QC parameters          ##
# ===================================
info_threshold <- 0.8
MAF_threshold <- 0.01

# ===================================
##   Number of SNPs Threshold      ##
# ===================================
group_clump_threshold <- 700

# ===================================
##    Pleiotropy Threshold         ##
# ===================================
theta_exploration <- 0.5

# =================================================
##   Decision for running LD Clumping Step       ##
# =================================================
# If LD_Clumping should be run or Not (TRUE or FALSE)
toclump <- FALSE

# =====================================================
# NOTE: If toclump==TRUE, these parameters will be used
# =====================================================
# directory path for the reference GWAS bfiles
ref_path_b_files <- "/PATH/"
# GWAS Reference Files Name (
# NOTE: All three "bed/bim/fam" files MUST have the same name
ref_b_files <- "EUR"
# LD threshold used for the LD clumping based on r2
clump_threshold_r2 <- 0.8
# Maximum distance in based pair used to consider clumping
clump_threshold_kb <- 10000
# Threshold based on p for ld clumping
clump_threshold_p <- 0.99
# P_value used for PLACO analysis
placo_pval_threshold <- 0.05
# ================================================================================
# ================================================================================
# ================================================================================
# ================================================================================

#####################################
##       HARD CODED SECTION        ##
##    NOTE: No Need to Change      ##
#####################################
number_of_studies <- 2
output_dir <- work_dir
script_dir <- work_dir
path_input_data  <- work_dir
path_output_data <- output_dir
path_gwas_annot  <- script_dir
path_annot     <- script_dir
output_step1   <- "step1_output"
output_step2   <- "step2_PLACO_output"
output_step3   <- "step3_ldclumping_output"
output_step4c  <- paste0("output_clump_", clump_threshold_r2)
output_step4wc <- "output_without_clump"
output_step5c  <- paste0("output_genes_longer_than_", group_clump_threshold)
output_step5wc <- paste0("output_genes_shorter_than_", group_clump_threshold)
output_step6c  <- paste0("output_GCPBayes_clump_", clump_threshold_r2)
output_step6wc <- "output_GCPBayes_without_clump"

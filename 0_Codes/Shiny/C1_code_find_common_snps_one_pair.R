# ================================================================================
# 					EXTRACTION OF COMMON SNPs + Z calculation
#			 			OUTPUT will be used for running PLACO 
# ================================================================================
# Summary: The code reads a pair of GWAS data and extract common SNPs and stores
# the rows in the separated output files
# ================================================================================
# Written first by: Yazdan Asgari & Pierre-Emmanuel Sugier
# Modified by: Yazdan Asgari
# Initial Creation Date: 12/2020
# Edited Date: 03/2023
# ================================================================================

# for calculation the running time of the program
start_time <- Sys.time()

# ================================================================================
# libraries used in the code
library(vroom)
library(tidyverse)
library(data.table)

# ================================================================================
# 						DEFINITION SECTION
# ================================================================================
nstud <- number_of_studies
output_names <- output_step1

# ================================================================================
# 						            RUNNING SECTIONS
# ================================================================================
# Summary: 
# Reading GWAS data, filtering some SNPs, calculation of Z based on "beta" and "se"
# Extraction of common SNPs between pair of GWAS data
# ================================================================================

# treat_1: very fast reading method
gwas_traits_1 <- vroom(file = paste0(path_input_data, input[1]))
gwas_traits_1 <- as.data.frame(gwas_traits_1)

# Definition of old and new header names
old_headers_1 <- c(g1_rsid, g1_chr, g1_pos, g1_EA, g1_nE_A, g1_beta, g1_se, g1_pval, g1_info, g1_EAF, g1_MAF)
new_headers_1 <- c("snp", "chr", "bp_hg19", "Effect_A", "Baseline_A", "beta", "se", "pval", "info", "EAF", "MAF")

# Finding indices of old headers in the treat_1
header_indices_1 <- match(old_headers_1, names(gwas_traits_1))

# Renaming the columns
names(gwas_traits_1)[header_indices_1] <- new_headers_1


# filtering step
gwas_traits_1 <- gwas_traits_1 %>% filter(gwas_traits_1$info > info_threshold & gwas_traits_1$MAF > MAF_threshold) %>% dplyr::select(c("snp","chr","bp_hg19","Effect_A","Baseline_A","beta","se","pval","info","EAF","MAF"))

# calculation of Z
# if se==0, put a VERY LARGE VALUE for Z (here 1.0e+08)
gwas_traits_1 %>% mutate(Z = case_when(gwas_traits_1$se == 0   ~ 1.0e+08 , gwas_traits_1$se != 0 ~ (gwas_traits_1$beta/gwas_traits_1$se) )) -> gwas_traits_1


################################  
# to read other studies to get pool of SNPs to keep

for (i in 2:nstud){
  # treat_2: very fast reading method
  gwas_traits_2 <- vroom(file = paste0(path_input_data,input[i]))
  gwas_traits_2 <- as.data.frame(gwas_traits_2)
  
  # Definition of old and new header names
  old_headers_2 <- c(g2_rsid, g2_chr, g2_pos, g2_EA, g2_nE_A, g2_beta, g2_se, g2_pval, g2_info, g2_EAF, g2_MAF)
  new_headers_2 <- c("snp", "chr", "bp_hg19", "Effect_A", "Baseline_A", "beta", "se", "pval", "info", "EAF", "MAF")
  
  # Finding indices of old headers in the treat_1
  header_indices_2 <- match(old_headers_2, names(gwas_traits_2))
  
  # Renaming the columns
  names(gwas_traits_2)[header_indices_2] <- new_headers_2
  
  # filtering step
  gwas_traits_2 <- gwas_traits_2 %>% filter(gwas_traits_2$info > info_threshold & gwas_traits_2$MAF > MAF_threshold) %>% dplyr::select(c("snp","chr","bp_hg19","Effect_A","Baseline_A","beta","se","pval","info","EAF","MAF"))
  
  # calculation of Z
  # if se==0, put a VERY LARGE VALUE for Z (here 1.0e+08)
  gwas_traits_2 %>% mutate(Z = case_when(gwas_traits_2$se == 0   ~ 1.0e+08 , gwas_traits_2$se != 0 ~ (gwas_traits_2$beta/gwas_traits_2$se) )) -> gwas_traits_2
  
  # merging two GWAS data based on "snp", "chr", and "bp_hg19" columns
  gwas_merge <- inner_join(gwas_traits_1,gwas_traits_2, by = c("snp", "chr", "bp_hg19"))
  
  # checking number of total, unique, and duplicated rsids
  length(gwas_merge$snp)
  length(unique(gwas_merge$snp))
  dim(gwas_merge[duplicated(gwas_merge$snp), ])[1]
  
  # removing potential duplicated entries
  gwas_merge <- gwas_merge[!duplicated(gwas_merge$snp), ]
  
  # checking data after removing the duplicated entries
  length(gwas_merge$snp)
  length(unique(gwas_merge$snp))
  dim(gwas_merge[duplicated(gwas_merge$snp), ])[1]
  
  # extracting common SNPs for the first trait
  gwas_traits_1 <- gwas_merge[, c("snp", "chr", "bp_hg19", 
                                  "Effect_A.x", "Baseline_A.x", "beta.x", "se.x", 
                                  "pval.x", "info.x", "EAF.x", "MAF.x", "Z.x") ]
  colnames(gwas_traits_1) <- c("snp","chr","bp_hg19","Effect_A","Baseline_A","beta","se","pval","info","EAF","MAF","Z")
}                              

writeLines("\n\n")
print(input_shortname[1])
writeLines("\n\n")

# sorting based on the ID column (almost fast: about 30 sec)
output_traits_1 <- as.data.table(gwas_traits_1)
output_traits_1 <- setorder(output_traits_1, snp)
# writing the output file 1 (very fast)
vroom_write(output_traits_1, file = paste0(path_output_data, output_names, "_", input_shortname[1], ".txt"))


################################                     
# To read other studies to get corresponding outputs

for (i in 2:nstud){
  
  print(input_shortname[i])
  writeLines("\n\n")
  # treat_2: very fast reading method
  gwas_traits_2 <- vroom(file = paste0(path_input_data, input[i]))
  gwas_traits_2 <- as.data.frame(gwas_traits_2)
  
  # Definition of old and new header names
  old_headers_2 <- c(g2_rsid, g2_chr, g2_pos, g2_EA, g2_nE_A, g2_beta, g2_se, g2_pval, g2_info, g2_EAF, g2_MAF)
  new_headers_2 <- c("snp", "chr", "bp_hg19", "Effect_A", "Baseline_A", "beta", "se", "pval", "info", "EAF", "MAF")
  
  # Finding indices of old headers in the treat_1
  header_indices_2 <- match(old_headers_2, names(gwas_traits_2))
  
  # Renaming the columns
  names(gwas_traits_2)[header_indices_2] <- new_headers_2
  
  # filtering step
  gwas_traits_2 <- dplyr::select(filter(gwas_traits_2, gwas_traits_2$info > info_threshold & gwas_traits_2$MAF > MAF_threshold),
                                 c("snp","chr","bp_hg19","Effect_A","Baseline_A","beta","se","pval","info","EAF","MAF"))
  
  # calculation of Z
  # if se==0, put a VERY LARGE VALUE for Z (here 1.0e+08)
  gwas_traits_2 %>% mutate(Z = case_when(gwas_traits_2$se == 0   ~ 1.0e+08 , gwas_traits_2$se != 0 ~ (gwas_traits_2$beta/gwas_traits_2$se) )) -> gwas_traits_2
  
  # merging two GWAS data based on "snp", "chr", and "bp_hg19" columns
  gwas_merge <- inner_join(gwas_traits_1,gwas_traits_2, by = c("snp", "chr", "bp_hg19"))
  
  # checking number of total, unique, and duplicated rsids
  length(gwas_merge$snp)
  length(unique(gwas_merge$snp))
  dim(gwas_merge[duplicated(gwas_merge$snp), ])[1]
  
  # removing potential duplicated entries
  gwas_merge <- gwas_merge[!duplicated(gwas_merge$snp), ]
  
  # checking data after removing the duplicated entries
  length(gwas_merge$snp)
  length(unique(gwas_merge$snp))
  dim(gwas_merge[duplicated(gwas_merge$snp), ])[1]
  
  # extracting common SNPs for the second trait
  gwas_traits_2 <- gwas_merge[, c("snp", "chr", "bp_hg19", 
                                  "Effect_A.y", "Baseline_A.y", "beta.y", "se.y", 
                                  "pval.y", "info.y", "EAF.y", "MAF.y", "Z.y") ]
  
  colnames(gwas_traits_2) <- c("snp", "chr", "bp_hg19", 
                               "Effect_A", "Baseline_A", "beta", "se", 
                               "pval", "info", "EAF", "MAF", "Z") 
  
  # sorting based on the ID column (almost fast: about 30 sec)
  output_traits_2 <- as.data.table(gwas_traits_2)
  output_traits_2 <- setorder(output_traits_2, snp)
  
  # writing the output files (very fast)
  vroom_write(output_traits_2, file = paste0(path_output_data, output_names, "_", input_shortname[i], ".txt")) 
}

writeLines("\n\n")
print(paste0('End of C1 Running'))
writeLines("\n\n")

# ================================================================================
# write the details information of the computer and softwares related to this code 
# please keep this info since it might be needed for any scientist in the future!
sessionInfo()
# calculation the running time of the program
end_time <- Sys.time()
end_time - start_time
# ================================================================================


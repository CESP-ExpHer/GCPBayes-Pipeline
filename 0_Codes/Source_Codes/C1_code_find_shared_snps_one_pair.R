# 					EXTRACTION OF SHARED SNPs + Z calculation
#			 			OUTPUT will be used for running PLACO 
# ================================================================================
# Summary: The code reads a pair of GWAS data and extract shared SNPs and stores
# the rows in the separated output files
# ================================================================================
# Written first by: Elise Lucotte
# Modified by: Yazdan Asgari
# Initial Creation Date: 12/2020
# Edited Date: 1/2022
# https://cesp.inserm.fr/en/equipe/exposome-and-heredity
# ================================================================================
# for calculation the running time of the program
start_time <- Sys.time()

# ================================================================================
# 						DEFINITION SECTION
# should be changed by a user
# ================================================================================

# directory in which input data exist
path_input_data_trait1 = "~/BCAC_OCAC/"
path_input_data_trait2 = "~/BCAC_OCAC/"

# the file names of input data
traits_1 <- "BCAC_2020_onco_ALL_reformatted.txt"

traits_2 <- "OCAC_BCAC_2020_onco_ALL_reformatted.txt"

# directory in which an output data would be written
path_output_data = "~/BCAC_OCAC/"

# the file names of the output data 
output_names <- c(  "BCAC_2020_ALL_Shared_OCAC_inc_Z",
                    "OCAC_Shared_BCAC_2020_ALL_inc_Z")

# used for filtering GWAS data              
info_threshold <- 0.8
MAF_threshold <- 0.01

# ================================================================================
# libraries used in the code
# please install the following packages if needed
#BiocManager::install("vroom")
#BiocManager::install("dplyr")
#BiocManager::install("data.table")
library(vroom)
library(dplyr)
library(data.table)

# ================================================================================
# 						            RUNNING SECTIONS
# ================================================================================
# Summary: 
# Reading GWAS data, filtering some SNPs, calculation of Z based on "beta" and "se"
# Extraction of shared SNPs between pair of GWAS data
# ================================================================================

writeLines("\n\n")
print(traits_1)
writeLines("\n\n")
print(traits_2)
writeLines("\n\n")

# treat_1: very fast reading method
gwas_traits_1 <- vroom(file=paste0(path_input_data_trait1, traits_1))
gwas_traits_1 <- as.data.frame(gwas_traits_1)

# filtering step
gwas_traits_1 <- select(filter(gwas_traits_1, gwas_traits_1$info > info_threshold & gwas_traits_1$MAF > MAF_threshold),
                        c("snp","chr","bp_hg19","Effect_A","nonEffect_A","beta","se","pval","info","EAF","MAF"))

# calculation of Z
# if se==0, put a VERY LARGE VALUE for Z (here 1.0e+08)
gwas_traits_1 %>% mutate(Z = case_when(gwas_traits_1$se == 0   ~ 1.0e+08 , gwas_traits_1$se != 0 ~ (gwas_traits_1$beta/gwas_traits_1$se) )) -> gwas_traits_1

# treat_2: very fast reading method
gwas_traits_2 <- vroom(file=paste0(path_input_data_trait2, traits_2))
gwas_traits_2 <- as.data.frame(gwas_traits_2)

# filtering step
gwas_traits_2 <- select(filter(gwas_traits_2, gwas_traits_2$info > info_threshold & gwas_traits_2$MAF > MAF_threshold),
                        c("snp","chr","bp_hg19","Effect_A","nonEffect_A","beta","se","pval","info","EAF","MAF"))

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

# extracting shared SNPs for the first trait
gwas_traits_1 <- gwas_merge[, c("snp", "chr", "bp_hg19", 
                                "Effect_A.x", "nonEffect_A.x", "beta.x", "se.x", 
                                "pval.x", "info.x", "EAF.x", "MAF.x", "Z.x") ]

colnames(gwas_traits_1) <- c("snp", "chr", "bp_hg19", 
                             "Effect_A", "nonEffect_A", "beta", "se", 
                             "pval", "info", "EAF", "MAF", "Z") 

# sorting based on the ID column (almost fast: about 30 sec)
gwas_traits_1 <- as.data.table(gwas_traits_1)
gwas_traits_1 <- setorder(gwas_traits_1, snp)

# extracting shared SNPs for the second  trait
gwas_traits_2 <- gwas_merge[, c("snp", "chr", "bp_hg19", 
                                "Effect_A.y", "nonEffect_A.y", "beta.y", "se.y", 
                                "pval.y", "info.y", "EAF.y", "MAF.y", "Z.y") ]

colnames(gwas_traits_2) <- c("snp", "chr", "bp_hg19", 
                             "Effect_A", "nonEffect_A", "beta", "se", 
                             "pval", "info", "EAF", "MAF", "Z") 

# sorting based on the ID column (almost fast: about 30 sec)
gwas_traits_2 <- as.data.table(gwas_traits_2)
gwas_traits_2 <- setorder(gwas_traits_2, snp)

# writing the output files (very fast)
vroom_write(gwas_traits_1, file = paste0(path_output_data, output_names[1], ".txt"))
vroom_write(gwas_traits_2, file = paste0(path_output_data, output_names[2], ".txt"))

writeLines("\n\n")
print(paste0('End of Running'))

# ================================================================================
# write the details information of the computer and softwares related to this code 
# please keep this info since it might be needed for any scientist in the future!
sessionInfo()
# calculation the running time of the program
end_time <- Sys.time()
end_time - start_time
# ================================================================================


#!/usr/bin/env Rscript

# 					EXTRACTION OF COMMON SNPs + Z calculation
#			 			OUTPUT will be used for running PLACO 
# ================================================================================
# Summary: The code reads a pair of GWAS data and extract common SNPs and stores
# the rows in the separated output files
# ================================================================================
# Written first by: Yazdan Asgari & Pierre-Emmanuel Sugier
# Modified by: 
# Initial Creation Date: 12/2020
# Edited Date: 09/2022
# ================================================================================

# for calculation the running time of the program
start_time <- Sys.time()

# ================================================================================
# libraries used in the code
# please install the following packages if needed
#BiocManager::install("vroom")
#BiocManager::install("dplyr")
#BiocManager::install("data.table")
library(vroom)
library(tidyverse)
library(data.table)
library("optparse")

# ================================================================================
# Read of command inputs
# ================================================================================

option_list = list(
  make_option(c("--path"), type="character", default=NULL, 
              help="directory path", metavar="character"),
  make_option(c("--inputs"), type="character", default=NULL, 
              help="name of file to read for inputs", metavar="character"),
  make_option(c("--out"), type="character", default="step1_out1", 
              help="prefix file name for output [default= %default]", metavar="character"),
  make_option(c("--info"), type="numeric", default=NULL, 
              help="threshold used for minimum imputation quality of SNPs", metavar="number"),
  make_option(c("--maf"), type="numeric", default=NULL, 
              help="threshold used for minimum MAF", metavar="number")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# opt <- data.frame(pathin=NA, inputs=NA, out=NA, info=NA, maf=NA)
# opt$pathin <- "/home1/6_AMLAP/PE/GCPBayes/scripts/Multipheno/test/"
# opt$inputs <- "readinputs.txt"
# opt$pathout <- "/home1/6_AMLAP/PE/GCPBayes/scripts/Multipheno/test/"
# opt$out <- "step1_output"
# opt$info <- 0.9
# opt$maf <- 0.05
#

# ================================================================================
# 						DEFINITION SECTION
# should be changed by a user
# ================================================================================

inputs <- read.table(file=paste0(opt$path,opt$inputs),header=TRUE)
nstud <- dim(inputs)[1]

# directory in which an output data would be written
path_output_data = opt$path

# the file names of the output data 
output_names <- opt$out

# used for filtering GWAS data              
info_threshold <- opt$info
MAF_threshold <- opt$maf



# ================================================================================
# 						            RUNNING SECTIONS
# ================================================================================
# Summary: 
# Reading GWAS data, filtering some SNPs, calculation of Z based on "beta" and "se"
# Extraction of common SNPs between pair of GWAS data
# ================================================================================


# treat_1: very fast reading method
gwas_traits_1 <- vroom(file=paste0(inputs$path[1],inputs$filename[1]))
gwas_traits_1 <- as.data.frame(gwas_traits_1)

# filtering step
gwas_traits_1 <- gwas_traits_1 %>% filter(gwas_traits_1$info > info_threshold & gwas_traits_1$MAF > MAF_threshold) %>% select(c("snp","chr","bp_hg19","Effect_A","Baseline_A","beta","se","pval","info","EAF","MAF"))

# calculation of Z
# if se==0, put a VERY LARGE VALUE for Z (here 1.0e+08)
gwas_traits_1 %>% mutate(Z = case_when(gwas_traits_1$se == 0   ~ 1.0e+08 , gwas_traits_1$se != 0 ~ (gwas_traits_1$beta/gwas_traits_1$se) )) -> gwas_traits_1


################################  
# to read other studies to get pool of SNPs to keep

for (i in 2:nstud){
# treat_2: very fast reading method
    gwas_traits_2 <- vroom(file=paste0(inputs$path[i],inputs$filename[i]))
    gwas_traits_2 <- as.data.frame(gwas_traits_2)

# filtering step
    gwas_traits_2 <- gwas_traits_2 %>% filter(gwas_traits_2$info > info_threshold & gwas_traits_2$MAF > MAF_threshold) %>% select(c("snp","chr","bp_hg19","Effect_A","Baseline_A","beta","se","pval","info","EAF","MAF"))

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
print(inputs$shortname[1])
writeLines("\n\n")
                                
# sorting based on the ID column (almost fast: about 30 sec)
output_traits_1 <- as.data.table(gwas_traits_1)
output_traits_1 <- setorder(output_traits_1, snp)
# writing the output file 1 (very fast)
vroom_write(output_traits_1, file = paste0(path_output_data, output_names, "_", inputs$shortname[1], ".txt"))
                          
                                
################################                     
# To read other studies to get corresponding outputs
                     
for (i in 2:nstud){

    print(inputs$shortname[i])
    writeLines("\n\n")
# treat_2: very fast reading method
    gwas_traits_2 <- vroom(file=paste0(inputs$path[i],inputs$filename[i]))
    gwas_traits_2 <- as.data.frame(gwas_traits_2)

# filtering step
    gwas_traits_2 <- select(filter(gwas_traits_2, gwas_traits_2$info > info_threshold & gwas_traits_2$MAF > MAF_threshold),
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
    vroom_write(output_traits_2, file = paste0(path_output_data, output_names, "_", inputs$shortname[i], ".txt")) 
}

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


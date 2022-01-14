# ================================================================================
#             *** SCRIPT TO CREATE GCPBAYES PIPELINE INPUT (5 STEPS) ***		 
# ================================================================================
# Summary: 
# Creates an input of the GCPBayes Package in 5 steps
# Starting from "reformatted GWAS Summary Statistics Data"
# Also needs "Annotation" and "LD Clumping" files as INPUT
# ================================================================================
# Written by: Yazdan Asgari
# Initial Creation Date: 07/2021
# Edited Date: 1/2022
# https://cesp.inserm.fr/en/equipe/exposome-and-heredity
# ================================================================================
# for calculation the running time of the program
start_time <- Sys.time()

# ================================================================================
# 						DEFINITION SECTION
# IMPORTANT NOTE: should be changed by a user
# ================================================================================
# directory in which your GWAS summary statistics data exist
path_gwas_sum_stat_bcac = "~/BCAC_OCAC/"
path_gwas_sum_stat_OCAC = "~/BCAC_OCAC/"

# the file names of your input data (used for step #1)
traits_1 <- "BCAC_2020_onco_ALL_reformatted.txt"
traits_2 <- "OCAC_BCAC_2020_onco_ALL_reformatted.txt"


# First GWAS Annotation file (used for step #2)
# directory in which the GWAS annotation data exist
path_annot_data = "~/BCAC_OCAC/"
# the file names of the GWAS annotation data
gwas_annot_file <- "Annot_BCAC_2020_onco_ALL_reformatted_coding.txt"

# Second GWAS Annotation file (used for step #5)
# For Adding Extra Information as new columns to GCPBayes input file
# directory in which the main annotation data exist
path_Annot = "~/BCAC_OCAC/"
# the file names of the annotation data
Annot <- "annot_gencode_v38lift37_modified_gene_class.txt"

# For LD Clumping (used for step #4)
# directory in which the LD data exist
LD_file_path = "~/BCAC_OCAC/"
# the file names of the LD data
LD_file <- "output_ld_clumping_08_BCAC_2020_ALL_OCAC.txt"


# Path for saving/reading files created by this script 
# IMPORTANT NOTE: These two paths ("path_input_data" & "path_output_data") MUST be defined as the same 
path_input_data = "~/BCAC_OCAC/"
path_output_data = "~/BCAC_OCAC/"

# Definition of some file names used during running the script in different steps
# used for steps #1 & #2 & #3 
output_shared <- c("output_pipeline_BCAC_ALL_Shared_OCAC_coding_clumping_08",
                   "output_pipeline_OCAC_Shared_BCAC_ALL_coding_clumping_08")

# used for steps #2 & #3
output_names <- "output_pipeline_BCAC_ALL_OCAC_coding_clumping_08"

# used for steps #3 & #4 & #5
output_names_summary <- c("output_pipeline_Summary_SNP_in_genes_BCAC_ALL_OCAC_coding_clumping_08",
                          "output_pipeline_Summary_SNP_in_genes_OCAC_BCAC_ALL_coding_clumping_08")


# used for filtering GWAS data (used for step #1)
info_threshold <- 0.9
MAF_threshold <- 0.05

# ================================================================================
# 						REQUIRED LIBRARIES
# Summary:
# libraries used in the code
# please install the following packages if needed
# ================================================================================
#BiocManager::install("vroom")
#BiocManager::install("dplyr")
#BiocManager::install("data.table")
#BiocManager::install("tidyr")
#BiocManager::install("tidyverse")
library(vroom)
library(dplyr)
library(data.table)
library(tidyr)
library(tidyverse)

# ================================================================================
# 						RUNNING SECTION
# ================================================================================
# 						   ### 1 ###
# Summary:
# keeping shared SNPs data between two traits
# NOTE: Please check column names to be sure if they are the same as the input data.
# ================================================================================
# reading first trait 
gwas_BCAC <- vroom(file=paste0(path_gwas_sum_stat_bcac, traits_1))
gwas_BCAC <- as.data.frame(gwas_BCAC)

# filtering step
gwas_BCAC <- select(filter(gwas_BCAC, gwas_BCAC$info > info_threshold & gwas_BCAC$MAF > MAF_threshold),
                    c("snp","chr","bp_hg19","Effect_A","nonEffect_A","beta","se","pval","info","MAF"))

# reading second trait 
gwas_OCAC <- vroom(file=paste0(path_gwas_sum_stat_OCAC, traits_2))
gwas_OCAC <- as.data.frame(gwas_OCAC)

# filtering step
gwas_OCAC <- select(filter(gwas_OCAC, gwas_OCAC$info > info_threshold & gwas_OCAC$MAF > MAF_threshold),
                    c("snp","chr","bp_hg19","Effect_A","nonEffect_A","beta","se","pval","info","MAF"))

# finding shared SNPs between two traits
SNP_shared <- intersect(gwas_BCAC$snp, gwas_OCAC$snp)

# extracting shared SNPs for the first trait
gwas_BCAC <- select(filter(gwas_BCAC, is.element(gwas_BCAC$snp, SNP_shared)),
                    c("snp","chr","bp_hg19","Effect_A","nonEffect_A","beta","se","pval","info","MAF"))

# sorting based on the ID column (almost fast: about 30 sec)
gwas_BCAC <- as.data.table(gwas_BCAC)
gwas_BCAC <- setorder(gwas_BCAC, snp)

# extracting shared SNPs for the second  trait
gwas_OCAC <- select(filter(gwas_OCAC, is.element(gwas_OCAC$snp, SNP_shared)),
                    c("snp","chr","bp_hg19","Effect_A","nonEffect_A","beta","se","pval","info","MAF"))

# sorting based on the ID column (almost fast: about 30 sec)
gwas_OCAC <- as.data.table(gwas_OCAC)
gwas_OCAC <- setorder(gwas_OCAC, snp)

# writing the output files 
vroom_write(gwas_BCAC, file = paste0(path_output_data, output_shared[1], ".txt"))
vroom_write(gwas_OCAC, file = paste0(path_output_data, output_shared[2], ".txt"))

writeLines("\n\n")
print("End of Running Section:#1")
writeLines("\n\n")

# ================================================================================
# 						   ### 2 ###
# Summary:
# Comparing an input file and the annotation file and adding the gene column corresponding to each SNP row 
# NOTE: Please check column names to be sure if they are the same as the input data.
# ================================================================================
# reading the annotation file 
gwas_annot <- vroom(file=paste0(path_annot_data, gwas_annot_file))
gwas_annot <- as.data.frame(gwas_annot)

# Because both files (output_shared[1] and output_shared[2]) have similar SNPs, choosing just one of them is enough for this step
gwas_input_file_2 <- paste0(output_shared[1], ".txt")

# reading input files one by one 
gwas_input <- vroom(file=paste0(path_input_data, gwas_input_file_2))
gwas_input <- as.data.frame(gwas_input)

# finding shared rows between an input file and the annotation file and adding the gene name correspond to the SNP
gwas_shared_gene_added <- inner_join(gwas_input, gwas_annot, by = c ("snp" = "SNP"))

# removing the rows in which no gene was assigned to
gwas_shared_gene_added <- filter(gwas_shared_gene_added, ANNOT != ".")

# extracting specific columns from the "gwas_shared_gene_added" 
gwas_shared_gene_added <- select(gwas_shared_gene_added, ANNOT, snp, chr, bp_hg19)

# changing the column names
colnames(gwas_shared_gene_added) <- c("gene", "SNP", "chr", "pos")

# writing the output file 
vroom_write(gwas_shared_gene_added, file = paste0(path_output_data, "output_pipeline_SNP_in_genes_", output_names, ".txt"), 
            col_names = TRUE, quote = "none")

writeLines("\n\n")
print("End of Running Section:#2")
writeLines("\n\n")

# ================================================================================
# 						   ### 3 ###
# Summary:
# Merging the Information From the two GWAS Summaries SNP files ("gwas_input_file_3[1]" & "gwas_input_file_3[2]") and 
# the file created at the end of the second section ("output_pipeline_SNP_in_genes_...txt") 
# ================================================================================
# the name of the files which created in the FIRST STEP OF THE GCPBAYES PIPELINE
gwas_input_file_3 = c()
gwas_input_file_3[1] <- paste0(output_shared[1], ".txt")
gwas_input_file_3[2] <- paste0(output_shared[2], ".txt")

# the name of the file created in the #2 step of this code
gwas_SNP_in_gene_iput_file <- paste0("output_pipeline_SNP_in_genes_", output_names, ".txt")

# reading the "gwas_shared_gene_added" file including gene name column (created in the #2 step of this code)
gwas_gene_input <- vroom(file=paste0(path_input_data, gwas_SNP_in_gene_iput_file))
gwas_gene_input <- as.data.frame(gwas_gene_input)

# ================================================================================
for (i in 1:length(gwas_input_file_3)) {
  
  # reading gwas summary statistics files 
  gwas_3 <- vroom(file=paste0(path_input_data, gwas_input_file_3[i]))
  gwas_3 <- as.data.frame(gwas_3)
  
  gwas_gene_added <- left_join(gwas_3, gwas_gene_input, by = c ("snp" = "SNP"))
  
  # remove the rows in which no gene was assigned to
  gwas_gene_added <- filter(gwas_gene_added, gene != "NA")
  
  # extracting specific columns from the "gwas_gene_added" 
  gwas_gene_added_select_col_3 <- select(gwas_gene_added, gene, snp, chr.x, pos, Effect_A, nonEffect_A, beta, se, pval)
  
  # writing the output files 
  if (i%%2 != 0){
    vroom_write(gwas_gene_added_select_col_3, file = paste0(path_output_data, output_names_summary[1], ".txt"), 
                col_names = TRUE, quote = "none")
  }
  if (i%%2 == 0){
    vroom_write(gwas_gene_added_select_col_3, file = paste0(path_output_data, output_names_summary[2], ".txt"), 
                col_names = TRUE, quote = "none")
  }
  
}

writeLines("\n\n")
print("End of Running Section:#3")
writeLines("\n\n")

# ================================================================================
# 						   ### 4 ###
# Summary:
# Selection of independant SNPs which are available in the LD file
# ================================================================================

# reading the LD file 
LD <- vroom(file=paste0(LD_file_path, LD_file), col_names = FALSE)
LD <- as.data.frame(LD)

# changing the column names and class of the data
colnames(LD) <- c("snp","chr","pos_hg19","Effect_A","nonEffect_A","pval","ldclump_id")

gwas_gene_added_select_col_4 = c()
gwas_gene_added_select_col_4[1] <- paste0(output_names_summary[1], ".txt")
gwas_gene_added_select_col_4[2] <- paste0(output_names_summary[2], ".txt")

for (i in 1:length(gwas_gene_added_select_col_4)) {
  
  # reading gwas summary statistics file
  gwas_4 <- vroom(file=paste0(path_output_data, gwas_gene_added_select_col_4[i]))
  gwas_4 <- as.data.frame(gwas_4)
  
  # finding shared rows between a gwas summary statistics file and the the LD file
  gwas_genes_pruned <- inner_join(gwas_4, LD, by = "snp")
  
  # writing the output files 
  if (i%%2 != 0){
    vroom_write(gwas_genes_pruned, file = paste0(path_output_data, output_names_summary[1],"_indep_SNPs.txt"), 
                col_names = TRUE, quote = "none")
  }
  if (i%%2 == 0){
    vroom_write(gwas_genes_pruned, file = paste0(path_output_data, output_names_summary[2],"_indep_SNPs.txt"), 
                col_names = TRUE, quote = "none")
  }
  
}

writeLines("\n\n")
print("End of Running Section:#4")
writeLines("\n\n")

# ================================================================================
# 						   ### 5 ###
# Summary:
# creation of the final files in ".Rdata" format which would be used as input data for running GCPBayes
# Reading two GWAS files and splitting rows in which more than one gene assigned to a SNP (separated by "|")
# and removing other information created during the annotation step (rows with "missense/nonsense/frameshift/splice")
# ================================================================================
# obtaining the file names of your input data (created at the end of fourth section)
gwas_5 = c()
gwas_5[1] <- paste0(path_input_data, output_names_summary[1], "_indep_SNPs.txt")
gwas_5[2] <- paste0(path_input_data, output_names_summary[2], "_indep_SNPs.txt")

# reading GWAS data for the first trait1
gwas5_1 <- vroom(file=gwas_5[1])
gwas5_1 <- as.data.frame(gwas5_1)

# removing potential duplicated entries
gwas5_1 <- gwas5_1[!duplicated(gwas5_1$snp), ]

# Splitting rows contain more than one gene
gwas5_1 <- separate_rows(gwas5_1,gene, sep="\\|")

# Delete all characters started with "(" and end with ")" 
gwas5_1$gene  <- gsub("\\(.*\\)","" , gwas5_1$gene , ignore.case = TRUE)

# Removing rows include "missense/nonsense/frameshift/splice" (lowercase or UPPERCASE) as gene name
gwas5_1 <- gwas5_1[-which(gwas5_1$gene=="=missense"), ]
gwas5_1 <- gwas5_1[-which(gwas5_1$gene=="=MISSENSE"), ]
gwas5_1 <- gwas5_1[-which(gwas5_1$gene=="=nonsense"), ]
gwas5_1 <- gwas5_1[-which(gwas5_1$gene=="=NONSENSE"), ]
gwas5_1 <- gwas5_1[-which(gwas5_1$gene=="=frameshift"), ]
gwas5_1 <- gwas5_1[-which(gwas5_1$gene=="=FRAMESHIFT"), ]
gwas5_1 <- gwas5_1[-which(gwas5_1$gene=="=splice"), ]
gwas5_1 <- gwas5_1[-which(gwas5_1$gene=="=SPLICE"), ]

# reading GWAS data for the second trait1
gwas5_2 <- vroom(file=gwas_5[2])
gwas5_2 <- as.data.frame(gwas5_2)

# removing potential duplicated entries
gwas5_2 <- gwas5_2[!duplicated(gwas5_2$snp), ]

# Splitting rows contain more than one gene
gwas5_2 <- separate_rows(gwas5_2,gene, sep="\\|")

# Delete all characters started with "(" and end with ")" 
gwas5_2$gene  <- gsub("\\(.*\\)","" , gwas5_2$gene , ignore.case = TRUE)

# Removing rows include "missense/nonsense/frameshift/splice" (lowercase or UPPERCASE) as gene name
gwas5_2 <- gwas5_2[-which(gwas5_2$gene=="=missense"), ]
gwas5_2 <- gwas5_2[-which(gwas5_2$gene=="=MISSENSE"), ]
gwas5_2 <- gwas5_2[-which(gwas5_2$gene=="=nonsense"), ]
gwas5_2 <- gwas5_2[-which(gwas5_2$gene=="=NONSENSE"), ]
gwas5_2 <- gwas5_2[-which(gwas5_2$gene=="=frameshift"), ]
gwas5_2 <- gwas5_2[-which(gwas5_2$gene=="=FRAMESHIFT"), ]
gwas5_2 <- gwas5_2[-which(gwas5_2$gene=="=splice"), ]
gwas5_2 <- gwas5_2[-which(gwas5_2$gene=="=SPLICE"), ]

# ================================================================================
# Creation of "extra_info" list including some columns from the "Annotation file"
# ================================================================================
# Reading the main "annotation file"
Annot <- read.csv(paste0(path_Annot, Annot), header = T, sep = "\t", dec = ".")

# removing first three characters from the column
x <- substring(Annot[,1], 4, nchar(as.character(Annot[,1])))
Annot$chr <- x

# replacing characters into numbers for "chromosome index"
Annot$chr <- gsub("X", "23", Annot$chr)
Annot$chr <- gsub("x", "23", Annot$chr)
Annot$chr <- gsub("Y", "24", Annot$chr)
Annot$chr <- gsub("y", "24", Annot$chr)
Annot$chr <- gsub("M", "25", Annot$chr)
Annot$chr <- gsub("m", "25", Annot$chr)

Annot$chr <- as.numeric(Annot$chr)

# merging gwas summary statistics of one of the traits with the "Annot" file
# Why one of the traits: Because both trats now have similar SNPs in row 
gwas5_extra_annot_columns <- inner_join(gwas5_1, Annot, by = c ("gene" = "gene_name", "chr"))

# picking of the selected columns
gwas5_extra_annot_columns_selected <- gwas5_extra_annot_columns[ , c("snp", "chr", "pos_hg19", 
                                                                     "Effect_A.y", "nonEffect_A.y", "beta", "se",
                                                                     "gene", "start", "end", "gene_type", "seqnames")]

colnames(gwas5_extra_annot_columns_selected) <- c("snp", "chr", "pos_hg19", 
                                                  "Effect_A", "nonEffect_A", "beta", "se",
                                                  "gene", "start", "end", "gene_type", "seqnames")

# sorting data based on "gene" column
gwas5_merge_extra_sorted <- setorder(gwas5_extra_annot_columns_selected, gene)
gwas5_merge_extra_sorted <- as.data.frame(gwas5_merge_extra_sorted)

# splitting the data frame based on the "gene" column
gwas5_merge_extra_sorted_tolist <- gwas5_merge_extra_sorted %>% split(f=gwas5_merge_extra_sorted$gene)

# creation of "gcpbayes_extra_input_data_final" (as a list) 
# It would be used for adding some info columns to the "result" file of the GCPBayes output
f_build_extra_inputs = function(data){
  gene <- as.character(data$gene[1])
  snp_number <- as.numeric(length(data$snp))
  chr <- as.numeric(data$chr[1])
  pos <- as.numeric(data$pos[1])
  start <- as.numeric(data$start[1])
  end <- as.numeric(data$end[1])
  gene_type <- as.character(data$gene_type[1])
  gene_length <- as.numeric(end - start)
  snp_number_gene_length_ratio <- as.numeric(snp_number / gene_length)
  extra_info <- list(gene=gene, snp_number=snp_number, chr=chr, pos=pos,
                     start=start, end=end, gene_type=gene_type,
                     gene_length=gene_length, snp_number_gene_length_ratio=snp_number_gene_length_ratio)
}

# applying the "f_build_extra_inputs" function for creation of "gcpbayes_extra_input_data_final" (as a list)
gcpbayes_extra_input_data_final <- gwas5_merge_extra_sorted_tolist %>% lapply(FUN=f_build_extra_inputs)

# sorting the "gcpbayes_extra_input_data_final" list based on length of SNPs in the Gene (in an increasing order)
gcpbayes_extra_input_data_final <- gcpbayes_extra_input_data_final[order(sapply(seq(1:length(gcpbayes_extra_input_data_final)), function(x) gcpbayes_extra_input_data_final[[x]]$snp_number), decreasing = FALSE)]

# saving the "gcpbayes_extra_input_data_final" output as a ".Rdata" file format
save(gcpbayes_extra_input_data_final,file=paste0(path_output_data, "Matrices_extra_info_", output_names, ".Rdata"))

# ================================================================================
# creation of "gcpbayes_input_data_final" (as a list)
# ================================================================================
# merging GWAS summary statistics of two traits
gwas5_merge <- inner_join(gwas5_1, gwas5_2, by = c("gene","snp","pos"))

# sorting data based on "gene" column
gwas5_merge_sorted <- setorder(gwas5_merge, gene)
gwas5_merge_sorted <- as.data.frame(gwas5_merge_sorted)

# extraction of columns needed for creation of "GCPBayes" input ("gene", "snp", "beta and se" for both traits)
gwas5_merge_sumstat <- gwas5_merge_sorted %>% select(c(gene, snp, beta.x, beta.y, se.x, se.y))

# renaming the extracted columns
colnames(gwas5_merge_sumstat) <- c("Gene", "SNP", "beta1", "beta2", "se1", "se2")

# splitting the data frame based on the "Gene" column
gwas5_merge_sumstat_tolist <- gwas5_merge_sumstat %>% split(f=gwas5_merge_sumstat$Gene)

# defining a function for creation of "gcpbayes_input_data_final" (as a list)
f_build_inputs = function(data){
  b1 <- data$beta1
  b2 <- data$beta2
  s1 <- diag(data$se1*data$se1, length(data$se1), length(data$se1))
  s2 <- diag(data$se2*data$se2, length(data$se2), length(data$se2))
  names(b1) <- data$SNP
  names(b2) <- data$SNP
  dimnames(s1) <- list(data$SNP, data$SNP)
  dimnames(s2) <- list(data$SNP, data$SNP)
  Betah <- list(b1, b2)
  Sigmah <- list(diag(data$se1), diag(data$se2))
  list(Betah=list(b1, b2), Sigmah=list(s1, s2))
}

# applying the "f_build_inputs" function for creation of "gcpbayes_input_data_final" (as a list)
gcpbayes_input_data_final <- gwas5_merge_sumstat_tolist %>% lapply(FUN=f_build_inputs)

# sorting the "gcpbayes_input_data_final" list based on length of SNPs in the Gene (in an increasing order)
gcpbayes_input_data_final <- gcpbayes_input_data_final[order(sapply(seq(1:length(gcpbayes_input_data_final)), function(x) length(gcpbayes_input_data_final[[x]]$Betah[[1]])), decreasing = FALSE)]

# saving the "gcpbayes_input_data_final" output as a ".Rdata" file format
save(gcpbayes_input_data_final, file=paste0(path_output_data, "Matrices_", output_names, ".Rdata"))

writeLines("\n\n")
print("End of Running Section:#5")

# ================================================================================
# write the details information of the computer and softwares related to this code 
# please keep this info since it might be needed for any scientist in the future!
sessionInfo()
# calculation the running time of the program
end_time <- Sys.time()
end_time - start_time
# ================================================================================

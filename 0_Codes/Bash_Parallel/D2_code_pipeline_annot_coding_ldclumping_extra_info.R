#!/usr/bin/env Rscript

#             *** SCRIPT TO CREATE GCPBAYES PIPELINE INPUT ***		 
# ================================================================================
# Summary: 
# Creates an input of the GCPBayes Package in 5 steps
# Starting from "reformatted GWAS Summary Statistics Data"
# Also needs "Annotation" and "LD Clumping" files as INPUT
# ================================================================================
# Written by: Yazdan Asgari
# Modified by: Yazdan Asgari
# Initial Creation Date: 07/2021
# Edited Date: 12/2022
# ================================================================================

# for calculation the running time of the program
start_time <- Sys.time()

# ================================================================================
# 						REQUIRED LIBRARIES
# Summary:
# libraries used in the code
# please install the following packages if needed
# ================================================================================
library(vroom)
library(data.table)
library(tidyverse)
library(optparse)

# ================================================================================
# Read of command inputs
# ================================================================================

option_list = list(
  make_option(c("--pathin"), type="character", default=NULL, 
              help="working directory", metavar="character"),
  make_option(c("--inputs"), type="character", default=NULL, 
              help="name of file to read for inputs", metavar="character"),
  make_option(c("--annot_gwas_path"), type="character", default=NULL, 
              help="directory in which the GWAS annotation data exist", metavar="character"),
  make_option(c("--annot_gwas_file"), type="character", default=NULL, 
              help="the file name of the GWAS annotation data", metavar="character"),			  
  make_option(c("--annot_path"), type="character", default=NULL, 
              help="directory path for the annotations", metavar="character"),
  make_option(c("--annot_file"), type="character", default=NULL, 
              help="annotation file name", metavar="character"),
  make_option(c("--ldclumpout"), type="character", default=NULL, 
              help="LD clumping output file", metavar="character"),
  make_option(c("--info"), type="numeric", default=NULL, 
              help="threshold used for minimum imputation quality of SNPs", metavar="number"),
  make_option(c("--maf"), type="numeric", default=NULL, 
              help="threshold used for minimum MAF", metavar="number"),
  make_option(c("--pathout"), type="character", default=NULL, 
              help="directory path for output", metavar="character"),
  make_option(c("--out"), type="character", default="step1_out1.txt", 
              help="generic output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# ================================================================================
# 						DEFINITION SECTION
# ================================================================================

# to read paths, names, and shortnames of the inputs / and number of studies to consider
inputs <- read.table(file=paste0(opt$pathin,opt$inputs),header=TRUE)
nstud <- dim(inputs)[1]

# IMPORTANT NOTE: These to PATHS MUST be defined as the same 
# ("path_input_data" & "path_output_data" for reading/writing the input and output files) 
path_input_data = inputs$path[1] %>% as.character
path_output_data = opt$pathout

# used for steps #1 & #2 & #3
output_common <- c()
for (i in 1:nstud){output_common[i] <- paste("D2_output_pipeline",inputs$shortname[i],"common2all",nstud,"studies","coding_clumping",sep="_")}

# used for steps #2 & #3
output_names <- paste("output_pipeline",opt$out,sep="_")

# used for steps #3 & #4
output_names_summary <- c()
for (i in 1:nstud){output_names_summary[i] <- paste("D2_Summary_SNP_in_genes","output_pipeline",inputs$shortname[i],"common2all",nstud,"studies","coding_clumping",sep="_")}

# For GWAS Annotation (used for step #2)
# directory in which the GWAS annotation data exist
path_annot_data = opt$annot_gwas_path
# the file name of the GWAS annotation data
gwas_annot_file <- opt$annot_gwas_file

# For LD Clumping (used for step #4)
# directory in which the LD data exist
LD_file_path = opt$pathin
# the file names of the LD data
LD_file <- opt$ldclumpout

# For Adding Extra Information as new columns to GCPBayes input file (used for step #5)
# directory in which the main annotation data exist
path_Annot = opt$annot_path
# the file names of the annotation data
file_Annot <- opt$annot_file

# used for extrating SNPs with value greater than these (used for step #1)
info_threshold <- opt$info
MAF_threshold <- opt$maf

# ================================================================================
# 						RUNNING SECTION
# ================================================================================
# 						   ### 1 ###
# Summary:
# keeping common SNPs data between two traits
# NOTE: Please check column names to be sure if they are the same as the input data.
# ================================================================================

# reading first trait / very fast reading method
gwas_traits_1 <- vroom(file=paste0(inputs$path[1],inputs$filename[1]))
gwas_traits_1 <- as.data.frame(gwas_traits_1)

# filtering step
gwas_traits_1 <- gwas_traits_1 %>% filter(gwas_traits_1$info > info_threshold & gwas_traits_1$MAF > MAF_threshold) %>% select(c("snp","chr","bp_hg19","Effect_A","Baseline_A","beta","se","pval","info","EAF","MAF"))

# reading other traits
for (i in 2:nstud){
# treat_2: very fast reading method
    gwas_traits_others <- vroom(file=paste0(inputs$path[i],inputs$filename[i]))
    gwas_traits_others <- as.data.frame(gwas_traits_others)
# filtering step
    gwas_traits_others <- gwas_traits_others %>% filter(gwas_traits_others$info > info_threshold & gwas_traits_others$MAF > MAF_threshold) %>% select(c("snp","chr","bp_hg19","Effect_A","Baseline_A","beta","se","pval","info","EAF","MAF"))
# merging two GWAS data based on "snp", "chr", and "bp_hg19" columns
    gwas_merge <- inner_join(gwas_traits_1,gwas_traits_others, by = c("snp", "chr", "bp_hg19"))
# checking number of total, unique, and duplicated rsids
    length(gwas_merge$snp)
    length(unique(gwas_merge$snp))
    dim(gwas_merge[duplicated(gwas_merge$snp), ])[1]
# removing potential duplicated entries
    gwas_merge <- gwas_merge[!duplicated(gwas_merge$snp), ] # NOTE: Only based on snp id? Maybe this would be better to delete all rows for the SNPs and not only duplicates to avoid merging non corresponding alleles?
# checking data after removing the duplicated entries
    length(gwas_merge$snp)
    length(unique(gwas_merge$snp))
    dim(gwas_merge[duplicated(gwas_merge$snp), ])[1]
# extracting common SNPs for the first trait
    gwas_traits_1 <- gwas_merge[, c("snp", "chr", "bp_hg19", 
                                "Effect_A.x", "Baseline_A.x", "beta.x", "se.x", 
                                "pval.x", "info.x", "EAF.x", "MAF.x") ]
    colnames(gwas_traits_1) <- c("snp","chr","bp_hg19","Effect_A","Baseline_A","beta","se","pval","info","EAF","MAF")
} 

# sorting based on the ID column (almost fast: about 30 sec)
output_traits_1 <- as.data.table(gwas_traits_1)
output_traits_1 <- setorder(output_traits_1, snp)
# writing the output file 1 (very fast)
vroom_write(output_traits_1, file = paste0(path_output_data, output_common[1],".txt"))

################################                     
# To read again other studies to get corresponding outputs

for (i in 2:nstud){
    print(inputs$shortname[i])
    writeLines("\n\n")
# treat_2: very fast reading method
    gwas_traits_others <- vroom(file=paste0(inputs$path[i],inputs$filename[i]))
    gwas_traits_others <- as.data.frame(gwas_traits_others)
# filtering step
    gwas_traits_others <- select(filter(gwas_traits_others, gwas_traits_others$info > info_threshold & gwas_traits_others$MAF > MAF_threshold),
                        c("snp","chr","bp_hg19","Effect_A","Baseline_A","beta","se","pval","info","EAF","MAF"))
# merging two GWAS data based on "snp", "chr", and "bp_hg19" columns
    gwas_merge <- inner_join(gwas_traits_1,gwas_traits_others, by = c("snp", "chr", "bp_hg19"))               
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
    gwas_traits_others <- gwas_merge[, c("snp", "chr", "bp_hg19", 
                                "Effect_A.y", "Baseline_A.y", "beta.y", "se.y", 
                                "pval.y", "info.y", "EAF.y", "MAF.y") ]
    colnames(gwas_traits_others) <- c("snp", "chr", "bp_hg19", 
                             "Effect_A", "Baseline_A", "beta", "se", 
                             "pval", "info", "EAF", "MAF") 
# sorting based on the ID column (almost fast: about 30 sec)
    output_traits_others <- as.data.table(gwas_traits_others)
    output_traits_others <- setorder(output_traits_others, snp)
# writing the output files (very fast)
    vroom_write(output_traits_others, file = paste0(path_output_data, output_common[i],".txt"))
}

writeLines("\n\n")
print("End of Running Section:#1")
writeLines("\n\n")

# ================================================================================
# 						   ### 2 ###
# Summary:
# Comparing an input file and the annotation file and adding the gene coulumn corresponding to each SNP row 
# NOTE: Please check column names to be sure if they are the same as the input data.
# ================================================================================
# reading the annotation file 
gwas_annot <- vroom(file=paste0(path_annot_data, gwas_annot_file,".txt")) %>% as.data.frame

# Because both files (output_common[1] and output_common[2]) have similar SNPs, choosing just one of them is enough for this step
gwas_input_step2 <- vroom(file=paste0(path_output_data, output_common[i],".txt"))
gwas_input_step2 <- as.data.frame(gwas_input_step2)

# merging input file to annotation file & adding the corresponding gene name, & removing rows not assigned to any gene, & extracting desired columns
gwas_common_gene_added <- inner_join(gwas_input_step2, gwas_annot, by = c ("snp" = "SNP")) %>% filter(ANNOT != ".") %>% select(ANNOT,snp,chr,bp_hg19)

# changing the column names
colnames(gwas_common_gene_added) <- c("gene", "SNP", "chr", "pos")

# writing the output file 
vroom_write(gwas_common_gene_added, file = paste0(path_output_data, "D2_output_pipeline_SNP_in_genes_", output_names, ".txt"), 
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
gwas_input_file_step3 = c()
for (i in 1:nstud){gwas_input_file_step3[i] <- paste0(output_common[i], ".txt")}

# reading the "gwas_common_gene_added" file including gene name column (created in the #2 step of this code)
gwas_gene_input <- vroom(file=paste0(path_output_data, "D2_output_pipeline_SNP_in_genes_", output_names, ".txt")) %>% as.data.frame

# ================================================================================
for (i in 1:nstud) {
  # reading gwas summary statistics files 
  gwas_step3 <- vroom(file=paste0(path_output_data, gwas_input_file_step3[i])) %>% as.data.frame
  
  # merging with annotated base, & removing the rows in which no gene was assigned to, & selection of desired columns
  gwas_gene_added_step3 <- left_join(gwas_step3, gwas_gene_input, by = c ("snp" = "SNP")) %>% filter(gene != "NA") %>% select(gene, snp, chr.x, pos, Effect_A, Baseline_A, beta, se, pval)
  colnames(gwas_gene_added_step3)[colnames(gwas_gene_added_step3) == "chr.x"] <- "chr"
  
  #writing the results
  vroom_write(gwas_gene_added_step3, file = paste0(path_output_data, output_names_summary[i], ".txt"), col_names = TRUE, quote = "none")
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
colnames(LD) <- c("snp","chr","pos_hg19","Effect_A","Baseline_A","pval","ldclump_id")

for (i in 1:nstud) {
  # reading gwas summary statistics file
  gwas_stepld <- vroom(file=paste0(path_output_data, output_names_summary[i], ".txt")) %>% as.data.frame
  # finding common rows between a gwas summary statistics file and the the LD file
  gwas_genes_pruned <- inner_join(gwas_stepld, LD, by = "snp")
  # writing the output files 
  vroom_write(gwas_genes_pruned, file = paste0(path_output_data, output_names_summary[i],"_indep_SNPs.txt"), 
                col_names = TRUE, quote = "none")
}

writeLines("\n\n")
print("End of Running Section:#3.5: LD step")
writeLines("\n\n")

# ================================================================================
# 						   ### 4 ###
# Summary:
# creation of the final files in ".Rdata" format which would be used as input data for running GCPBayes
# Reading two GWAS files and splitting rows in which more than one gene assigned to a SNP (separated by "|")
# and removing other information created during the annotation step (rows with "missense/nonsense/frameshift/splice")
# ================================================================================
# obtaining the file names of your input data (created at the end of fourth section)
gwas_step4_names = c()
for (i in 1:nstud){gwas_step4_names[i] <- paste0(path_output_data, output_names_summary[i], "_indep_SNPs.txt")}

# reading GWAS data for each trait
gwas_step4 <-list()
for (i in 1:nstud){
    gwas_step4[[i]] <- vroom(file=gwas_step4_names[i]) %>% as.data.frame
# removing potential duplicated entries & splitting rows that contain more than one gene
    gwas_step4[[i]] <- gwas_step4[[i]][!duplicated(gwas_step4[[i]]$snp), ]
# Splitting rows contain more than one gene
    gwas_step4[[i]] <- separate_rows(gwas_step4[[i]], gene, sep="\\|")
# Delete all characters started with "(" and end with ")" 
    gwas_step4[[i]]$gene <- gsub("\\(.*\\)", "", gwas_step4[[i]]$gene, ignore.case = TRUE)
# Removing rows including "missense/nonsense/frameshift/splice" (lowercase or UPPERCASE) as gene name
    if(length(which(gwas_step4[[i]]$gene=="=missense"))!=0){gwas_step4[[i]] <- gwas_step4[[i]][-which(gwas_step4[[i]]$gene=="=missense"), ]}
    if(length(which(gwas_step4[[i]]$gene=="=MISSENSE"))!=0){gwas_step4[[i]] <- gwas_step4[[i]][-which(gwas_step4[[i]]$gene=="=MISSENSE"), ]}
    if(length(which(gwas_step4[[i]]$gene=="=nonsense"))!=0){gwas_step4[[i]] <- gwas_step4[[i]][-which(gwas_step4[[i]]$gene=="=nonsense"), ]}
    if(length(which(gwas_step4[[i]]$gene=="=NONSENSE"))!=0){gwas_step4[[i]] <- gwas_step4[[i]][-which(gwas_step4[[i]]$gene=="=NONSENSE"), ]}
    if(length(which(gwas_step4[[i]]$gene=="=frameshift"))!=0){gwas_step4[[i]] <- gwas_step4[[i]][-which(gwas_step4[[i]]$gene=="=frameshift"), ]}
    if(length(which(gwas_step4[[i]]$gene=="=FRAMESHIFT"))!=0){gwas_step4[[i]] <- gwas_step4[[i]][-which(gwas_step4[[i]]$gene=="=FRAMESHIFT"), ]}
    if(length(which(gwas_step4[[i]]$gene=="=splice"))!=0){gwas_step4[[i]] <- gwas_step4[[i]][-which(gwas_step4[[i]]$gene=="=splice"), ]}
    if(length(which(gwas_step4[[i]]$gene=="=SPLICE"))!=0){gwas_step4[[i]] <- gwas_step4[[i]][-which(gwas_step4[[i]]$gene=="=SPLICE"), ]}
    colnames(gwas_step4[[i]]) <- c("gene", "snp", "chr", "pos", "Effect_A.x", "Baseline_A.x", "beta", "se",
                                                  "pval.x", "chr.y", "pos_hg19", "Effect_A.y", "Baseline_A.y","pval.y", "ldclump_id")
    gwas_step4[[i]] <- gwas_step4[[i]][ , c("snp", "chr.y", "pos_hg19", "Effect_A.y", "Baseline_A.y", "beta", "se", "pval.x",
                                                                     "gene", "pval.y", "ldclump_id")]
    colnames(gwas_step4[[i]]) <- c("snp", "chr", "pos", "Effect_A", "Baseline_A", "beta", "se", "pval",
                                                                     "gene", "pval_placo", "ldclump_id")
}

# ================================================================================
# Creation of "extra_info" list including some columns from the "Annotation file"
# ================================================================================
# Reading the main "annotation file"
Annot <- read.csv(paste0(path_Annot, file_Annot), header = T, sep = "\t", dec = ".")

# removing first three characters from the column
Annot$chr <- substring(Annot[,1], 4, nchar(as.character(Annot[,1])))

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
gwas_step4_extra_annot_columns <- inner_join(gwas_step4[[1]], Annot, by = c ("gene" = "gene_name", "chr"))

# picking of the selected columns
gwas_step4_extra_annot_columns_selected <- gwas_step4_extra_annot_columns[ , c("snp", "chr", "pos", 
                                                                     "Effect_A", "Baseline_A", "beta", "se",
                                                                     "gene", "start", "end", "gene_type", "seqnames")]
colnames(gwas_step4_extra_annot_columns_selected) <- c("snp", "chr", "pos_hg19", 
                                                                     "Effect_A", "Baseline_A", "beta", "se",
                                                                     "gene", "start", "end", "gene_type", "seqnames")

# sorting data based on "gene" column
gwas_step4_merge_extra_sorted <- setorder(gwas_step4_extra_annot_columns_selected, gene) %>% as.data.frame

# splitting the data frame based on the "gene" column
gwas_step4_merge_extra_sorted_tolist <- gwas_step4_merge_extra_sorted %>% split(f=gwas_step4_merge_extra_sorted$gene)

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
gcpbayes_extra_input_data_final <- gwas_step4_merge_extra_sorted_tolist %>% lapply(FUN=f_build_extra_inputs)

# sorting the "gcpbayes_extra_input_data_final" list based on length of SNPs in the Gene (in an increasing order)
gcpbayes_extra_input_data_final <- gcpbayes_extra_input_data_final[order(sapply(seq(1:length(gcpbayes_extra_input_data_final)), function(x) gcpbayes_extra_input_data_final[[x]]$snp_number), decreasing = FALSE)]

# saving the "gcpbayes_extra_input_data_final" output as a ".Rdata" file format
save(gcpbayes_extra_input_data_final,file=paste0(path_output_data, "D2_Matrices_extra_info_", output_names, ".Rdata"))

# ================================================================================
# creation of "gcpbayes_input_data_final" (as a list)
# ================================================================================
# merging GWAS summary statistics of all the traits
gwas_step5_merge <- gwas_step4[[1]] %>% select(c(gene, snp, pos, beta, se))
colnames(gwas_step5_merge) <- c("gene","snp","pos","beta1","se1") 

for (i in 2:nstud){
	tempo_clean_gwas <- gwas_step4[[i]] %>% select(c(gene, pos, snp, beta, se))
	gwas_step5_merge <- inner_join(gwas_step5_merge, tempo_clean_gwas, by = c("gene","snp","pos"))
	# ordering and renaming columns
	head_beta=list()
	head_se=list()
	for (j in 1:i){
		head_beta[[j]] <- paste0('beta',j)
		head_se[[j]] <- paste0('se',j)
	}
	gwas_step5_merge <- gwas_step5_merge[,c("gene","snp","pos",unlist(head_beta)[-i],"beta",unlist(head_se)[-i],"se")]
	colnames(gwas_step5_merge) <- c("gene","snp","pos",unlist(head_beta),unlist(head_se))
}

# sorting data based on "gene" column & extraction of columns needed for creation of "GCPBayes" input ("gene", "snp", "beta and se" for all traits)
gwas_step5_merge_sumstat <- setorder(gwas_step5_merge, gene) %>% select(-pos) %>% as.data.frame
# renaming the extracted columns
head_beta=list()
head_se=list()
for (i in 1:nstud){
	head_beta[[i]] <- paste0('beta',i)
	head_se[[i]] <- paste0('se',i)
}
colnames(gwas_step5_merge_sumstat) <- c("Gene", "SNP", unlist(head_beta),unlist(head_se)) 

# splitting the data frame based on the "Gene" column
gwas_step5_merge_sumstat_tolist <- gwas_step5_merge_sumstat %>% split(f=gwas_step5_merge_sumstat$Gene)

# defining a function for creation of "gcpbayes_input_data_final" (as a list)
f_build_inputs = function(data){
	n_studies <- (dim(data)[2]-2)/2 # to retrieve nstud
	Betah <- list()
	Sigmah <- list()
	for (i in 1:n_studies){
		n_col_beta <- 2 + i
		n_col_se <- 2 + n_studies + i
		Betah[[i]] <- data[,n_col_beta] %>% as.numeric
		names(Betah[[i]]) <- data$SNP
		Sigmah[[i]] <- diag(data[,n_col_se]*data[,n_col_se], length(data[,n_col_se]), length(data[,n_col_se]))
		dimnames(Sigmah[[i]]) <- list(data$SNP, data$SNP)
	}
	list(Betah=Betah, Sigmah=Sigmah)
}

# applying the "f_build_inputs" function for creation of "gcpbayes_input_data_final" (as a list)
gcpbayes_input_data_final <- gwas_step5_merge_sumstat_tolist %>% lapply(FUN=f_build_inputs)

# sorting the "gcpbayes_input_data_final" list based on length of SNPs in the Gene (in an increasing order)
gcpbayes_input_data_final <- gcpbayes_input_data_final[order(sapply(seq(1:length(gcpbayes_input_data_final)), function(x) length(gcpbayes_input_data_final[[x]]$Betah[[1]])), decreasing = FALSE)]

# saving the "gcpbayes_input_data_final" output as a ".Rdata" file format
save(gcpbayes_input_data_final, file=paste0(path_output_data, "D2_Matrices_", output_names, ".Rdata"))

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

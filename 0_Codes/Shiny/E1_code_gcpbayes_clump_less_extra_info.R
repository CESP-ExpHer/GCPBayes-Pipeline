# =====================================================================
#   ***Running GCPBayes for a pair of traits***
# =====================================================================
# Summary:
# Running GCPBayes for Genes which includes SNPs less than a Threshold
# It creates 3 different output files: "results", "pleiotropy", and "errors" 
# ================================================================================
# Initially Written by: PE Sugier, Elise Lucotte
# Modified by Yazdan Asgari
# Modified by: Yazdan Asgari
# Initial Creation Date: 11/2021
# Edited Date: 03/2023
# ================================================================================

# for calculation the running time of the program
start_time <- Sys.time()

# ================================================================================
# libraries used in the code
library(tictoc)
library(GCPBayes)
# library(BhGLM)

# ================================================================================
# 						DEFINITION SECTION
# ================================================================================
# PATH in which a GCPBayes input file exists
path_input_data_read <- path_input_data

# GCPBayes input file names
Matrix_input_file <- paste("D3_Matrices_output_pipeline_",output_step5c,".Rdata",sep="")
Matrix_extra_info_input_file <- paste("D3_Matrices_extra_info_output_pipeline_",output_step5c,".Rdata",sep="")

# PATH for writing output files
path_data_write = path_output_data

output_names <- c(paste(output_step6c,"less_threshold",group_clump_threshold,"errors",sep="_"),
				  paste(output_step6c,"less_threshold",group_clump_threshold,"results",sep="_"),
				  paste(output_step6c,"less_threshold",group_clump_threshold,"pleiotropy",sep="_"))

# Parameters used in this script
# A threshold for SNP numbers to split Genes based on that 
SNP_number_threshold <- group_clump_threshold
# A threshold for considering a Gene with a potential pleiotropic effect 
theta_threshold <- theta_exploration


# ================================================================================
# 						RUNNING SECTION
# ================================================================================
# Summary:
# Reading the loaded data
# Definition of parameter
# Running GCPBayes function for every Gene
# Saving 3 output files: "results", "pleiotropy", and "errors" 
# ================================================================================
# Loading Matrix data created by the pipeline (including genes and their SNPs data)
load(file=paste0(path_input_data_read, Matrix_input_file[1]))

# reading data from the loaded Matrix and preparing it for the splitting (next step)
len_Matrix <- seq(1, length(gcpbayes_input_data_final))
final_data_new <-as.data.frame(do.call(rbind, lapply(len_Matrix, function(x) length(gcpbayes_input_data_final[[x]]$Betah[[1]]))))
summary(unlist(final_data_new))

# splitting Genes based on their number of SNPs using a threshold (less/greater than threshold)
S_less_count <- sum(final_data_new$V1 <= SNP_number_threshold)
S_less <- which((final_data_new$V1 <= SNP_number_threshold))
Genes_with_SNP_number_less_than_threshold <- sapply(S_less, function(x) gcpbayes_input_data_final[x])

S_greater_count <- sum(final_data_new$V1 > SNP_number_threshold)
S_greater <- which((final_data_new$V1 > SNP_number_threshold))
Genes_with_SNP_number_greater_than_threshold <- sapply(S_greater, function(x) gcpbayes_input_data_final[x])


# Loading Matrix data including Extra info data
load(file=paste0(path_input_data_read, Matrix_extra_info_input_file[1]))

# reading data from the loaded Matrix and preparing it for the splitting (next step)
len_Matrix_extra_info <- seq(1, length(gcpbayes_extra_input_data_final))
final_data_new_extra_info <- as.data.frame(do.call(rbind, lapply(len_Matrix_extra_info, function(x) gcpbayes_extra_input_data_final[[x]]$snp_number)))
summary(unlist(final_data_new_extra_info))

# splitting "Matrix_extra_info" data based on their number of SNPs using a threshold (less/greater than threshold)
S_less_count_extra_info <- sum(final_data_new_extra_info$V1 <= SNP_number_threshold)
S_less_extra_info <- which((final_data_new_extra_info$V1 <= SNP_number_threshold))
Genes_with_SNP_number_less_than_threshold_extra_info <- sapply(S_less_extra_info, function(x) gcpbayes_extra_input_data_final[x])

S_greater_count_extra_info <- sum(final_data_new_extra_info$V1 > SNP_number_threshold)
S_greater_extra_info <- which((final_data_new_extra_info$V1 > SNP_number_threshold))
Genes_with_SNP_number_greater_than_threshold_extra_info <- sapply(S_greater_extra_info, function(x) gcpbayes_extra_input_data_final[x])


#------------------------------------------------
ngenes <- length(Genes_with_SNP_number_less_than_threshold_extra_info)
DS_gene_bglm_gaussian_table <- NULL
#------------------------------------------------

# This index uses for printing output of genes with errors
j <- 1

# This index uses for printing output of pleiotropic genes
jj <- 1

# Start running GCPBayes for each Gene (one by one)
for(i in 1:ngenes){
  tic()
  
  mg <- length(Genes_with_SNP_number_less_than_threshold[[i]]$Betah[[1]])
  genename <- Genes_with_SNP_number_less_than_threshold_extra_info[[i]]$gene[[1]]
  snpnames <- names(Genes_with_SNP_number_less_than_threshold[[i]]$Betah[[1]])
  list_betah <- Genes_with_SNP_number_less_than_threshold[[i]]$Betah
  list_sigmah <- Genes_with_SNP_number_less_than_threshold[[i]]$Sigmah
  
  # Launch of GCPBayes function (here DS)
  Resume_test = try(DS(Betah = list_betah, Sigmah = list_sigmah, 
                       kappa0 = 0.5, 
                       sigma20 = 1,
                       m = mg,
                       K = 2,
                       niter = 3000, 
                       burnin = 1000, 
                       nthin = 2,
                       nchains = 1, 
                       a1 = 0.1, 
                       a2 = 0.1, 
                       d1= 0.1, 
                       d2= 0.1,
                       snpnames = snpnames, 
                       genename = genename)
  )
  
  # skipping the gene if GCPBayes calculation got "error" (and saving the index of the gene in the "error" output file)
  if ( class(Resume_test) == "try-error") {
    print("error?")
    if (j==1) {
      write.table(i, file=paste0(path_data_write,output_names[1], ".txt"), 
                  quote = F, row.names = F, col.names = F)
      j=j+1
    } else {
      write.table(i, file=paste0(path_data_write,output_names[1], ".txt"), 
                  append = T, quote = F, row.names = F, col.names = F)
    }
    next
  }
  
  # storing and saving results
  # also saving runtime for each Gene pleiotropy calculation
  exectime <- toc()
  exectime <- exectime$toc - exectime$tic
  res_line <- c(i,
             Resume_test$Criteria$`Name of Gene`,
             Resume_test$Criteria$log10BF,
             Resume_test$Criteria$lBFDR,
             Resume_test$Criteria$theta,
             Resume_test$Criteria$PPA[[1]],
             Resume_test$Criteria$PPA[[2]],
             genename,
             mg,
             exectime,
             Genes_with_SNP_number_less_than_threshold_extra_info[[i]][1]$gene,
             Genes_with_SNP_number_less_than_threshold_extra_info[[i]][2]$snp_number,
             Genes_with_SNP_number_less_than_threshold_extra_info[[i]][3]$chr,
             Genes_with_SNP_number_less_than_threshold_extra_info[[i]][5]$start,
             Genes_with_SNP_number_less_than_threshold_extra_info[[i]][6]$end,
             Genes_with_SNP_number_less_than_threshold_extra_info[[i]][7]$gene_type,
             Genes_with_SNP_number_less_than_threshold_extra_info[[i]][8]$gene_length,
             Genes_with_SNP_number_less_than_threshold_extra_info[[i]][9]$snp_number_gene_length_ratio
             )
  res_line_df <- t(as.data.frame(res_line))
  colnames(res_line_df) <- c("gene_index", "gene","log10BF","lBFDR","theta","PPA1","PPA2","gene","SNP_numbers","run_time",
                          "gene", "snp_number", "chr", "start", "end", "gene_type", "gene_length", "snp_number_gene_length_ratio")

  # saving results for the gene
  if (i==1) {
    write.table(res_line_df, file=paste0(path_data_write,output_names[2],".txt"),
                col.names = T, row.names = F, quote = F)
  } else {
    write.table(res_line_df, file=paste0(path_data_write,output_names[2],".txt"),
                append=T, col.names = F, row.names = F, quote = F)
  }
  
  # saving results if the gene has potential pleiotropic effects (based on its theta value)
  if (Resume_test$Criteria$theta > theta_threshold) {
    if (jj==1) {
      Plei_col_names <- c("GeneIndex", "GeneName", "#SNPs", "Theta","Pleiotropic_Effect_based_on_CI", "Pleiotropic_Effect_based_on_Median",
                          "gene", "snp_number", "chr", "start", "end", "gene_type", "gene_length", "snp_number_gene_length_ratio")
      
      write.table(t(Plei_col_names), file=paste0(path_data_write,output_names[3],".txt"), col.names = F, row.names = F, quote = F)
      
      if (sum(Resume_test$Indicator$`Significant studies and Pleiotropic effect based on CI`$`Pleiotropic effect`=="Yes") != 0 ) {
        Pleiotropic_effect_based_on_CI <- "Yes"
      }      else {
        Pleiotropic_effect_based_on_CI <- "No"
      }
      if (sum(Resume_test$Indicator$`Significant studies and Pleiotropic effect based on median thresholding`$`Pleiotropic effect`=="Yes") != 0 ) {
        Pleiotropic_effect_based_on_median <- "Yes"
      }      else {
        Pleiotropic_effect_based_on_median <- "No"
      }
      
      Gene_plei <- data.frame(c(i, genename, mg, Resume_test$Criteria$theta, Pleiotropic_effect_based_on_CI, Pleiotropic_effect_based_on_median,
                                Genes_with_SNP_number_less_than_threshold_extra_info[[i]][1]$gene,
                                Genes_with_SNP_number_less_than_threshold_extra_info[[i]][2]$snp_number,
                                Genes_with_SNP_number_less_than_threshold_extra_info[[i]][3]$chr,
                                Genes_with_SNP_number_less_than_threshold_extra_info[[i]][5]$start,
                                Genes_with_SNP_number_less_than_threshold_extra_info[[i]][6]$end,
                                Genes_with_SNP_number_less_than_threshold_extra_info[[i]][7]$gene_type,
                                Genes_with_SNP_number_less_than_threshold_extra_info[[i]][8]$gene_length,
                                Genes_with_SNP_number_less_than_threshold_extra_info[[i]][9]$snp_number_gene_length_ratio))
      
      write.table(t(Gene_plei), file=paste0(path_data_write,output_names[3], ".txt"),
                  append=T, col.names = F, row.names = F, quote = F)
      jj=jj+1
    } else { 
      if (sum(Resume_test$Indicator$`Significant studies and Pleiotropic effect based on CI`$`Pleiotropic effect`=="Yes") != 0 ) {
        Pleiotropic_effect_based_on_CI <- "Yes"
      }      else {
        Pleiotropic_effect_based_on_CI <- "No"
      }
      if (sum(Resume_test$Indicator$`Significant studies and Pleiotropic effect based on median thresholding`$`Pleiotropic effect`=="Yes") != 0 ) {
        Pleiotropic_effect_based_on_median <- "Yes"
      }      else {
        Pleiotropic_effect_based_on_median <- "No"
      }
      
      Gene_plei <- data.frame(c(i, genename, mg, Resume_test$Criteria$theta, Pleiotropic_effect_based_on_CI, Pleiotropic_effect_based_on_median,
                                Genes_with_SNP_number_less_than_threshold_extra_info[[i]][1]$gene,
                                Genes_with_SNP_number_less_than_threshold_extra_info[[i]][2]$snp_number,
                                Genes_with_SNP_number_less_than_threshold_extra_info[[i]][3]$chr,
                                Genes_with_SNP_number_less_than_threshold_extra_info[[i]][5]$start,
                                Genes_with_SNP_number_less_than_threshold_extra_info[[i]][6]$end,
                                Genes_with_SNP_number_less_than_threshold_extra_info[[i]][7]$gene_type,
                                Genes_with_SNP_number_less_than_threshold_extra_info[[i]][8]$gene_length,
                                Genes_with_SNP_number_less_than_threshold_extra_info[[i]][9]$snp_number_gene_length_ratio))
      
      write.table(t(Gene_plei), file=paste0(path_data_write,output_names[3], ".txt"),
                  append=T, col.names = F, row.names = F, quote = F)
      
    }
  }
  
  # printing the index of the gene in the command line output
  print(i)
}

# ================================================================================
# write the details information of the computer and softwares related to this code 
# please keep this info since it might be needed for any scientist in the future!
sessionInfo()
# calculation the running time of the program
end_time <- Sys.time()
end_time - start_time
# ================================================================================

#!/usr/bin/env Rscript

# 		                Script for running PLACO 
# ================================================================================
# Summary: 
# running PLACO based on GWAS summary statistic files which contains common SNPs
# ================================================================================
# Written first by: Elise Lucotte
# Modified by: Yazdan Asgari
# Initial Creation Date: 12/2020
# Edited Date: 12/2022
# ================================================================================

# for calculation the running time of the program
start_time <- Sys.time()

# ================================================================================
# libraries used in the code
# please install the following packages if needed
require(devtools)
source_url("https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.1.1.R?raw=TRUE")
require(MASS)
library(dplyr)
library(tidyr)
library(vroom)
library(optparse)

# ================================================================================
# Read of command inputs
# ================================================================================

option_list = list(
  make_option(c("--path"), type="character", default=NULL, 
              help="directory path", metavar="character"),
  make_option(c("--inputs"), type="character", default=NULL, 
              help="name of file to read for inputs", metavar="character"),
  make_option(c("--prefix_in"), type="character", default="step1_out1", 
              help="prefix file name for inputs [default= %default]", metavar="character"),
  make_option(c("--prefix_out"), type="character", default="step1_out1", 
              help="prefix file name for output [default= %default]", metavar="character"),
  make_option(c("--pval"), type="numeric", default=NULL, 
              help="threshold used for decorrelating the Z-scores", metavar="number")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# ================================================================================
# 						DEFINITION SECTION
# ================================================================================

# working directory
path = opt$path

# directory in which input data exist
inputs <- read.table(file=paste0(opt$path,opt$inputs),header=TRUE)
nstud <- dim(inputs)[1]

# the file names of input data (which must include "snp", "Z", and "pval" columns)
traits_1 <- paste(opt$prefix_in, inputs$shortname[1], sep="_")
traits_2 <- paste(opt$prefix_in, inputs$shortname[2], sep="_")

# the file names of the output data 
output_name <- opt$prefix_out

# used for decorrelating the Z-scores
pval_threshold <-opt$pval



# ================================================================================
# 						            RUNNING SECTIONS
# ================================================================================
# 						               ### 1 ###
# Summary:
# Reading input data, checking SNPs order in both traits, creation of Z and P matrices
# performing a total correlation between two traits to check whether a "decorrelation step" is needed or not
# ================================================================================
for (i in 1:length(traits_1)) {
  
  # reading first trait (very fast)
  t1 <- vroom(file=paste0(path, traits_1[i],'.txt'))
  t1 <- as.data.frame(t1)
  
  # reading second trait (very fast)
  t2 <- vroom(file=paste0(path, traits_2[i],'.txt'))
  t2 <- as.data.frame(t2)
  
  # Checking if all the SNPs are in the same order
  writeLines("\n\n")
  print('Check if all the SNPs are in the same order')
  print(all(t1$snp==t2$snp))
  
  # creation of Z and P matrices
  Z.matrix <- cbind(t1$Z, t2$Z)
  P.matrix <- cbind(t1$pval, t2$pval)
  k <- 2
  colnames(Z.matrix) <- paste("Z", 1:k, sep="")
  colnames(P.matrix) <- paste("P", 1:k, sep="")
  
  writeLines("\n\n")
  print('End of creation of Z and P matrices')
  
  # performing a total correlation between two traits
  print(cor.test(t1$Z, t2$Z, method='pearson'))
  A <- cor.test(t1$Z, t2$Z, method='pearson')
  
  writeLines("\n\n")
  print('End of performing a total correlation between two traits')

  # deciding whether a "Decorrelate step" is needed or not  
  # based on comparing the obtained p_value between two traits and the "pval_threshold"
  if (A$p.value < pval_threshold) {
    writeLines("\n\n")
    print('It needs decorrelating the Z-scores')
    # Decorrelate the matrix of Z-scores
    R <- cor.pearson(Z.matrix, P.matrix, p.threshold=1e-4)
    # function for raising matrix to any power
    "%^%" <- function(x, pow)
      with(eigen(x), vectors %*% (values^pow * t(vectors)))
    Z.matrix.decor <- Z.matrix %*% (R %^% (-0.5))
    colnames(Z.matrix.decor) <- paste("ZD", 1:k, sep="")
  } else {
    writeLines("\n\n")
    print('No need to decorrelate the Z-scores')
    Z.matrix.decor <- Z.matrix
  }
  
  writeLines("\n\n")
  print('End of Running Section:#1')

  # ================================================================================
  # 						               ### 2 ###
  # Summary: 
  # Obtaining the variance parameter estimates
  # ================================================================================  
  VarZ <- var.placo(Z.matrix.decor, P.matrix, p.threshold=1e-4)
  writeLines("\n\n")
  print('End of Running Section:#2')

  # ================================================================================
  # 						               ### 3 ###
  # Summary: 
  # Applying test of pleiotropy for each variant
  # ================================================================================  
  out_1 <- sapply(1:nrow(t1), function(i) placo(Z=Z.matrix.decor[i,], VarZ=VarZ))
  out_2 <- t(out_1)
  out_2 <- as.data.frame(out_2)
  out_2$T.placo <- as.numeric(out_2$T.placo)
  out_2$p.placo <- as.numeric(out_2$p.placo)
  
  # calculation of Adjusted p_value for each variant
  q.value <- p.adjust(out_2$p.placo, method='fdr', n=nrow(out_2))
  out_2 <- mutate(out_2, FDR.placo=q.value)
  t3 <- merge(t1,t2, by=c('snp','chr','bp_hg19'), sort=FALSE)
  all_t <- cbind(t3, out_2)
  
  writeLines("\n\n")
  print('End of Running Section:#3')

  # writing the output files (very fast)
  vroom_write(all_t, file = paste0(path, output_name[i],'.txt'))
}

# ================================================================================
# write the details information of the computer and softwares related to this code 
# please keep this info since it might be needed for any scientist in the future!
sessionInfo()
# calculation the running time of the program
end_time <- Sys.time()
end_time - start_time
# ================================================================================

# ================================================================================
# 		                Script for running ld_clumping locally
# ================================================================================
# Summary:
# In this program, ld_clumping is run locally based on recommendation in "ieugwasr" library manual
# NOTE1:
# The script assumes there is a "plink" which its PATH could be recognized by "genetics.binaRies::get_plink_binary()" command
# If the above command returns an empty path, please install "plink" using the following command:
# devtools::install_github("explodecomputer/genetics.binaRies")
# NOTE2:
# The script needs reference bfiles (bed/bim/fam) to be located in the path that you defined in the script
# ================================================================================
# Written first by: Yazdan Asgari
# Initial Creation Date: 10/2021
# Modified by: Yazdan Asgari
# Edited Date: 03/2023
# ================================================================================

# for calculation the running time of the program
start_time <- Sys.time()

# ================================================================================
# libraries used in the code
library(vroom)

# ================================================================================
# 						PATH DEFINITION SECTION
# ================================================================================
# directory path for the input data
path <- path_input_data

# input file name
placo_input_file <- output_step2

# output file name
output_file <- output_step3

# ================================================================================
# Definition of functions needed for running ld_clumping locally
# ================================================================================
ld_clump_local <- function(dat, clump_kb, clump_r2, clump_p, bfile, plink_bin)
{
  # Make textfile
  shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
  fn <- tempfile()
  write.table(data.frame(SNP=dat[["rsid"]], P=dat[["pval"]]), file=fn, row.names=F, col.names=T, quote=F)

  fun2 <- paste0(
    shQuote(plink_bin, type=shell),
    " --bfile ", shQuote(bfile, type=shell),
    " --clump ", shQuote(fn, type=shell),
    " --clump-p1 ", clump_p,
    " --clump-r2 ", clump_r2,
    " --clump-kb ", clump_kb,
    " --allow-extra-chr ",
    " --out ", shQuote(fn, type=shell)
  )

  system(fun2)
  res <- read.table(paste(fn, ".clumped", sep=""), header=T)
  unlink(paste(fn, "*", sep=""))
  y <- subset(dat, !dat[["rsid"]] %in% res[["SNP"]])
  if(nrow(y) > 0)
  {
    message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), " variants due to LD with other variants or absence from LD reference panel")
  }
  return(subset(dat, dat[["rsid"]] %in% res[["SNP"]]))
}
# ================================================================================
random_string <- function(n=1, len=6)
{
  randomString <- c(1:n)
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    len, replace=TRUE),
                             collapse="")
  }
  return(randomString)
}
# ================================================================================
ld_clump <- function(dat=NULL, clump_kb=10000, clump_r2=0.001, clump_p=0.99, pop = "EUR", access_token=NULL, bfile=NULL, plink_bin=NULL)
{

  stopifnot("rsid" %in% names(dat))
  stopifnot(is.data.frame(dat))

  if(is.null(bfile))
  {
    message("Please look at vignettes for options on running this locally if you need to run many instances of this command.")
  }

  if(! "pval" %in% names(dat))
  {
    if( "p" %in% names(dat))
    {
      warning("No 'pval' column found in dat object. Using 'p' column.")
      dat[["pval"]] <- dat[["p"]]
    } else {
      warning("No 'pval' column found in dat object. Setting p-values for all SNPs to clump_p parameter.")
      dat[["pval"]] <- clump_p
    }
  }

  if(! "id" %in% names(dat))
  {
    dat$id <- random_string(1)
  }

  if(is.null(bfile))
  {
    access_token = check_access_token()
  }

  ids <- unique(dat[["id"]])
  res <- list()
  for(i in 1:length(ids))
  {
    x <- subset(dat, dat[["id"]] == ids[i])
    if(nrow(x) == 1)
    {
      message("Only one SNP for ", ids[i])
      res[[i]] <- x
    } else {
      message("Clumping ", ids[i], ", ", nrow(x), " variants, using ", pop, " population reference")
      if(is.null(bfile))
      {
        res[[i]] <- ld_clump_api(x, clump_kb=clump_kb, clump_r2=clump_r2, clump_p=clump_p, pop=pop, access_token=access_token)
      } else {
        res[[i]] <- ld_clump_local(x, clump_kb=clump_kb, clump_r2=clump_r2, clump_p=clump_p, bfile=bfile, plink_bin=plink_bin)
      }
    }
  }
  res <- dplyr::bind_rows(res)
  return(res)
}

# ================================================================================
# 						RUNNING SECTION
# ================================================================================
# Reading PLACO output file
# "rsid" and "pval" columns should be available
# "pval" is the P-value calculated by PLACO

SNPs_communs_gene_Z_P_Placo <- vroom(file=paste0(path, placo_input_file, ".txt"))
SNPs_communs_gene_Z_P_Placo <- as.data.frame(SNPs_communs_gene_Z_P_Placo)

# renaming the columns names
colnames(SNPs_communs_gene_Z_P_Placo) <- c("rsid","chr","pos_hg19",
                                           "Effect_A_bcac","nonEffect_A_bcac","beta_bcac","se_bcac","p_bcac","info_bcac","EAF_bcac","MAF_bcac","Z_bcac",
                                           "Effect_A_OCAC","nonEffect_A_OCAC","beta_OCAC","se_OCAC","p_OCAC","info_OCAC","EAF_OCAC","MAF_OCAC","Z_OCAC",
                                           "T_placo","pval","FDR_placo")

# extraction of the selected columns
SNPs_communs_gene_Z_P_Placo <- SNPs_communs_gene_Z_P_Placo[, c("rsid","chr","pos_hg19",
                                                               "Effect_A_bcac","nonEffect_A_bcac","pval")]
# renaming the columns names
colnames(SNPs_communs_gene_Z_P_Placo) <- c("rsid","chr","pos",
                                           "Effect_A","nonEffect_A","pval")
SNPs_communs_gene_Z_P_Placo$pval[which(SNPs_communs_gene_Z_P_Placo$pval>1)]<-1

# splitting data based on chromosome numbers
SNPs_communs_gene_Z_P_Placo_split <- split(SNPs_communs_gene_Z_P_Placo, SNPs_communs_gene_Z_P_Placo[,"chr"])

# extraction of the chromosome names
chromosomes <- names(SNPs_communs_gene_Z_P_Placo_split)

# LD Calculation using PLINK clumping method for chromosomes 1-22
for (i in 1:22) {
  cat("\n")
  cat("Chromosome",i)
  cat("\n")
  chr_i = chromosomes[i]

  result_ld <- ld_clump(
    dat = SNPs_communs_gene_Z_P_Placo_split[[chr_i]],
    clump_kb = clump_threshold_kb,
    clump_r2 = clump_threshold_r2,
    clump_p = clump_threshold_p,
    pop = "EUR",
    access_token = NULL,
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = paste0(ref_path_b_files, ref_b_files)
  )

  if (i==1){
    write.table(result_ld, file = paste0(path, output_file, ".txt"),
                col.names = F,  row.names = F, quote = F)
  }
  else {
    write.table(result_ld, file = paste0(path, output_file, ".txt"),
                append = T, col.names = F,  row.names = F, quote = F)
  }

}

# ================================================================================
# write the details information of the computer and softwares related to this code
# please keep this info since it might be needed for any scientist in the future!
sessionInfo()
# calculation the running time of the program
end_time <- Sys.time()
end_time - start_time
# ================================================================================


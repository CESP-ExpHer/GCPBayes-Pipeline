# 		              Script for checking of BCAC GWAS input file
# ================================================================================
# Summary: 
# it extracts some SNPs from the BCAC GWAS summary statistics file and compare them with GWAS database (such as GWAS catalog and 1K Genome)  
# ================================================================================
# Written by: CheckSumStats developers (https://github.com/MRCIEU/CheckSumStats)
# modified Date: 1/2022
# https://cesp.inserm.fr/en/equipe/exposome-and-heredity
# ================================================================================
# for calculation the running time of the program
start_time <- Sys.time()

# ================================================================================
# 						DEFINITION SECTION
# should be changed by a user
# ================================================================================
path_input_data = "~/BCAC_OCAC/"
input_data = "BCAC_2020_onco_ALL_reformatted.txt"

# ================================================================================
# libraries used in the code
# please install the following packages if needed
# BiocManager::install("ieugwasr")
# BiocManager::install("gwasrapidd")
# BiocManager::install("devtools")
# devtools::install_github("MRCIEU/CheckSumStats")
library(ieugwasr)
library(gwasrapidd)
library(CheckSumStats)

# ================================================================================
# 						RUNNING SECTIONS
# ================================================================================
# Extract SNPs list based on EFO ID of Breast Carcinoma
EFO <- get_efo(trait="Breast Carcinoma")
EFO
snplist <- make_snplist(efo_id = EFO$efo_id, trait = "Breast Carcinoma", ref1000G_superpops = TRUE)

# ========================================================================================
# reading BCAC GWAS data and extracting just SNPs that created from the previous step ("snplist")
bc <- extract_snps(snplist = snplist, path_to_target_file = paste0(path_input_data, input_data))

# 1
# checking allele frequency
Dat <- format_data(dat = bc, outcome = "Breast Carcinoma", population = "European", rsid = "snp", effect_allele = "Effect_A", other_allele = "nonEffect_A", beta = "beta", se = "se", eaf = "EAF", p = "pval", efo = "Breast Carcinoma")
Plot1 <- make_plot_maf(ref_1000G = c("AFR", "AMR", "EAS", "EUR", "SAS", "ALL"), target_dat = Dat)
Plot1

# checking allele frequency output IF we select the columns WRONG
Dat_wrong <- format_data(dat = bc, outcome = "Breast Carcinoma", population = "European", rsid = "snp", effect_allele = "nonEffect_A", other_allele = "Effect_A", beta = "beta", se = "se", eaf = "EAF", p = "pval", efo = "Breast Carcinoma")
Plot1_wrong <- make_plot_maf(ref_1000G = c("AFR", "AMR", "EAS", "EUR", "SAS", "ALL"), target_dat = Dat_wrong)
Plot1_wrong

# ========================================================================================
# 2
# checking the effect allele by comparing Z scores in the test and GWAS catalog datasets
# it needs "lnor" instead of "beta" as column name
Dat1 <- format_data(dat = bc, outcome = "Breast Carcinoma", population = "European", rsid = "snp", effect_allele = "Effect_A", other_allele = "nonEffect_A", lnor = "beta", se = "se", eaf = "EAF", p = "pval", efo = "Breast Carcinoma")
Plot2 <- make_plot_gwas_catalog(dat = Dat1, efo_id = EFO$efo_id, trait = "Breast Carcinoma")
Plot2

# ========================================================================================
# 3
# checking the effect allele by comparing effect allele frequency (EAF) between the test dataset and the GWAS catalog
Plot3 <- make_plot_gwas_catalog(dat = Dat1, plot_type = "plot_eaf", efo = unique(Dat1$efo), trait = unique(Dat1$outcome))
Plot3

# ========================================================================================
# 4
# checking for errors or analytical issues in the summary data
bc_4 <- extract_snps(snplist = snplist, path_to_target_file = paste0(path_input_data, input_data), get_sig_snps = TRUE, p_val_col_number = 8)
# The "predict_lnor_sh" function needs number of "cases" and "controls" samples
# So, here we added these two columns to the input data
ncases <- rep(72000, length(bc_4$snp))
ncontrols <- rep(58383, length(bc_4$snp))
bc_4$ncases <- ncases
bc_4$ncontrols <- ncontrols
Dat_4 <- format_data(dat = bc_4, outcome = "Breast Carcinoma", population = "European", rsid = "snp", effect_allele = "Effect_A", other_allele = "nonEffect_A", lnor = "beta", se = "se", eaf = "EAF", p = "pval", efo = "Breast Carcinoma")

# if #SNPs > 100, it is recommended to run ld_clumping to shrink the data, because "predict_lnor_sh" function is very slow to run
Clump <- ieugwasr::ld_clump(clump_r2 = 0.01, clump_p = 1e-8, dplyr::tibble(rsid = Dat_4$rsid, pval = Dat_4$p, id = Dat_4$id), pop = "EUR")
Dat_5 <- Dat_4[Dat_4$rsid %in% Clump$rsid,]
Pred_clump <- predict_lnor_sh(dat=Dat_5)

# comparing the expected and reported effect sizes.
Plot4 <- make_plot_pred_effect(dat = Pred_clump)

# plotting the relative bias, i.e. the percentage deviation of the expected from the reported effect size
Plot4_2 <- make_plot_pred_effect(dat = Pred_clump, bias = TRUE)
Plot4_2

# ========================================================================================
# 5
# checking whether the reported P values correspond to the reported effect sizes in the dataset
# it needs "lnor" instead of "beta" as name
Dat5 <- format_data(dat = bc, outcome = "Breast Carcinoma", population = "European", rsid = "snp", effect_allele = "Effect_A", other_allele = "nonEffect_A", lnor = "beta", se = "se", eaf = "EAF", p = "pval", efo = "Breast Carcinoma")
Plot5 <- zz_plot(dat = Dat5)

# ========================================================================================
# Saving all outputs in one file
setwd("~/BCAC2020_OCAC/")
Plot_list2 <- list("Plot1", "Plot1_wrong", "Plot2", "Plot3", "Plot4", "Plot4_2", "Plot5")
Plot_list <- lapply(1:length(Plot_list2), FUN = function(x) eval(parse(text = Plot_list2[x])))
combine_plots(Plot_list = Plot_list, out_file = "checksumstats_BCAC_2020_report.png")

# ================================================================================
# write the details information of the computer and softwares related to this code 
# please keep this info since it might be needed for any scientist in the future!
sessionInfo()
# calculation the running time of the program
end_time <- Sys.time()
end_time - start_time
# ================================================================================



# 		                Script for Analysis of PLACO Outputs
# ================================================================================
# Summary: 
# Analysis of PLACO outputs and creation of a Manhattan Plot for chr 1-22
# And showing Significant SNPs (positive and negative) in the plots
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
# libraries used in the code
# please install the following packages if needed
#BiocManager::install("ggplot2")
#BiocManager::install("patchwork")
#BiocManager::install("tidyr")
#BiocManager::install("dplyr")
#BiocManager::install("vroom")
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(vroom)

# ================================================================================
# DEFINITION SECTION
# ================================================================================
# directory in which input data exist
path_input_data = "~/BCAC_OCAC/"

# the file name of input data
placo_result_name <- "output_PLACO_BCAC_2020_ALL_OCAC.txt"

# directory in which an output data would be written
path_output_data = "~/BCAC_OCAC/"

# the file names of the output data
output_name <- "output_PLACO_BCAC_2020_ALL_OCAC_Sig.txt"

# ================================================================================
# RUNNING SECTION
# ================================================================================
# reading the input data
placo_result <- vroom(file=paste0(path_input_data, placo_result_name))
placo_result <- as.data.frame(placo_result)

# performing a correlation analysis between two traits (using Z columns) 
# correlation test is based on two methods: spearman and pearson
A <- cor.test(placo_result$Z.x, placo_result$Z.y, method='spearman')
A
B <- cor.test(placo_result$Z.x, placo_result$Z.y, method='pearson')
B

# writing correlation test output into a text file
capture.output(A, file = paste0(path_output_data, substr(output_name, 1, nchar(output_name)-8), "_cor_test.txt") )
capture.output(B, file = paste0(path_output_data, substr(output_name, 1, nchar(output_name)-8), "_cor_test.txt"), append = TRUE )

# finding significant SNPs (p < 5*10^-8) and divide them into "positive" and "negative" effect groups
ds <- placo_result[which(placo_result$p.placo<5*10^-8),]
dpos <- ds[which(ds$T.placo>0),]
dneg <- ds[which(ds$T.placo<0),]

# writing output data
vroom_write(ds, file = paste0(path_output_data, output_name))

# creation of a Manhattan Plot
dpos <- mutate(dpos, pos2=round(bp_hg19/1000000, digits=2))
dneg <- mutate(dneg, pos2=round(bp_hg19/1000000, digits=2))

dall <- mutate(placo_result, pos2=round(bp_hg19/1000000, digits=2))

both <- ggplot() +
  geom_point(data=dall, aes(x=pos2, y=-log10(p.placo)), col="grey") +
  geom_point(data=dpos, aes(x=pos2, y=-log10(p.placo)), col="cornflowerblue") +
  geom_point(data=dneg, aes(x=pos2, y=-log10(p.placo)), col="brown2") +
  theme_bw() +
  ggtitle(substr(output_name, 1, nchar(output_name)-8)) +
  geom_hline(yintercept=-log10(5E-8), linetype='dashed', color='black') +
  labs(x ="Position (Mb)") +
  theme(axis.text.x = element_text(angle=45)) +
  facet_wrap(~as.numeric(chr), scales='free')

ggsave(both, filename=paste0(path_output_data, substr(output_name, 1, nchar(output_name)-8), 
                             '_Manhattan_Plot.png', sep=''), device='png', w=8, h=8 )

# ================================================================================
# write the details information of the computer and softwares related to this code 
# please keep this info since it might be needed for any scientist in the future!
sessionInfo()
# calculation the running time of the program
end_time <- Sys.time()
end_time - start_time
# ================================================================================


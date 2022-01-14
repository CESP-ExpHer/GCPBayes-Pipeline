# 		              Script for ploting a Pie chart of Annotation file
# ================================================================================
# Summary: 
# creation of a pie-chart for an overview of genes (here just coding-genes) distribution among all chromosomes 
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
# should be changed by a user
# ================================================================================
path_input_data = "~/BCAC_OCAC/"
input_data = "annot_gencode_v38lift37_modified_gene_class_coding.txt"

# ================================================================================
# libraries used in the code
# please install the following packages if needed
#BiocManager::install("ggplot2")
#BiocManager::install("plotly")
library(ggplot2)
library(plotly)

# ================================================================================
# 						RUNNING SECTIONS
# ================================================================================
annot_coding_final <- read.table(file = paste0(path_input_data, input_data), header = TRUE)
  
summary(annot_coding_final)
summary(annot_coding_final$seqnames)
sort(summary(annot_coding_final$seqnames),decreasing = TRUE)

# creation of pie chart for number of protein-coding genes in each chromosome
A <- summary(annot_coding_final$seqnames)
B <- paste0('chr', 1:22)
B <- c(B,'chrX','chrY','chrM')
gene_df <- data.frame(A,B)


C <- plot_ly(gene_df, labels = ~B, values = ~A, type = 'pie',textposition = 'outside',textinfo = 'label+percent') %>%
  layout(
    #title = 'Protein_coding Genes Distributions (hg19)',
    xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
    yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

# ================================================================================
# Saving plot using 
# https://chart-studio.plotly.com/
# https://plotly.com/r/chart-studio-image-export/

# Let the R session know about your Chart Studio authorization credentials by setting environment variables using Sys.setenv()
Sys.setenv("plotly_username"="yazdan59")
Sys.setenv("plotly_api_key" = "1egJuyytCqJ9Qi0uftMU")

Jpeg <- plotly_IMAGE(C, format = "jpeg", out_file = paste0(path_input_data, "Protein-Coding_Genes_Distribution_Pie_Chart.jpeg"))
Pdf <- plotly_IMAGE(C, format = "pdf",  out_file = paste0(path_input_data, "Protein-Coding_Genes_Distribution_Pie_Chart.pdf"))
Svg <- plotly_IMAGE(C, format = "svg",  out_file = paste0(path_input_data, "Protein-Coding_Genes_Distribution_Pie_Chart.svg"))

# # static
# orca(C, file=paste0(path_input_data, "Protein-Coding_Genes_Distribution_Pie_Chart.pdf"))
# 
# # in jpeg format
# jpeg(file=paste0(path_input_data, "Protein-Coding_Genes_Distribution_Pie_Chart.jpeg"))
# C
# dev.off()
# 
# # in pdf format
# pdf(file=paste0(path_input_data, "Protein-Coding_Genes_Distribution_Pie_Chart.pdf"))
# C
# dev.off()

# ================================================================================
# write the details information of the computer and softwares related to this code 
# please keep this info since it might be needed for any scientist in the future!
sessionInfo()
# calculation the running time of the program
end_time <- Sys.time()
end_time - start_time
# ================================================================================

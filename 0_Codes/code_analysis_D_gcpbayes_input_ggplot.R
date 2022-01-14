# ================================================================================
# Summary:
# Extraction some information about GCPBayes input data such as:
# Total Gene Numbers, SNPs Distribution among the Genes
# ================================================================================
# Initially Written by: Yazdan Asgari
# Edited Date: 1/2022
# https://cesp.inserm.fr/en/equipe/exposome-and-heredity
# ================================================================================
# 						DEFINITION SECTION
# should be changed by a user
# ================================================================================
# PATH in which a GCPBayes input file exists
gcpbayes_input_path <- "~/BCAC_OCAC/"

# GCPBayes input file name
gcpbayes_input <- "Matrices_output_pipeline_BCAC_ALL_OCAC_coding_clumping_08.Rdata"

# A threshold for dividing the results into two plots in one page  
threshold <- 50
# =================================================================================

# Loading the GCPBayes input data 
load(file=paste0(gcpbayes_input_path, gcpbayes_input))

library('ggplot2')
library('patchwork')

# Histogram ALL Genes
N <- sapply( seq(1:length(gcpbayes_input_data_final)), function(x) length(gcpbayes_input_data_final[[x]]$Betah[[1]]) )
N1 <- as.data.frame(N)
length(N1$N)
breaks_N <- c(1,2,3,4,5,6,7,8,9,10,100,500,1000,max((N1$N)))

p0 <- ggplot (N1, aes(N)) + 
  geom_histogram(colour="darkgreen", size=1, fill="green", bins = 17 ) + 
  stat_bin(geom="text", aes(label=..count..) , hjust=-0.4, vjust = 0.4, bins = 17, angle = 90, size = 2.5) + 
  scale_x_log10('Numbers of SNPs', breaks = breaks_N, labels = breaks_N, expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous('Freq', expand = expansion(mult = c(0, 0.1))) +
  ggtitle("SNP Distribution _ Histogram ALL Genes") +
  theme(axis.text.x = element_text(size = 8, angle = 60, vjust = 0.8, hjust=0.5), plot.title = element_text(hjust = 0.5)) 

p0

#=====================================================================
# Two Plots in one page dividing by a threshold (Bar Plot + Histogram)

H <- sapply( seq(1:length(gcpbayes_input_data_final)), function(x) length(gcpbayes_input_data_final[[x]]$Betah[[1]]) )
class(H)
A <- as.data.frame(table(H))
colnames(A) <- c("snp_number", "snp_count")
class(A$snp_number)
A$snp_number<-as.numeric(as.character(A$snp_number))
class(A$snp_number)
class(A$snp_count)
A$snp_count<-as.numeric(as.character(A$snp_count))
class(A$snp_count)



p1 <- ggplot( subset(A, A$snp_number < threshold), aes(x=snp_number, y=snp_count, fill = snp_number )) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = snp_count), size = 2, vjust = 0.5, hjust=-0.2, angle = 90, colour = "black") +
  theme(axis.text.x = element_text(size = 5, angle = 60, vjust = 0.8, hjust=1), legend.position="none") +
  scale_x_continuous('#SNPs less than the 50', expand = expansion(mult = c(0, 0.01)) ) +
  scale_y_continuous('Count', expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(size = 8, angle = 60, vjust = 0.8, hjust=1), plot.title = element_text(hjust = 0.5)) +
  ggtitle("SNP Distribution _ Bar Plot")

# p1

breaks_def <- floor((max(A$snp_number) - threshold)/15)
breaks <- seq(threshold, max(A$snp_number)-200, by = breaks_def)
breaks <- c(breaks,max(A$snp_number))
p2 <- ggplot( subset(A, A$snp_number > threshold), aes(x=snp_number) ) +
  geom_histogram(colour="darkgreen", size=1, fill="green") + 
  stat_bin(geom="text", aes(label=..count..), hjust=-0.8, vjust = 0.4, angle = 90, size = 2.5 ) + 
  scale_x_log10('#SNPs greater than the 50', breaks = breaks, labels = breaks, expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous('Count', expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(size = 8, angle = 60, vjust = 0.8, hjust=1), plot.title = element_text(hjust = 0.5)) +
  ggtitle("SNP Distribution _ Histogram")

# p2
p1 + p2

#===============================================================================

breaks_p0 <- c(100,200,300,400,500,600,700,800,900,1000, max(A$snp_number))
p0 <- ggplot( subset(A, A$snp_number > threshold), aes(x=snp_number) ) +
  geom_histogram(colour="darkgreen", size=1, fill="green", bins = 20, binwidth = 0.01) + 
  stat_bin(geom="text", aes(label=..count..) , hjust=-0.4, vjust = 0.4, bins = 20, binwidth = 0.1 , angle = 90, size = 2.5) + 
  scale_x_log10('#SNPs greater than the Threshold', breaks = breaks_p0, labels = breaks_p0, expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous('Count', expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(size = 8, angle = 60, vjust = 0.8, hjust=1))
p0
#=====================================================================
M <- sapply( seq(1:length(gcpbayes_input_data_final)), function(x) length(gcpbayes_input_data_final[[x]]$Betah[[1]]) )
M1 <- as.data.frame(M)
length(M1$M)

x <- unique(round(1.5^(seq(1:length(M1$M)))))

ggplot (M1, aes(M)) + 
  geom_histogram(colour="darkgreen", size=1, fill="green", binwidth = 0.1) + 
  stat_bin(geom="text", aes(label=..count..) , hjust=-0.4, vjust = 0.4, binwidth = 0.1 , angle = 90, size = 2.5) + 
  scale_x_log10('Integer Data', breaks = x, labels = x, expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous('Count', expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(size = 8, angle = 60, vjust = 0.8, hjust=1))
#=====================================================================




# GCPBayes Pipeline
Created by: Yazdan Asgari<br>
Creation date: 14 Jan 2022<br>
Update: Jul 2022<br>
https://cesp.inserm.fr/en/equipe/exposome-and-heredity
<br>
<br>

### Here we have provided detailed information about how to work with the *GCPBayes pipeline* including tips and suggestions. 
<br>

## NOTES
**IMPORTANT NOTE 1:** All SNPs and genes positions in the GWAS and annotation data are based on GRCh37 (hg19) Human Genome Assembly.
<br><br>
**IMPORTANT NOTE 2:** The Section names that are mentioned in this page are correspondence to Figure provided in [An Overview of the Pipeline](#An-Overview-of-the-Pipeline) section in the current page.
<br>
<br>
## Contents
This file includes **THREE** major sections:
- [**Description of Bash Files**](1) (plain text files that contain a series of commands) for running the whole procedure with a series of options.
**IMPORTANT NOTE:** The Bash file was tested on a Unix-based server with CentOS 7.

- [**A Step-by-Step Tutorial**](2) for how to use the GCPBayes pipeline to explore genes with potential pleiotropic effects on Breast and Ovarian cancers using GWAS summary statistics data of Breast (BCAC) and Ovarian (OCAC) cancers (which their results were presented in the manuscript).
**IMPORTANT NOTE:** This tutorial could be run on all Operating Systems (OS) including Windows, Linux, and Mac.

- [**GCPBayes Pipeline Wiki**](3) for All Scripts which includes detailed information for all scripts used in the pipeline.
**IMPORTANT NOTE:** This part includes all concepts and information for each scripts used in the GCPBayes pipeline in order to provide an option to modify/add/remove any part by a user based on its strategy.

## An Overview of the Pipeline
A general overview of the main sections of the GCPBayes pipeline are as follow:
<br></br>
<kbd> <img src="0_Images/Fig1_v3.jpg"/> </kbd>
<br></br>

## Test Dataset
Here, we provide a small dataset for testing the Pipeline. The data are GWAS summary statistics for The Breast Cancer Association Consortium (BCAC) and The Ovarian Cancer Association Consortium (OCAC) chromosome #5 and we want to run the Pipeline (without LD clumping method) and GCPBayes at a gene-level. For running the pipeline in the test set, please perform the following steps:
- Download input files: [Download](....)
  - BCAC and OCAC GWAS data on chromosome #5 (gwas_BCAC_chr5.txt, gwas_OCAC_chr5.txt)
  - An annotation file including all coding genes (annot_gencode_v38lift37_modified_gene_class.txt)
  - BCAC GWAS file with a gene column (Annot_BCAC_2020_onco_ALL_reformatted_coding.txt)
- Download the scripts and put them in the same folder as input data

## How to Cite
Asgari et al., "GCPBayes Pipeline: a tool for exploring pleiotropy at gene-level", xxxx. xxx x;x(x), doi:xxx [Link](https://..../)
<br>
<br>

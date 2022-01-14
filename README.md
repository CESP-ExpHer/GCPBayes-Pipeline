# GCPBayes Pipeline
Created by: Yazdan Asgari<br>
Creation date: 14 Jan 2022<br>
Update: Jan 2022<br>
https://cesp.inserm.fr/en/equipe/exposome-and-heredity
<br>
<br>
## How to Cite
Asgari et al., "GCPBayes Pipeline: a tool for exploring pleiotropy at gene-level", xxxx. xxx x;x(x):x, doi: xxx[Paper_link](https://..../)
<br>
<br>
## Some NOTES
Here we have provided detailed information about the GCPBayes pipeline including tips and suggestions for running various steps. 
<br>
**IMPORTANT NOTE 1:** All SNPs and genes positions in the GWAS and annotation data are based on GRCh37 (hg19) Human Genome Assembly.
<br>
**IMPORTANT NOTE 2:** The Section names that are mentioned in this page are correspondence to Figure provided in the [An Overview of the Pipeline](#An-Overview-of-the-Pipeline).
<br>
<br>
## Contents
This file includes **THREE** major sections:
- Description of a Bash file (a plain text file that contains a series of commands) for running the whole procedure with a series of options.
**IMPORTANT NOTE:** The Bash file was tested on a Unix-based server with CentOS 7.

- A step-by-step tutorial for how to use the GCPBayes pipeline to explore genes with potential pleiotropic effects on Breast and Ovarian cancers using GWAS summary statistics data of Breast (BCAC) and Ovarian (OCAC) cancers (which their results were presented in the manuscript).
**IMPORTANT NOTE:** This tutorial could be run on all Operating Systems (OS) including Windows, Linux, and Mac.

- Detailed information for all scripts in the pipeline.
**IMPORTANT NOTE:** This part includes all concepts and information for each scripts used in the GCPBayes pipeline in order to provide an option to modify/add/remove any part by a user based on its strategy.

## An Overview of the Pipeline
A general overview of the main sections of the GCPBayes pipeline are as follow:
<br></br>
<kbd> <img src="0_Images/Fig1_v3.jpg"/> </kbd>
<br></br>

## Acknowledgements 
We gratefully acknowledge the following packages which we used throughout our pipeline:
```
ggplot2
gwasrapidd
ieugwasr
biomaRt
dplyr
tibble
```

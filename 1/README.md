# Description of Bash file (a plain text file that contains a series of commands) for running the whole procedure with a series of options
Created by: Pierre-Emmanuel Sugier, Yazdan Asgari<br>
Creation date: 14 Jan 2022<br>
Update: Jan 2022<br>
https://cesp.inserm.fr/en/equipe/exposome-and-heredity
<br>
<br>


**IMPORTANT NOTE 1:** The Bash file was tested on a Unix-based server with CentOS 7.
<br>
<br>
**IMPORTANT NOTE 2:** We recommend a user to edit **"Global Bash file"** and use it for running any GWAS summary statistics data two traits. 
<br>
<br>
**IMPORTANT NOTE 3:** We have provided **THREE** different **Bash** files. All files are available in the **"Bash Folder"**. 
<br>
<br>

## Table of Contents
- [Bash file without LD Clumping](#bash-file-without-ld-clumping)
- [Bash file with LD Clumping](#bash-file-with-ld-clumping)
- [Bash file - Global](#bash-file---global)
  * [Description of required parameters](#description-of-required-parameters)
  * [User specifications](#user-specifications)
  * [Different steps run by the bash file](#different-steps-run-by-the-bash-file)
- [Acknowledgements](#acknowledgements)
- [How to Cite](#how-to-cite)


## Bash file without LD Clumping
**Bash file Name:** 00_Global_run_GCPBayes_Strategy_Without_LDClumping_sizemax700.sh
<br><br>
This bash file runs the *GCPBayes pipeline* used in the paper (without LD clumping at all) on Breast (BCAC) and Ovarian (OCAC) GWAS summary statistics data.
<br>

## Bash file with LD Clumping 
**Bash file Name:** 00_Global_run_GCPBayes_Strategy_With_LDClumping_for_all_sizemax700.sh
<br><br>
This bash file runs the *GCPBayes pipeline* used in the paper (with LD clumping step) on Breast (BCAC) and Ovarian (OCAC) GWAS summary statistics data.
<br>

## Bash file - Global 
**Bash file Name:** 00_Global_run_GCPBayes.sh
<br><br>
We provided an overall bash file for performing analyses on Breast (BCAC) and Ovarian (OCAC) GWAS summary statistics data. However, a user could edit this Bash file and use it for running any GWAS summary statistics data two traits. 
<br>
This code is running the GCPBayes pipeline from **Section C** to **Section E**, by calling for different *R scripts* in a simple and automated way. 
<br>
The overall bash file and all other scripts of the pipeline, from **Section C** to **Section E**, have to be in the same folder. 
<br>
Then, a user can use this bash file to run the overall analysis once their inputs have been formatted. 
<br>
A user can use annotated files we already prepared or use his own annotation files.

### Description of required parameters
A short description of each parameter can be found at the beginning of the script, in the **“PARAMETERS”** section, as the default values of the parameters we propose to use when it is relevant. Also, a short description for inputs and outputs can be found in the **“INTPUTS”** and **“OUTPUTS”** sections.

### User specifications
At first, a user needs to specify different parameter inside the script that can be found in the **“USER SPECIFICATIONS”** section of the script. 
<br>
Here are some examples of parameters the user needs to fill (without default value):
-	Path of the directories the user wants to user as working directory and directory for final outputs of GCPBayes
-	Paths and file names for inputs of each dataset
-	Short names defining the dataset 1 and the dataset 2
-	Path and file names for annotated data that will be used for annotation of groups (genes or pathways) for GCPBayes
-	Path for folder containing the bfiles for reference data

Here are some examples of parameters with defaults values:
-	Quality control parameters: thresholds for quality of imputation (=0.9) and MAF (=0.05) to consider SNPs in the analysis
-	Threshold for theta values to consider further exploration with HS (=0.1), or to consider pleiotropy (=0.5)
-	Threshold for group length to consider LD clumping (=500), or to not consider for analysis (=1500)
-	Boolean value to consider LD Clumping or not (default value=”TRUE”)
-	Threshold for p-value in PLACO (=0.05)
-	Different parameters used for LD Clumping

**NOTE:** The values of parameters coded in the HARD CODED section of the script do not need to be modified.

### Different steps run by the bash file
**Bash file Name:** "00_Global_run_GCPBayes.sh” will run the different steps of the pipeline by using corresponding scripts. 
<br><br>
In summary, these steps are presented below:
- **Step 1:** To extract common SNPs between dataset by using “C1_code_find_common_snps_one_pair.R”
- **Step 2:** To run PLACO on each dataset by using “C2_code_run_PLACO_decor_one_pair.R”. This step is run only if Clumping=”TRUE”. Also, this step will not be performed if a corresponding output already exist in the working directory.
- **Step 3:** To perform LD Clumping (locally) by using “C3_code_ldclumping_local.R”. This step is run only if Clumping=”TRUE”.
- **Step 4:** To prepare the GCPBayes inputs in the right format by using Rscript “D1_code_pipeline_annot_coding_withoutldclumping_extra_info.R” and Rscript “D2_code_pipeline_annot_coding_ldclumping_extra_info.R”. The second code will be run only if Clumping=”TRUE”.
- **Step 5:** To get correct list of groups for further analysis according to wanted threshold (length of groups) by using “D3_code_separate_groups_length_threshold.R”
- **Step 6:** To run GCPBayes (DS Method) by using “E1_code_gcpbayes_less_extra_info.R”. If Clumping=”TRUE”, this step will be done separately in two different lists of groups (according to the lengths of the groups).
- **Step 7:** To run GCPBayes (HS Method) on groups with “theta > theta_exploration” only, by using “E2_code_gcpbayes_HS_less_extra_info.R” (work in progress)
- **Step 8:** To plot figures (work in progress)

## Acknowledgements 
We gratefully acknowledge the following packages/softwares which we used throughout our pipeline:
```
vroom
dplyr
data.table
devtools
MASS
tidyr
PLACO
genetics.binaRies
tidyverse
tictoc
GCPBayes
BhGLM
splitstackshape
PLINK
ieugwasr
gwasrapidd
CheckSumStats
ggplot2
plotly
patchwork
biomaRt
regioneR
karyoploteR
readxl
defaultdict
stats
datetime
```
## How to Cite
Asgari et al., "GCPBayes Pipeline: a tool for exploring pleiotropy at gene-level", xxxx. xxx x;x(x):x, doi: xxx[Paper_link](https://..../)

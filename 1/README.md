# Description of Bash file (a plain text file that contains a series of commands) for running the whole procedure with a series of options
Created by: Pierre-Emmanuel Sugier, Yazdan Asgari<br>
Creation date: 14 Jan 2022<br>
Update: Dec 2022<br>
https://cesp.inserm.fr/en/equipe/exposome-and-heredity
<br>
<br>


**IMPORTANT NOTE 1:** The Bash file was tested on a Unix-based server with CentOS 7.
<br>
<br>
**IMPORTANT NOTE 2:** All a user needs is to define parameters and Paths in the [**parameters.ini**](../0_Codes/Bash) and [**readinputs.txt**](../0_Codes/Bash) file and run a bash file for running on GWAS summary statistics data. 
<br>
<br>

## Table of Contents


## Bash file - General
**Bash file Name:** [00_Global_run_GCPBayes.sh](../0_Codes/Bash)
<br><br>
We developed a Bash file for running the GCPBayes pipeline from **Section C** to **Section E**, by calling for different *R scripts* in a simple and automated way. 
<br>
This bash file uses the scripts in the "Bash" folder. In order to run it, the input GWAS summary statistics data **MUST** have the following coulumns with **the same header names**:
| snp	| chr	| bp_hg19	| Effect_A | Baseline_A | beta | se | pval | info | EAF | MAF |
| -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- |

**snp** = SNP rs id, **chr** = chromosome, **bp_hg19** = base pair position in hg19 assembly, **Effect_A** = Effect Allele, **Baseline_A** = Baseline Allele, **beta** = beta value, **se** = standard error, **pval** = P-value, **EAF** = Effect Allele Frequency, **MAF** = Minor Allele Frequency
<br><br>
**NOTE:** **chr** column should be in a **numeric** format.
<br><br>

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
<br><br>

For running the Bash file, put all scripts (in this folder) and input GWAS data in a folder, then you need to make changes in two files (**parameters.ini** and **readinputs.txt**). Define the parameters and file paths in the **parameters.ini** file and also define the input files and their paths in the **readinputs.txt** file.
<br><br>
Then, to run, simply type in the terminal:
~~~
$ ./00_Global_run_GCPBayes.sh parameters.ini readinputs.txt
~~~

## Bash file - Parallel
**Bash file Name:** [00_Global_run_GCPBayes_Parallel.sh](../0_Codes/Bash_Parallel)
<br><br>

We also developed a Bash file which could run GCPBayes for more than one genes at the same time using different CPU cores. A user needs to define the **number of CPUs** in the **parameters.ini** file. All other options are the same as the **General Bash** file.
<br><br>
Then, to run, simply type in the terminal:
~~~
$ ./00_Global_run_GCPBayes_Parallel.sh parameters.ini readinputs.txt
~~~

**NOTE:** Using more CPUs also needs more RAMs for running the GCPBayes.

## Bash file for running on Breast and Ovarian Cancer data
For running the *GCPBayes pipeline* on Breast (BCAC) and Ovarian (OCAC) GWAS summary statistics data (appeared in our manuscript), all a user needs is to use the **"parameters_Strategy_bcac_ocac_manuscript.ini"** file.
<br><br>
In summary, it runs *GCPBayes package* **without ld clumping step** for genes with **less than 700 SNPs**, then runs *GCPBayes package* **with ld clumping step** for genes with **more than 700 SNPs**.
<br><br>
To run it, simply type in the terminal:
~~~
$ ./00_Global_run_GCPBayes.sh parameters_Strategy_bcac_ocac_manuscript.ini
~~~

## Acknowledgements 
We gratefully acknowledge the following packages/softwares which we used throughout our pipeline:
```
BiocManager
vroom
data.table
devtools
MASS
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

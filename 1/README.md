# Description of Bash file (a plain text file that contains a series of commands) for running the whole procedure with a series of options
Created by: Pierre-Emmanuel Sugier, Yazdan Asgari<br>
Creation date: 14 Jan 2022<br>
Update: May 2022<br>
https://cesp.inserm.fr/en/equipe/exposome-and-heredity
<br>
<br>


**IMPORTANT NOTE 1:** The Bash file was tested on a Unix-based server with CentOS 7.
<br>
<br>
**IMPORTANT NOTE 2:** All a user needs is to define parameters and Paths in the [**parameters.ini**](../0_Bash) file and run the bash file [**"00_Global_run_GCPBayes.sh"**](../0_Bash) for running any GWAS summary statistics data two traits. 
<br>
<br>

## Table of Contents
- [Bash file - Global](#bash-file---global)
  * [Description of required parameters](#description-of-required-parameters)
  * [Different steps run by the bash file](#different-steps-run-by-the-bash-file)
- [Bash file for running on Breast and Ovarian Cancer data](#bash-file-for-running-on-breast-and-ovarian-cancer-data)
- [Acknowledgements](#acknowledgements)
- [How to Cite](#how-to-cite)


## Bash file - Global 
**Bash file Name:** [00_Global_run_GCPBayes.sh](../0_Bash)
<br><br>
We provided an overall bash file for performing analyses on GWAS summary statistics data. A user could edit this **parameters.ini** file and use it for running any GWAS summary statistics data two traits. 
<br>
This code is running the GCPBayes pipeline from **Section C** to **Section E**, by calling for different *R scripts* in a simple and automated way. 
<br>
The overall bash file and all other scripts of the pipeline, from **Section C** to **Section E**, have to be in the same folder. 
<br>
Then, a user can use this bash file to run the overall analysis once their inputs have been formatted. 
<br>
A user can use annotated files we already prepared or use his own annotation files.
<br><br>
To run it, simply type in the terminal:
~~~
$ ./00_Global_run_GCPBayes.sh parameters.ini
~~~
### Description of required parameters
A description of each parameter can be found in [**parameters.ini**](../0_Bash) file. A user could change them based on own relevance. 
<br><br>
**NOTE:** The values of parameters coded in the HARD CODED section of the script do not need to be modified.

### Different steps run by the bash file
**Bash file Name:** "00_Global_run_GCPBayes.sh??? will run the different steps of the pipeline by using corresponding scripts. 
<br><br>
In summary, these steps are presented below:
- **Step 1:** To extract common SNPs between dataset by using ???C1_code_find_common_snps_one_pair.R???
- **Step 2:** To run PLACO on each dataset by using ???C2_code_run_PLACO_decor_one_pair.R???. This step is run only if Clumping=???TRUE???. Also, this step will not be performed if a corresponding output already exist in the working directory.
- **Step 3:** To perform LD Clumping (locally) by using ???C3_code_ldclumping_local.R???. This step is run only if Clumping=???TRUE???.
- **Step 4:** To prepare the GCPBayes inputs in the right format by using Rscript ???D1_code_pipeline_annot_coding_withoutldclumping_extra_info.R??? and Rscript ???D2_code_pipeline_annot_coding_ldclumping_extra_info.R???. The second code will be run only if Clumping=???TRUE???.
- **Step 5:** To get correct list of groups for further analysis according to wanted threshold (length of groups) by using ???D3_code_separate_groups_length_threshold.R???
- **Step 6:** To run GCPBayes (DS Method) by using ???E1_code_gcpbayes_less_extra_info.R???. If Clumping=???TRUE???, this step will be done separately in two different lists of groups (according to the lengths of the groups).
- **Step 7:** To run GCPBayes (HS Method) on groups with ???theta > theta_exploration??? only, by using ???E2_code_gcpbayes_HS_less_extra_info.R??? *(work in progress)*
- **Step 8:** To plot figures *(work in progress)*


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

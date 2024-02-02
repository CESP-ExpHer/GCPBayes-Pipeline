Created by: Yazdan Asgari<br>
Creation date: 07 Mar 2023<br>
Update: Mar 2023<br>
https://cesp.inserm.fr/en/equipe/exposome-and-heredity
<br>
<br>

# - How to use R script to run GCPBayes Pipeline
**NOTE 1:** It is needed to put all files and scripts in the same folder. 
<br><br>
**NOTE 2:** GWAS input files should be standardized (or harmonized) first. If you need to know how to standardize GWAS data, see [**Standardization Part**](../2) in the **"Tutorial - Wiki"** file.
<br><br>
**NOTE 3:** R script was tested on a Windows OS and Unix-based server with CentOS 7.
<br>
<br>

## 1- Checking required packages
First, you need to run the ([*"GCPBayes_pipeline_check_packages.R"*](../0_Codes/R)) R script to check a list of required packages and install them if they are not available in your system. It will also print a warning message if any of the packages could not be installed.  


## 2- Definition of parameters
Second, you need to change the parameters on the ([*"GCPBayes_pipeline_parameters.R"*](../0_Codes/R)) file. You just easily open this file in RStudio or a text editor and change the required parameters. Here are some recommendations:
- In the **"PATH SPECIFICATIONS"** section, replace **/PATH/** regarding *"working directory"* with the path where you put all downloaded data and scripts in your system.
- In the **"INPUT GWAS FILES NAMES"** section, replace *input* and *input_shortname* with the file names of your GWAS data. Put all names in a quotation and divide them by a comma.
- In the **"GWAS HEADERS NAMES"** section, you need to specify the headers of your GWAS file and the GCPBayes Pipeline will use them during the running. So, you do not need to make any changes to your GWAS input data anymore.
- In the **"Decision for running LD Clumping Step"** section, determine whether you want to perform the LD Clumping step or not (TRUE or FALSE).
- You could change all other parameters based on the definitions in their comments part. 

## 3- Running GCPBayes Pipeline
Finally, all you need is to run the ([*"GCPBayes_pipeline.R"*](../0_Codes/R)) R script. The script reads all parameters from the *"GCPBayes_pipeline_parameters.R"* file and runs every step of the Pipeline. You could run this script using RStudio. If you are using a Unix-based OS, you could run the script using the following command:
~~~
$ Rscript GCPBayes_pipeline.R
~~~

**GCPBayes Pipeline Steps in a summary** (considering LD Clumping step as well): 
<br>
- **Step 1:** Extracting common SNPs between two GWAS datasets. 
- **Step 2:** Running PLACO on both GWAS datasets.
- **Step 3:** Performing LD Clumping (locally).
- **Step 4:** Switch from SNP-level into Gene-level for GWAS data and preparation of GCPBayes inputs.
- **Step 5:** Separation of the list of groups for GCPBayes analysis according to a selected threshold.
- **Step 6:** Running GCPBayes (DS Method).
<br><br>

## 4- Visualization of the GCPBayes Results
After running the Pipeline, you could use one of the visualization tools (online or local shiny App) that we developed in order to create different plots from the output of the GCPBayes method. More information is available on the [**Visualization Part**](../README.md#visualization)


# - How to use Bash file to run GCPBayes Pipeline
**NOTE 1:** It is recommended to put the Bash file, parameter file, input files, and all scripts in the same folder. 
<br><br>
**NOTE 2:** GWAS input files should be standardized (or harmonized) first. If you need to know how to standardize GWAS data, see [**Standardization Part**](../2) in the **"Tutorial - Wiki"** file.
<br><br>
**NOTE 3:** The Bash file was tested on a Unix-based server with CentOS 7.
<br>
<br>

## Bash file - General
We developed a Bash file for running the GCPBayes pipeline from **Section C** to **Section E**, by calling for different *R scripts* in a simple and automated way. 
<br>
This bash file uses the scripts in the "Bash" folder. In order to run it, the input GWAS summary statistics data **MUST** have the following columns with **the same header names**:
| snp	| chr	| bp_hg19	| Effect_A | Baseline_A | beta | se | pval | info | EAF | MAF |
| -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- |

**snp** = SNP rs id, **chr** = chromosome, **bp_hg19** = base pair position in hg19 assembly, **Effect_A** = Effect Allele, **Baseline_A** = Baseline Allele, **beta** = beta value, **se** = standard error, **pval** = P-value, **EAF** = Effect Allele Frequency, **MAF** = Minor Allele Frequency
<br><br>
**NOTE:** **chr** column should be in a **numeric** format.
<br><br>

For running the Bash file, first, put [**all scripts**](../0_Codes/Bash) and input GWAS data in the same folder. Then, you need to make changes in two files (**parameters.ini** and **readinputs.txt**) as explained below. 
<br><br>

### Parameter file ([parameter.ini](../0_Codes/Bash))
All parameters for running the pipeline are written in this file. So, a user just needs to change some parts of the **“USER SPECIFICATIONS”** section. 
- You need to replace **“/PATH/”** with the path where you put all data and scripts (you need to change these three parts: working directory, output directory, and directory for scripts). 
- It is not needed to change **“Input files - datasets”** section. Because in the **“readinputs.txt”** file, you define paths and file names for GWAS input files.
- Then, you need to change the path and name of the **“Input files - annotation”** part based on the annotation files you are going to use in the analyses.
- If an LD clumping step is used, you need to determine a directory path for the reference GWAS data (*bfiles*) in the **“Input files - reference data”** part.
- In **“QC parameters”** and **“parameters”** parts, a user could change different parameters used while running the pipeline.
  - Quality control parameters: thresholds for quality of imputation (=0.9) and MAF (=0.05) to consider SNPs in the analysis
  - Threshold for theta values to consider further exploration with HS (=0.1), or to consider pleiotropy (=0.5)
  - Threshold for group length to consider LD clumping (=1500), or to not consider for running GCPBayes (=1500)
  - Threshold for p-value in PLACO (=0.05)
  - Boolean value to consider LD Clumping or not (default value=false)
  - Different parameters used for LD Clumping (r2, kb, and p)
**NOTE:** You do NOT need to change the **“HARD CODED”** part.

### Input file ([readinputs.txt](../0_Codes/Bash))
In this file, a user should define the paths and names of GWAS data used in the pipeline. It is also needed to define a short name for each GWAS data which would be used for the creation of output file names during the pipeline. So, in **“readinputs.txt”** file, the first row is the header of the file and a user needs to change other rows while each row belongs to one GWAS dataset.

### Main Bash file ([00_Global_run_GCPBayes.sh](../0_Codes/Bash))
To run the Bash file, a user just needs to run the following command in the terminal:

~~~
$ ./00_Global_run_GCPBayes.sh parameters.ini readinputs.txt
~~~

**NOTE:** It might need to change the permission of the Bash file in order to be executed.
~~~
$ chmod 777 00_Global_run_GCPBayes.sh
~~~

**Bash file Name** *"00_Global_run_GCPBayes.sh”* will run 6 steps in the case of running all sections available in the script. Here, we provide all the steps in summary: 
<br>
- **Step 1:** Extracting common SNPs between two GWAS datasets. 
- **Step 2:** Running PLACO on both GWAS datasets.
- **Step 3:** Performing LD Clumping (locally).
- **Step 4:** Switch from SNP-level into Gene-level for GWAS data and preparation of GCPBayes inputs.
- **Step 5:** Separation of the list of groups for GCPBayes analysis according to a selected threshold.
- **Step 6:** Running GCPBayes (DS Method).
<br><br>

## Bash file - Parallel
We also developed a Bash file ([**00_Global_run_GCPBayes_Parallel.sh**](../0_Codes/Bash_Parallel)) which could run GCPBayes for more than one gene at the same time using different CPU cores. 
<br>
A user needs to define the **number of CPUs** in the [**parameters_parallel.ini**](../0_Codes/Bash_Parallel) file. All other options are the same as the **General Bash** file (explained above).
<br><br>
Then, to run, simply type in the terminal:
~~~
$ ./00_Global_run_GCPBayes_Parallel.sh parameters_parallel.ini readinputs.txt
~~~

**NOTE:** It might need to change the permission of the Bash file in order to be executed.
~~~
$ chmod 777 00_Global_run_GCPBayes_Parallel.sh
~~~


**IMPORTNAT NOTE:** Using **more CPUs** also needs **more RAMs** for running the GCPBayes. So, a user MUST consider this before setting the number of CPUs to a large number.

## Visualization of the GCPBayes Results
After running the Pipeline, you could use one of the visualization tools (online or local shiny App) that we developed in order to create different plots from the output of the GCPBayes method. More information is available on the [**Visualization Part**](../README.md#visualization)

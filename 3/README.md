# How to use R script to run GCPBayes Pipeline
Created by: Yazdan Asgari<br>
Creation date: 07 Mar 2023<br>
Update: Mar 2023<br>
https://cesp.inserm.fr/en/equipe/exposome-and-heredity
<br>
<br>

**NOTE 1:** It is needed to put all files and scripts in the same folder. 
<br><br>
**NOTE 2:** GWAS input files should be standardized (or harmonized) first. If you need to know how to standardize GWAS data, see [**Standardization Part**](../2) in the **"Tutorial - Wiki"** file.
<br><br>
**NOTE 3:** R script was tested on a Windows OS and Unix-based server with CentOS 7.
<br>
<br>

## 1- Checking required packages
First, you need to run the *"GCPBayes_pipeline_check_packages.R"* [(Link)](../0_Codes/R) R script to check a list of required packages and install them if they are not available in your system. It will also print a warning message if any of the packages could not be installed.  


## 2- Definition of parameters
Second, you need to change the parameters on the *"GCPBayes_pipeline_parameters.R"* [(Link)](../0_Codes/R) file. You just easily open this file in RStudio or a text editor and chenge the required parameters. Here are some recommendations:
- In the **"PATH SPECIFICATIONS"** section, replace **/PATH/** regarding *"working directory"* with the path where you put all downloaded data and scripts in your system.
- In the **"INPUT GWAS FILES NAMES"** section, replace *input* and *input_shortname* with the file names of your GWAS data. Put all names in a quotation and devide them by comma.
- In the **"GWAS HEADERS NAMES"** section, you need to specify the headers of your GWAS file and the GCPBayes Pipeline will use them during the running. So, you do not need to make any change to your GWAS input data any more.
- In the **"Decision for running LD Clumping Step"** section, determine whether you want to perform LD Clumping step or not (TRUE or FALSE).
- You could change all other parameters based on the definitions in their comments part. 

## 3- Running GCPBayes Pipeline
Finally, all you need is to run the *"GCPBayes_pipeline.R"* [(Link)](../0_Codes/R) R script. The script read all parameters from the *"GCPBayes_pipeline_parameters.R"* file and run every steps of the Pipeline. You could run this script using RStudio. If you are using a Unix-based OS, you could run the script using the following command:
~~~
$ Rscript GCPBayes_pipeline.R
~~~

**GCPBayes Pipeline Steps in a summary** (considering LD Clumping step as well): 
<br>
- **Step 1:** Extracting common SNPs between two GWAS datasets. 
- **Step 2:** Running PLACO on both GWAS datasets.
- **Step 3:** Performing LD Clumping (locally).
- **Step 4:** Switch from SNP-level into Gene-level for GWAS data and preparation of GCPBayes inputs.
- **Step 5:** Separation of list of groups for GCPBayes analysis according to a selected threshold.
- **Step 6:** Running GCPBayes (DS Method).
<br><br>

## 4- Visualization of the GCPBayes Results
After running the Pipeline, you could use one of the visualization tools (online or local shiny App) that we developed in order to create different plots from the output of the GCPBayes method. More information is available on the [**Visualization Part**](../README.md#visualization)


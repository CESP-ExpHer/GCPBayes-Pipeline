# How to use Shiny App to run GCPBayes Pipeline
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
**NOTE 3:** The Shiny App was tested on a Windows OS and Unix-based server with CentOS 7.
<br>
<br>

## 1- Downloading and extraction of the pipeline
For non-computer scientists, we have designed a Shiny App. All you need to do is download a folder called [**"Shiny_zipped"**](../0_Codes/Shiny_zipped). To download, click on the **"shiny_files_GitHub.zip"** file and then, click on **"View raw"** or click on the **"Raw"** button at the top-right of the page. After extraction of the folder, put your GWAS input and annotation files in the same folder. 

## 2-	Running “GCPBayes_pipeline_shiny_v1.2.R” file
In the **"shiny_files_GitHub"** folder, open **"GCPBayes_pipeline_shiny_v1.2.R"** file using RStudio.
<br>
<br>
**NOTE:** When you run the **"GCPBayes_pipeline_shiny_v1.2.R"** file, it uses the path that the **"shiny_files_GitHub"** folder exists as a default path of the analysis. 

## 3- Checking required packages
By clicking on the **"Check Required Packages"** button, the App checks a list of required packages for the pipeline. Missed packages will be installed automatically and it also prints a warning message if any of the packages could not be installed. 

## 4- Definition of parameters
You could easily change any parameter (such as working directory path, GWAS column names, threshold values, etc.) by seeing its current value in a GUI (Graphical user interface). And just need to click **“Update Parameters File”** button to apply all changes 

## 5- Running GCPBayes Pipeline
By clicking on the **"Run GCPBayes Pipeline"** button, the App reads all parameters from the **"GCPBayes_pipeline_parameters.R"** file and runs every step of the Pipeline.

**GCPBayes Pipeline Steps in a summary** (considering LD Clumping step as well): 
<br>
- **Step 1:** Extracting common SNPs between two GWAS datasets. 
- **Step 2:** Running PLACO on both GWAS datasets.
- **Step 3:** Performing LD Clumping (locally).
- **Step 4:** Switch from SNP-level into Gene-level for GWAS data and preparation of GCPBayes inputs.
- **Step 5:** Separation of the list of groups for GCPBayes analysis according to a selected threshold.
- **Step 6:** Running GCPBayes (DS Method).
<br><br>

## 6- Visualization of the GCPBayes Results
After running the Pipeline, a user could use one of the visualization tools (online or local shiny App) that we developed in order to create different plots from the output of the GCPBayes method. More information is available on the [**Visualization Part**](../README.md#visualization)




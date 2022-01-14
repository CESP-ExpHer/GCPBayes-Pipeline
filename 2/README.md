# A Step-by-Step Tutorial for Analyses of BCAC and OCAC GWAS Summary Statistics Data
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
**IMPORTANT NOTE 1:** All scripts all available in the [**"Codes folder"**](../0_Codes)
<br><br>
**IMPORTANT NOTE 2:** In order to follow all procedure easier, we considered all paths in the same directory (for inputs and outputs) throughout the pipeline (in our example: “~/BCAC_OCAC/”)
<br><br>
**IMPORTANT NOTE 3:** The Section names that are mentioned in this page are correspondence to Figure provided in the [An Overview of the Pipeline](/).
<br>

## Table of Contents
- [Summary](#summary)
  * [Running the GCPBayes pipeline without LD Clumping](#running-the-gcpbayes-pipeline-without-ld-clumping)
  * [Running the GCPBayes pipeline with LD Clumping](#running-the-gcpbayes-pipeline-with-ld-clumping)
- [Standardization-Section A](#standardization-section-a)
  * [First Step](#first-step)
  * [Second Step](#second-step)
- [Annotation-Section B](#annotation-section-b)
  * [First Step](#first-step-1)
  * [Second Step](#second-step-1)
- [LD Clumping-Section C](#ld-clumping-section-c)
  * [First Step](#first-step-2)
  * [Second Step](#second-step-2)
  * [Third Step](#third-step)
- [Core and Running GCPBayes without LD Clumping](#core-and-running-gcpbayes-without-ld-clumping)
  * [Section D](#section-d)
  * [Section E](#section-e)
- [Core and Running GCPBayes with LD Clumping](#core-and-running-gcpbayes-with-ld-clumping)
  * [Section D](#section-d-1)
  * [Section E](#section-e-1)
- [Visualization](#visualization)
  * [Section A](#section-a)
  * [Section B](#section-b)
  * [Section C](#section-c)
  * [Section D](#section-d-2)
  * [Section E](#section-e-2)
- [References](#references)
- [Acknowledgements](#acknowledgements)


## Summary
An overall summary of running the GCPBayes pipeline for BCAC and OCAC GWAS summary statistics data is provided in the following **TWO Tables**:
<br>
### Running the GCPBayes pipeline without LD Clumping
| No	| Running the Scripts (in order)	| Section Name	| Program |
| -- | -- | -- | -- | 
| 1	| A1_code_reformatting_file_bcac_2020_all.py	| Standardization (A)	| Python |
| 2	| A2_code_reformatting_file_ocac_bcac_2020_all.py	| Standardization (A)	| Python |
| 3	| Annotation Step based on the explanations in the GitHub page	| Annotation (B1)	| R |
| 4	| Annotation Step based on the section 3.2.1	| Annotation (B2)	| R, PLINK |
| 5	| D1_code_pipeline_annot_coding_withoutldclumping_extra_info.R	| Core (D)	| R |
| 6	| E1_code_gcpbayes_less_extra_info.R	| Running GCPBayes (E)	| R |
| 7	| E2_code_gcpbayes_greater_extra_info.R	| Running GCPBayes (E)	| R |
| 8	| code_analysis_A_checksumstats_BCAC.R	| Visualization (A)	| R |
| 9	| code_analysis_A_checksumstats_OCAC.R	| Visualization (A)	| R |
| 10	| code_analysis_B_annotation_file_coding.R	| Visualization (B)	| R |
| 11	| code_analysis_D_gcpbayes_input_ggplot.R	| Visualization (D)	| R |
| 12	| code_analysis_E_gcpbayes_output_karyotype.R	| Visualization (E)	| R |
| 13	| code_analysis_E_gcpbayes_output_table_overview.R	| Visualization (E)	| R |

### Running the GCPBayes pipeline with LD Clumping
| No	| Running the Scripts (in order)	| Section Name	| Program |
| -- | -- | -- | -- | 
| 1	| A1_code_reformatting_file_bcac_2020_all.py	| Standardization (A)	| Python |
| 2	| A2_code_reformatting_file_ocac_bcac_2020_all.py	| Standardization (A)	| Python |
| 3	| Annotation Step based on the explanations in the GitHub page	| Annotation (B1)	| R |
| 4	| Annotation Step based on the section 3.2.1	| Annotation (B2)	| R, PLINK |
| 5	| C1_code_find_shared_snps_one_pair.R	| LD Clumping (C)	| R |
| 6	| C2_code_run_PLACO_decor_one_pair.R	| LD Clumping (C)	| R |
| 7	| C3_code_ldclumping_local.R	| LD Clumping (C)	| R |
| 8	| D2_code_pipeline_annot_coding_ldclumping_extra_info.R	| Core (D)	| R |
| 9	| E1_code_gcpbayes_less_extra_info.R	| Running GCPBayes (E)	| R |
| 10	| E2_code_gcpbayes_greater_extra_info.R	| Running GCPBayes (E)	| R |
| 11	| code_analysis_A_checksumstats_BCAC.R	| Visualization (A)	| R |
| 12	| code_analysis_A_checksumstats_OCAC.R	| Visualization (A)	| R |
| 13	| code_analysis_B_annotation_file_coding.R	| Visualization (B)	| R |
| 14	| code_analysis_C_PLACO_results_one_pair.R	| Visualization (C)	| R |
| 15	| code_analysis_D_gcpbayes_input_ggplot.R	| Visualization (D)	| R |
| 16	| code_analysis_E_gcpbayes_output_karyotype.R	| Visualization (E)	| R |
| 17	| code_analysis_E_gcpbayes_output_table_overview.R	| Visualization (E)	| R |

Here are detailed information for running each step which a user should RESPECTIVELY run in order to get the results provided in the manuscript: 
## Standardization-Section A 
***(Section A) (Python)***
### First Step 
***(First GWAS Summary Statistics Data called as a Reference file):***
<br>
**Script:** A1_code_reformatting_file_bcac_2020_all.py
<br>
**Input**
- Input files needed for running this step:
  - icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt

Consider that BCAC GWAS summary statistics data (v. 2020) (Zhang et al., 2020) downloaded from the [Consortium Webpage](https://bcac.ccge.medschl.cam.ac.uk/bcacdata/oncoarray/). We used the python script to extract the following columns and rename them as follows:
<br>
- File name: icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt
  - column 24: snp (RS ID)
  - column 25: chr (chromosome)
  - column 26: bp_hg19 (base pair position)
  - column 27: Effect_A (Effect Allele)
  - column 28: nonEffect_A (non-Effect Allele)
  - column 29: EAF (Effect Allele Frequency)
  - column 31: info (r2 value)
  - column 33: beta (beta value)
  - column 34: se (standard error)
  - column 38: pval (P-value)

**NOTE:** A user MUST modify the **“DEFINITION SECTION”** which includes the following options:
- Path for the input file: (in our example “~/BCAC_OCAC/”)
- Input file name: (in our example “icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt”)
- Path for the output file: (in our example “~/BCAC_OCAC/”)

**Output**
- Output files created after running this step:
  - BCAC_2020_onco_ALL_reformatted.txt
  - BCAC_2020_onco_ALL_SNP_ambiguous.txt
  - BCAC_2020_onco_ALL_SNP_dupicated.txt
  - BCAC_2020_onco_ALL_SNP_duplication_set.txt
  - BCAC_2020_onco_ALL_SNP_weird_alleles.txt
  - BCAC_2020_onco_ALL_Summary.txt

### Second Step 
***(Second GWAS Summary Statistics Data):***
<br>
**Script:** A2_code_reformatting_file_ocac_bcac_2020_all.py
<br>
**Input**
- Input files needed for running this step:
  - BCAC_2020_onco_ALL_reformatted.txt (created from the First Step)
  - extraction_OCAC.txt

Consider that OCAC GWAS summary statistics data (Phelan et al., 2017). We used the python script to extract the following columns and rename them as follows:
<br>
- File name: extraction_OCAC.txt
  - column 2: chr (chromosome)
  - column 3: bp_hg19 (base pair position)
  - column 4: Effect_A (Effect Allele)
  - column 5: nonEffect_A (non-Effect Allele)
  - column 6: EAF (Effect Allele Frequency)
  - column 7: nEAF (non-Effect Allele Frequency)
  - column 8: beta (beta value)
  - column 9: se (standard error)
  - column 10: pval (P-value)
  - column 11: info (r2 value)

**NOTE:** A user MUST modify the **“DEFINITION SECTION”** which includes the following options:
- Path for the GWAS reference file: (in our example “~/BCAC_OCAC/”)
- Reference file name: (in our example “BCAC_2020_onco_ALL_reformatted.txt”)
- Path for the input file: (in our example “~/BCAC_OCAC/”)
- Input file name: (in our example “extraction_OCAC.txt”)
- Path for the output file: (in our example “~/BCAC_OCAC/”)
- First part of the output file names: (in our example “OCAC_BCAC_2020_onco_ALL”)

**Output**
- Output files created after running this step:
  - OCAC_BCAC_2020_onco_ALL_reformatted.txt
  - OCAC_BCAC_2020_onco_ALL_SNP_ambiguous.txt
  - OCAC_BCAC_2020_onco_ALL_SNP_dupicated.txt
  - OCAC_BCAC_2020_onco_ALL_SNP_duplication_set.txt
  - OCAC_BCAC_2020_onco_ALL_SNP_removed_from_file.txt
  - OCAC_BCAC_2020_onco_ALL_SNP_removed_from_ref.txt
  - OCAC_BCAC_2020_onco_ALL_SNP_weird_alleles.txt
  - OCAC_BCAC_2020_onco_ALL_Summary_SNP.txt

## Annotation-Section B 
***(Section B) (R, PLINK)***
<br>
This section creates two files (one would be used in Section D4 and the other used in Section D2).
### First Step 
To create an annotation file that would be used in the Section D4, a user should follow the procedure explained in details in our [Annotation GitHub page](https://github.com/CESP-ExpHer/Gene_Annotation/tree/main/1_hg19)
<br>
We used the annotation file downloaded from GENCODE webpage (Frankish et al., 2019) and just used protein-coding genes for this study.
<br>
**Input**
- Input files needed for this step:
  - REMAP_gencode.v28lift37.annotation.gtf

**Output**
- Output file created after running this step:
  - annot_gencode_v38lift37_modified_gene_class.txt
### Second Step
To create an annotation file that would be used in the Section D2, a user should run the procedures mentioned in the section XXX of our [GitHub page](XXX) (which consists of some commands in R for preparation of input files, one annotation command in PLINK, and finally some post-processing commands in R on the PLINK annotation output file).
<br>
**Input**
- Input files needed for this step:
  - BCAC_2020_onco_ALL_reformatted.txt (created from Section A)
  - annot_gencode_v38lift37_modified_gene_class.txt (created from First Step)

**Output**
- Output files created after running this step:
  - annot_gencode_v38lift37_modified_gene_class_coding.txt
  - annot_gencode_v38lift37_modified_gene_class_coding_chr_num_plink_input.txt
  - BCAC_2020_onco_ALL_reformatted.assoc
  - plink.annot
  - Annot_BCAC_2020_onco_ALL_reformatted_coding.txt 

## LD Clumping-Section C 
***(Section C) (R)***
<br>
**NOTE:** If a user do not want to consider LD clumping through the pipeline, it could ignore this section and move to the next section.
### First Step 
***(finding shared SNPs between two traits)***
<br>
**Script:** C1_code_find_shared_snps_one_pair.R
<br>
**Input**
- Input files needed for running this step:
  - BCAC_2020_onco_ALL_reformatted.txt (created from Section A)
  - OCAC_BCAC_2020_onco_ALL_reformatted.txt (created from Section A)

**NOTE:** If a user performed the standardization (Section A) using our scripts, the column names are correct. Otherwise, a user MUST consider that the following column names are needed to be available during the running process. So, a user might need rename them in GWAS data or make change in the script. The column names used in the script are as follows:
  - snp (RS ID)
  - chr (chromosome)
  - bp_hg19 (base pair position)
  - Effect_A (Effect Allele)
  - nonEffect_A (non-Effect Allele)
  - beta (beta value)
  - se (standard error)
  - pval (P-value)
  - info (r2 value)
  - EAF (Effect Allele Frequency)
  - MAF (Minor Allele Frequency)

**NOTE:** A user MUST modify the **“DEFINITION SECTION”** which includes the following options:
- Path for two GWAS input files: (in our example “~/BCAC_OCAC/”)
- Input file names: (in our example “BCAC_2020_onco_ALL_reformatted.txt” and “OCAC_BCAC_2020_onco_ALL_reformatted.txt”)
- Path for the output file: (in our example “~/BCAC_OCAC/”)
- Output file names: (in our example “BCAC_2020_ALL_Shared_OCAC_inc_Z.txt” and “OCAC_Shared_BCAC_2020_ALL_inc_Z.txt”)
- info_threshold (used for filtering GWAS data) (in our example = 0.8)
- MAF_threshold (used for filtering GWAS data) (in our example = 0.01)

**Output**
- Output files created after running this step:
  - BCAC_2020_ALL_Shared_OCAC_inc_Z.txt
  - OCAC_Shared_BCAC_2020_ALL_inc_Z.txt

### Second Step 
***(running PLACO)***
<br>
**Script:** C2_code_run_PLACO_decor_one_pair.R
<br>
**Input**
- Input files needed for running this step:
  - BCAC_2020_ALL_Shared_OCAC_inc_Z.txt (created from Section C1.3.1)
  - OCAC_Shared_BCAC_2020_ALL_inc_Z. txt (created from Section C1.3.1)

**NOTE:** A user MUST modify the **“DEFINITION SECTION”** which includes the following options:
- Path for two input files: (in our example “~/BCAC_OCAC/”)
- Input file names: (in our example “BCAC_2020_ALL_Shared_OCAC_inc_Z.txt” and “OCAC_Shared_BCAC_2020_ALL_inc_Z.txt”)
- Path for the output file: (in our example “~/BCAC_OCAC/”)
- Output file names: (in our example “output_PLACO_BCAC_2020_ALL_OCAC”)
- pval_threshold (used for decorrelating the Z-scores) (in our example = 0.05)

**Output**
- Output file created after running this step:
  - output_PLACO_BCAC_2020_ALL_OCAC.txt 

### Third Step 
***(running LD Clumping)***
<br>
**Script:** C3_code_ldclumping_local.R
<br>
**Input**
- Input files needed for running this step:
  - output_PLACO_BCAC_2020_ALL_OCAC.txt (created from Section C1.3.2)

**NOTE:** A user MUST modify the **“DEFINITION SECTION”** which includes the following options:
- Path for the input file: (in our example “~/BCAC_OCAC/”)
- PLACO Input file name: (in our example “output_PLACO_BCAC_2020_ALL_OCAC.txt”) (created from Section C1.3.2)
- Path for the output file: (in our example “~/BCAC_OCAC/”)
- Output file name: (in our example “output_ld_clumping_08_BCAC_2020_ALL_OCAC.txt”)

**NOTE:** A user MUST also modify the LD Calculation parameters in the **“RUNNING SECTION”** which includes the following options:
- clump_r2 (Clumping r2 threshold) (in our example = 0.8)
- pop (super-population to use as reference panel) (in our example = “EUR”)
- bfile (the PATH where the super-population files exists) (in our example “~/BCAC_OCAC/1000K_ref_data/EUR”)

**Output**
- Output file created after running this step:
- output_ld_clumping_08_BCAC_2020_ALL_OCAC.txt 

## Core and Running GCPBayes without LD Clumping
***(Section D without LD Clumping + Section E) (R)***
### Section D
**Script:** D1_code_pipeline_annot_coding_withoutldclumping_extra_info.R
<br>
**Input**
- Input files needed for running this step:
  - BCAC_2020_onco_ALL_reformatted.txt (created from Section A)
  - OCAC_BCAC_2020_onco_ALL_reformatted.txt (created from Section A)
  - Annot_BCAC_2020_onco_ALL_reformatted_coding.txt (created from Section B)
  - annot_gencode_v38lift37_modified_gene_class.txt (created from Section B)

**NOTE:** A user MUST modify the **“DEFINITION SECTION”** which includes the following options:
- Path for two GWAS input files: (in our example “~/BCAC_OCAC/”)
- GWAS input file names: (in our example “BCAC_2020_onco_ALL_reformatted.txt” and “OCAC_BCAC_2020_onco_ALL_reformatted.txt”)
- First GWAS Annotation file path (in our example “~/BCAC_OCAC/”)
- First GWAS Annotation file name (in our example “Annot_BCAC_2020_onco_ALL_reformatted_coding.txt”)
- Second GWAS Annotation file path (in our example “~/BCAC_OCAC/”)
- Second GWAS Annotation file name (in our example “annot_gencode_v38lift37_modified_gene_class.txt”)
- Path for saving/reading files created by this script (called “path_input_data” and “path_output_data”) (in our example “~/BCAC_OCAC/”)
- Definition of some file names used during running the script in different steps. In our example are:
  - output_pipeline_BCAC_ALL_Shared_OCAC_coding_withoutclumping
  - output_pipeline_OCAC_Shared_BCAC_ALL_coding_withoutclumping
  - output_pipeline_BCAC_ALL_OCAC_coding_withoutclumping
  - output_pipeline_Summary_SNP_in_genes_BCAC_ALL_OCAC_coding_withoutclumping
  - output_pipeline_Summary_SNP_in_genes_OCAC_BCAC_ALL_coding_withoutclumping
- info_threshold (used for filtering GWAS data) (in our example = 0.9)
- MAF_threshold (used for filtering GWAS data) (in our example = 0.05)

**Output**
- Output files created after running this step:
  - output_pipeline_BCAC_ALL_Shared_OCAC_coding_withoutclumping.txt
  - output_pipeline_OCAC_Shared_BCAC_ALL_coding_withoutclumping.txt
  - output_pipeline_BCAC_ALL_OCAC_coding_withoutclumping.txt
  - output_pipeline_Summary_SNP_in_genes_BCAC_ALL_OCAC_coding_withoutclumping.txt
  - output_pipeline_Summary_SNP_in_genes_OCAC_BCAC_ALL_coding_withoutclumping.txt
  - Matrices_output_pipeline_BCAC_ALL_OCAC_coding_withoutclumping.Rdata
  - Matrices_extra_info_output_pipeline_BCAC_ALL_OCAC_coding_withoutclumping.Rdata

### Section E 
***(for genes with number of SNPs less than a threshold)***
<br>
**Script:** E1_code_gcpbayes_less_extra_info.R
<br>
**Input**
- Input files needed for running this step:
  - Matrices_output_pipeline_BCAC_ALL_OCAC_coding_withoutclumping.Rdata (created from Section 1.4.1)
  - Matrices_extra_info_output_pipeline_BCAC_ALL_OCAC_coding_withoutclumping.Rdata (created from Section 1.4.1)

**NOTE:** A user MUST modify the **“DEFINITION SECTION”** which includes the following options:
- Path for two GCPBayes input files: (in our example “~/BCAC_OCAC/”)
- GCPBayes input file names: (in our example “Matrices_output_pipeline_BCAC_ALL_OCAC_coding_withoutclumping.Rdata” and “Matrices_extra_info_output_pipeline_BCAC_ALL_OCAC_coding_withoutclumping.Rdata”)
- Path for the output files: (in our example “~/BCAC_OCAC/”)
- Output file names. In our example are:
  - output_GCPBayes_BCAC_All_OCAC_coding_withoutclumping_less_threshold_500_results
  - output_GCPBayes_BCAC_All_OCAC_coding_withoutclumping_less_threshold_500_ pleiotropy
  - output_GCPBayes_BCAC_All_OCAC_coding_withoutclumping_less_threshold_500_ errors
- SNP_number_threshold (A threshold for splitting genes based on SNP numbers) (in our example = 500)
- theta_threshold (for considering genes with potential pleiotropic effect) (in our example = 0.5)

**Output**
- Output files created after running this step:
  - output_GCPBayes_BCAC_All_OCAC_coding_withoutclumping_less_threshold_500_results.txt
  - output_GCPBayes_BCAC_All_OCAC_coding_withoutclumping_less_threshold_500_ pleiotropy.txt (if any genes with pleiotropic signals found)
  - output_GCPBayes_BCAC_All_OCAC_coding_withoutclumping_less_threshold_500_ errors.txt (if any error occurred during calculation for each gene)

***(for genes with number of SNPs greater than a threshold)***
<br>
**Script:** E2_code_gcpbayes_greater_extra_info.R
<br>
**Input**
- Input files needed for running this step:
  - Matrices_output_pipeline_BCAC_ALL_OCAC_coding_withoutclumping.Rdata (created from Section 1.4.1)
  - Matrices_extra_info_output_pipeline_BCAC_ALL_OCAC_coding_withoutclumping.Rdata (created from Section 1.4.1)

**NOTE:** A user MUST modify the **“DEFINITION SECTION”** which includes the following options:
- Path for two GCPBayes input files: (in our example “~/BCAC_OCAC/”)
- GCPBayes input file names: (in our example “Matrices_output_pipeline_BCAC_ALL_OCAC_coding_withoutclumping.Rdata” and “Matrices_extra_info_output_pipeline_BCAC_ALL_OCAC_coding_withoutclumping.Rdata”)
- Path for the output files: (in our example “~/BCAC_OCAC/”)
- Output file names. In our example are:
  - output_GCPBayes_BCAC_All_OCAC_coding_withoutclumping_greater_threshold_500_results
  - output_GCPBayes_BCAC_All_OCAC_coding_withoutclumping_greater_threshold_500_pleiotropy
  - output_GCPBayes_BCAC_All_OCAC_coding_withoutclumping_greater_threshold_500_errors
- SNP_number_threshold (A threshold for splitting genes based on SNP numbers) (in our example = 500)
- theta_threshold (for considering genes with potential pleiotropic effect) (in our example = 0.5)

**Output**
- Output files created after running this step:
  - output_GCPBayes_BCAC_All_OCAC_coding_withoutclumping_greater_threshold_500_results.txt
  - output_GCPBayes_BCAC_All_OCAC_coding_withoutclumping_greater_threshold_500_pleiotropy.txt (if any genes with pleiotropic signals found)
  - output_GCPBayes_BCAC_All_OCAC_coding_withoutclumping_greater_threshold_500_errors.txt (if any error occurred during calculation for each gene)

## Core and Running GCPBayes with LD Clumping
***(Section D with LD Clumping + Section E) (R)***
### Section D
**Script:** D2_code_pipeline_annot_coding_ldclumping_extra_info.R
<br>
**Input**
- Input files needed for running this step:
  - BCAC_2020_onco_ALL_reformatted.txt (created from Section A)
  - OCAC_BCAC_2020_onco_ALL_reformatted.txt (created from Section A)
  - Annot_BCAC_2020_onco_ALL_reformatted_coding.txt (created from Section B)
  - annot_gencode_v38lift37_modified_gene_class.txt (created from Section B)
  - output_ld_clumping_08_BCAC_2020_ALL_OCAC.txt (created from Section C)

**NOTE:** A user MUST modify the **“DEFINITION SECTION”** which includes the following options:
- Path for two GWAS input files: (in our example “~/BCAC_OCAC/”)
- GWAS input file names: (in our example “BCAC_2020_onco_ALL_reformatted.txt” and “OCAC_BCAC_2020_onco_ALL_reformatted.txt”)
- First GWAS Annotation file path (in our example “~/BCAC_OCAC/”)
- First GWAS Annotation file name (in our example “Annot_BCAC_2020_onco_ALL_reformatted_coding.txt”)
- Second GWAS Annotation file path (in our example “~/BCAC_OCAC/”)
- Second GWAS Annotation file name (in our example “annot_gencode_v38lift37_modified_gene_class.txt”)
- Path for LD Clumping input file: (in our example “~/BCAC_OCAC/”)
- LD Clumping file names: (in our example “output_ld_clumping_08_BCAC_2020_ALL_OCAC.txt”
- Path for saving/reading files created by this script (called “path_input_data” and “path_output_data”) (in our example “~/BCAC_OCAC/”)
- Definition of some file names used during running the script in different steps. In our example are:
  - output_pipeline_BCAC_ALL_Shared_OCAC_coding_clumping_08 
  - output_pipeline_OCAC_Shared_BCAC_ALL_coding_clumping_08 
  - output_pipeline_BCAC_ALL_OCAC_coding_clumping_08 
  - output_pipeline_Summary_SNP_in_genes_BCAC_ALL_OCAC_coding_clumping_08 
  - output_pipeline_Summary_SNP_in_genes_OCAC_BCAC_ALL_coding_clumping_08 
- info_threshold (used for filtering GWAS data) (in our example = 0.9)
- MAF_threshold (used for filtering GWAS data) (in our example = 0.05)

**Output**
- Output files created after running this step:
  - output_pipeline_BCAC_ALL_Shared_OCAC_coding_clumping_08.txt 
  - output_pipeline_OCAC_Shared_BCAC_ALL_coding_clumping_08.txt 
  - output_pipeline_BCAC_ALL_OCAC_coding_clumping_08.txt 
  - output_pipeline_Summary_SNP_in_genes_BCAC_ALL_OCAC_coding_clumping_08.txt 
  - output_pipeline_Summary_SNP_in_genes_OCAC_BCAC_ALL_coding_clumping_08.txt 
  - Matrices_output_pipeline_BCAC_ALL_OCAC_coding_clumping_08.Rdata
  - Matrices_extra_info_output_pipeline_BCAC_ALL_OCAC_coding_clumping_08.Rdata

### Section E 
***(for genes with number of SNPs less than a threshold)***
<br>
**Script:** E1_code_gcpbayes_less_extra_info.R
<br>
**Input**
- Input files needed for running this step:
  - Matrices_output_pipeline_BCAC_ALL_OCAC_coding_clumping_08.Rdata (created from Section 1.4.1)
  - Matrices_extra_info_output_pipeline_BCAC_ALL_OCAC_coding_clumping_08.Rdata (created from Section 1.4.1)

**NOTE:** A user MUST modify the **“DEFINITION SECTION”** which includes the following options:
- Path for two GCPBayes input files: (in our example “~/BCAC_OCAC/”)
- GCPBayes input file names: (in our example “Matrices_output_pipeline_BCAC_ALL_OCAC_coding_clumping_08.Rdata” and “Matrices_extra_info_output_pipeline_BCAC_ALL_OCAC_coding_clumping_08.Rdata”)
- Path for the output files: (in our example “~/BCAC_OCAC/”)
- Output file names. In our example are:
  - output_GCPBayes_BCAC_All_OCAC_coding_ldclumping_08_less_threshold_500_results
  - output_GCPBayes_BCAC_All_OCAC_coding_ldclumping_08_less_threshold_500_pleiotropy
  - output_GCPBayes_BCAC_All_OCAC_coding_ldclumping_08_less_threshold_500_errors
- SNP_number_threshold (A threshold for splitting genes based on SNP numbers) (in our example = 500)
- theta_threshold (for considering genes with potential pleiotropic effect) (in our example = 0.5)

**Output**
- Output files created after running this step:
  - output_GCPBayes_BCAC_All_OCAC_coding_ldclumping_08_less_threshold_500_results.txt
  - output_GCPBayes_BCAC_All_OCAC_coding_ldclumping_08_less_threshold_500_pleiotropy.txt (if any genes with pleiotropic signals found)
  - output_GCPBayes_BCAC_All_OCAC_coding_ldclumping_08_less_threshold_500_errors.txt (if any error occurred during calculation for each gene)

***(for genes with number of SNPs greater than a threshold)***
**Script:** E2_code_gcpbayes_greater_extra_info.R
<br>
**Input**
- Input files needed for running this step:
  - Matrices_output_pipeline_BCAC_ALL_OCAC_coding_clumping_08.Rdata (created from Section 1.4.1)
  - Matrices_extra_info_output_pipeline_BCAC_ALL_OCAC_coding_clumping_08.Rdata (created from Section 1.4.1)

**NOTE:** A user MUST modify the **“DEFINITION SECTION”** which includes the following options:
- Path for two GCPBayes input files: (in our example “~/BCAC_OCAC/”)
- GCPBayes input file names: (in our example “Matrices_output_pipeline_BCAC_ALL_OCAC_coding_clumping_08.Rdata” and “Matrices_extra_info_output_pipeline_BCAC_ALL_OCAC_coding_clumping_08.Rdata”)
- Path for the output files: (in our example “~/BCAC_OCAC/”)
- Output file names. In our example are:
  - output_GCPBayes_BCAC_All_OCAC_coding_ldclumping_08_greater_threshold_500_results
  - output_GCPBayes_BCAC_All_OCAC_coding_ldclumping_08_greater_threshold_500_pleiotropy
  - output_GCPBayes_BCAC_All_OCAC_coding_ldclumping_08_greater_threshold_500_errors
- SNP_number_threshold (A threshold for splitting genes based on SNP numbers) (in our example = 500)
- theta_threshold (for considering genes with potential pleiotropic effect) (in our example = 0.5)

**Output**
- Output files created after running this step:
  - output_GCPBayes_BCAC_All_OCAC_coding_ldclumping_08_greater_threshold_500_results.txt
  - output_GCPBayes_BCAC_All_OCAC_coding_ldclumping_08_greater_threshold_500_pleiotropy.txt (if any genes with pleiotropic signals found)
  - output_GCPBayes_BCAC_All_OCAC_coding_ldclumping_08_greater_threshold_500_errors.txt (if any error occurred during calculation for each gene)

## Visualization
### Section A
**Script:** code_analysis_A_checksumstats_BCAC.R
<br>
**Input**
- Input files needed for this step:
  - BCAC_2020_onco_ALL_reformatted.txt (created from Section A)

**Output**
- Output files created from this step:
  - Some Plots compared to known SNPs in the public databases

The script uses a R package recently introduced (called “CheckSumStats”) (Haycock et al., 2021) to  extracts some SNPs from the BCAC GWAS summary statistics input file and compare them with GWAS database (such as GWAS catalog and 1000 Genome). Here is the steps available in the script:
1.  Extracting SNPs list based on EFO ID of Breast Carcinoma
2.  Checking allele frequency (NOTE: We also checked an output if a user selects the “effect” and “non-effect” alleles columns WRONG!)
3.  Checking the effect allele by comparing Z scores in the test and GWAS catalog datasets
4.  Checking the effect allele by comparing effect allele frequency (EAF) between the test dataset and the GWAS catalog
5.  Checking for errors or analytical issues in the summary data
6.  Comparing the expected and reported effect sizes
7.  Plotting the relative bias, i.e. the percentage deviation of the expected from the reported effect size
8.  Checking whether the reported P values correspond to the reported effect sizes in the dataset
9.  Saving all outputs in one file

The output derived after running the script on BCAC GWAS data (created from Section A) is: 

</br>
<kbd> <img src="Images/checksumstats_BCAC_2020_report.png"/> </kbd>
<br></br>

Also, we performed similar checking on OCAC data as well: 
<br>
**Script:** code_analysis_A_checksumstats_OCAC.R
<br>
**Input**
- Input files needed for this step:
  - OCAC_BCAC_2020_onco_ALL_reformatted.txt (created from Section A)
- Output files created from this step:
  - Some Plots compared to known SNPs in the public databases

Here is the steps available in the script:
1.  Extracting SNPs list based on EFO ID of Ovarian Carcinoma
2.  Checking allele frequency 
3.  Checking the effect allele by comparing Z scores in the test and GWAS catalog datasets
4.  Checking the effect allele by comparing effect allele frequency (EAF) between the test dataset and the GWAS catalog
5.  Checking for errors or analytical issues in the summary data
6.  Comparing the expected and reported effect sizes
7.  Plotting the relative bias, i.e. the percentage deviation of the expected from the reported effect size
8.  Checking whether the reported P values correspond to the reported effect sizes in the dataset
9.  Saving all outputs in one file

The output derived after running the script on OCAC GWAS data (created from Section A) is:

</br>
<kbd> <img src="Images/checksumstats_OCAC_BCAC2020_report.png"/> </kbd>
<br></br>

### Section B
**Script:** code_analysis_B_annotation_file_coding.R
<br>
**Input**
- Input files needed for this step:
  - annot_gencode_v38lift37_modified_gene_class_coding.txt (created from Section B)
- Output files created from this step:
  - Pie Chart

A user could have an overview from the annotation file created in the Section B. 
<br>
There are two possibilities for the visualization of the file:
1.	Running “code_analysis_B_annotation_file_coding.R” script using RStudio. After running, an interactive pie chart would be created in which a user could find various information regarding each part of the plot (total number of genes, percentage, etc.) by moving mouse cursor on the section: 

</br>
<kbd> <img src="Images/pie_chart_interactive_1_1.JPG"/> </kbd>
<br></br>

In addition, by clicking on every items of the legend, it is possible to show/hide that part on the pie chart:

</br>
<kbd> <img src="Images/pie_chart_interactive_2_1.JPG"/> </kbd>
<br></br>

2.	Running “code_analysis_B_annotation_file_coding.R” script using R command line (in UNIX or other Operating Systems (OS)). A pie chart including an overview of genes distribution (here just coding-genes) among all chromosomes would be created in “.jpeg” and “.pdf” format. **NOTE:** since “plotly” package creates an interactive plot, for saving such plots in static file formats (like “.jpeg” or “.pdf”), a user needs to install “orca” on the OS (see https://github.com/plotly/orca#installation for more details).

An example of a pie chart created from the protein-coding genes annotation file (“annot_gencode_v38lift37_modified_gene_class_coding.txt”) is shown:

</br>
<kbd> <img src="Images/pie_chart_protein_coding_genes_chromosome.JPG"/> </kbd>
<br></br>

### Section C
**Script:** code_analysis_C_PLACO_results_one_pair.R 
<br>
**Input**
- Input files needed for this step:
  - output_PLACO_BCAC_2020_ALL_OCAC.txt (created from Section C)
- Output files created from this step:
  - Manhattan Plot 
  - Significant SNPs from the PLACO result (p < 5×10-8)

A user could run the “code_analysis_C_PLACO_results_one_pair.R” script which perform three analyses as follow:
1.	Performing a global correlation analysis between two traits (using Z columns) based on two methods (spearman and pearson) and save the results in “…cor_test.txt” file.
2.	Finding significant SNPs from the PLACO result (p < 5×10-8) and save the results in “output_PLACO_BCAC_2020_ALL_OCAC_Sig.txt” file.
3.	Creation of a Manhattan Plot for all chromosomes and drawing a horizontal dashed line in black (indicating p-value threshold 5×10-8). In addition, SNPs with positive effect are colored in red while SNPs with negative effect are in blue. The plot would be saved as a “.png” file format. An example of the Manhattan Plot for PLACO output for BCAC and OCAC GWAS data (used in our study) is shown:

</br>
<kbd> <img src="Images/output_PLACO_BCAC_2020_ALL_OCAC_Manhattan_Plot.png"/> </kbd>
<br></br>

### Section D
**Script:** code_analysis_D_gcpbayes_input_ggplot.R
<br>
**Input**
- Input files needed for this step:
  - Matrices_output_pipeline_BCAC_ALL_OCAC_coding_withoutclumping.Rdata (created from Section D without LD clumping)
  - *OR*
  - Matrices_output_pipeline_BCAC_ALL_OCAC_coding_clumping_08.Rdata (created from Section D with LD clumping)
- Output files created from this step:
  - Histogram and Bar Plots

A user could run the “code_analysis_D_gcpbayes_input_ggplot.R” script and obtain an overview about distribution of genes with different number of SNPs. We designed two different overviews for a user as follow:
- A histogram (based on number of SNPs available for each gene) for all genes. For example, the histogram for GCPBayes input file for BCAC and OCAC GWAS data (used in our study) while the pipeline run without LD clumping step is shown:

</br>
<kbd> <img src="Images/gcpbayes_input_1.JPG"/> </kbd>
<br></br>

- Based on our experience, in a GWAS data, most of the genes contain a number of SNPs less than 50 (due to our practical experience with various GWAS datasets, especially when working with only protein-coding genes). Therefore, in order to have a deeper overview, a user could define a threshold value (Default = 50) to have a plot with a more detailed view (bar plot) to see the number of genes with SNPs less than the threshold value, as well as a histogram for genes in which the number of SNPs are more than the threshold value. For instance, two plots for GCPBayes input file for BCAC and OCAC GWAS data (used in our study) while the pipeline run without LD clumping step are shown:

</br>
<kbd> <img src="Images/gcpbayes_input_2.JPG"/> </kbd>
<br></br>

### Section E
**Script:** code_analysis_E_gcpbayes_output_karyotype.R
<br>
**Input**
- Input files needed for this step:
  - gene_list.txt (list of pleiotropic genes symbols in one column extracted from the output of the Section E) 
- Output files created from this step:
  - Pleiotropic Genes Karyotype Plot

A user could run the “code_analysis_E_gcpbayes_output_karyotype.R” script which reads the GCPBayes output contains pleiotropic genes information, then show them based on their positions on a chromosome overview (called Karyotype plot). A Karyotype plot for GCPBayes output of BCAC and OCAC GWAS data (used in our study) while the pipeline run without LD clumping step would be as follow:

</br>
<kbd> <img src="Images/gcpbayes_output_karyotype_both.jpg"/> </kbd>
<br></br>

**NOTE:** Since for drawing a Karyotype plot, the script uses the HUGO Gene Nomenclature Committee (HGNC) symbol (for each input gene symbol) extracted from Ensembl database, it would **NOT** be able to show the genes that their HGNC symbol did not find.  
<br>
<br>
**Script:** code_analysis_E_gcpbayes_output_table_overview.R
<br>
**Input**
- Input files needed for this step:
  - gene_list.txt (list of pleiotropic genes symbols with “chr”, “gene_length”, and “snp_number” columns extracted from the output of the Section E) 
- Output files created from this step:
  - gcpbayes_pleiotropy_summary_table.txt (Similar to Table 1 in the Manuscript)

There is also a script called “code_analysis_E_gcpbayes_output_table_overview.R” which a user could perform on the GCPBayes pleiotropic genes output and has an overall information about potential pleiotropic genes. Here is an example:

| chr	| #genes	| Min gene length (bp)	| Max gene length (bp)	| Min SNP numbers	| Max SNP numbers |
| -- | -- | -- | -- | -- | -- |
| 1	| 5	| 4,959	| 57,942	| 3	| 98 |
| 2	| 3	| 54,268	| 160,447	| 82	| 214 |
| 3	| 11	| 33,195	| 363,895	| 18	| 524 |
| 4	| 9	| 19,285	| 148,512	| 32	| 549 |
| 5	| 6	| 16,272	| 96,657	| 37	| 472 |
| 6	| 8	| 40,931	| 353,488	| 55	| 696 |
| 7	| 3	| 37,608	| 261,613	| 71	| 278 |
| 8	| 3	| 97,883	| 158,920	| 178	| 447 |
| 9	| 5	| 10,883	| 226,958	| 18	| 320 |
| 10	| 10	| 12,893	| 354,629	| 3	| 446 |
| 11	| 7	| 11,010	| 255,035	| 11	| 530 |
| 12	| 14	| 5,150	| 328,502	| 2	| 392 |
| 15	| 9	| 8,255	| 225,330	| 18	| 364 |
| 16	| 8	| 11,508	| 88,046	| 12	| 197 |
| 17	| 21	| 444	| 371,278	| 1	| 609 |
| 19	| 14	| 3,130	| 123,106	| 1	| 143 |
| 20	| 7	| 6,807	| 178,396	| 10	| 162 |
| 21	| 2	| 49,232	| 104,701	| 82	| 151 |
| 22	| 3	| 15,476	| 701,851	| 30	| 652 |

## References
- Frankish,A. et al. (2019) GENCODE reference annotation for the human and mouse genomes. Nucleic Acids Res., 47, D766–D773.
- Haycock,P.C. et al. (2021) Design and quality control of large-scale two-sample Mendelian randomisation studies. medRxiv, 2021.07.30.21260578.
- Phelan,C.M. et al. (2017) Identification of 12 new susceptibility loci for different histotypes of epithelial ovarian cancer. Nat. Genet., 49, 680–691.
- Zhang,H. et al. (2020) Genome-wide association study identifies 32 novel breast cancer susceptibility loci from overall and subtype-specific analyses. Nat. Genet. 2020 526, 52, 572–581.

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

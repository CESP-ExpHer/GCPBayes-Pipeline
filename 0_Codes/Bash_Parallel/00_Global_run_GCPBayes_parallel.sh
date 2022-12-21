#===============================================================================
#           USAGE:  ./00_Global_run_GCPBayes_parallel.sh parameters_parallel.ini readinputs.txt
#
#     DESCRIPTION:  global codes to run the GCPBayes pipeline as parallel: this code encapsulate all the different procedures in one running code, 
#                                                               with all the parameters to be specified by the user in a separated .ini file
#                                                               and run GCPBayes in parallel which means multiple groups run in each for loop
#
#===============================================================================
#          AUTHOR:  Pierre-Emmanuel Sugier
#         VERSION:  Yazdan Asgari
#         CREATED:  2021-11-23
#        REVISION:  2022-12-21
# ================================================================================                 
#          DIFFERENT STEPS:
#                    
#                    Step 1 : To extract common SNPs between dataset
#                    Step 2 : To run PLACO on each dataset
#                    Step 3 : To perform LD Clumping (locally)
#                    Step 4 : To prepare the GCPBayes inputs in the right format
#                    Step 5 : To get correct list of groups for further analysis according to wanted threshold (length of groups)
#                    Step 6 : To run GCPBayes(DS) as Parallel
#                    
# ================================================================================                 
#    REQUIREMENTS:  R
# 	 REQUIRED LIBRARIES
# please install the following packages if needed
#	 BhGLM, BioCircos, BiocManager, biomaRt, CheckSumStats, data.table, datasets
#	 devtools, GCPBayes, genetics.binaRies, ggpubr, gridExtra, gwasrapidd, ieugwasr
#	 karyoploteR, MASS, optparse, patchwork, PLACO, plotly, readxl, regioneR
#	 shiny, splitstackshape, tictoc, tidyverse, vroom, doParallel, foreach, plyr, doSNOW
# ================================================================================

## REQUIRED INFORMATIONS
##########################################################################################################################################
start=$(date +%s.%N)

ref_file_parameters="$1"

source $ref_file_parameters

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##       START RUN       ##
###########################

# STEP 1 - To extract common SNPs between dataset
#####################################################

if $toclump ; then
Rscript ${script_dir}/C1_code_find_common_snps_one_pair.R \
                    --path ${work_dir} \
                    --inputs ${Inputs_file_to_read} \
                    --out ${output_step1} \
                    --info ${info_threshold} \
                    --maf ${MAF_threshold}
fi

# Step 2 - To run PLACO on each dataset
#####################################################

if  $toclump ; then
  if  -e ${work_dir}/${output_step2}.txt  # this step does not run if a corresponding output file already exist
  then
    echo "A corresponding PLACO output is present in the work directory and is going to be used."
  else
    Rscript ${script_dir}/C2_code_run_PLACO_decor_one_pair.R \
        --path ${work_dir} \
        --inputs ${Inputs_file_to_read} \
        --prefix_in ${output_step1} \
        --prefix_out ${output_step2} \
        --pval ${placo_pval_threshold}
  fi
fi

# Step 3 - To perform LD Clumping (locally)
#####################################################

if  $toclump ; then
  Rscript ${script_dir}/C3_code_ldclumping_local.R \
        --path ${work_dir} \
        --file ${output_step2} \
        --out ${output_step3}_${clump_threshold_r2}.txt \
        --ref ${ref_path_b_files} \
        --r2 ${clump_threshold_r2} \
        --kb ${clump_threshold_kb} \
        --p ${clump_threshold_p}
fi

# Step 4 - To prepare the GCPBayes inputs in the right format
#####################################################

Rscript ${script_dir}/D1_code_pipeline_annot_coding_withoutldclumping_extra_info.R \
        --pathin ${work_dir} \
        --inputs ${Inputs_file_to_read} \
        --annot_gwas_path ${path_gwas_annot} \
        --annot_gwas_file ${file_gwas_annot} \
        --annot_path ${path_annot} \
        --annot_file ${file_annot} \
        --info ${info_threshold} \
        --maf ${MAF_threshold} \
        --pathout ${work_dir} \
        --out ${output_step4wc}

if $toclump ; then
  Rscript ${script_dir}/D2_code_pipeline_annot_coding_ldclumping_extra_info.R \
        --pathin ${work_dir} \
        --inputs ${Inputs_file_to_read} \
        --annot_gwas_path ${path_gwas_annot} \
        --annot_gwas_file ${file_gwas_annot} \
        --annot_path ${path_annot} \
        --annot_file ${file_annot} \
        --ldclumpout ${output_step3}_${clump_threshold_r2}.txt \
        --info ${info_threshold} \
        --maf ${MAF_threshold} \
        --pathout ${work_dir} \
        --out ${output_step4c}
fi

# Step 5 - To separate list of groups for GCPBayes analysis according to selected threshold (length of groups)
#####################################################

if $toclump ; then
Rscript ${script_dir}/D3_code_separate_groups_length_threshold.R \
        --path ${work_dir} \
        --file_noclump ${output_step4wc} \
        --file_clump ${output_step4c} \
        --out_noclump ${output_step5wc} \
        --out_clump ${output_step5c} \
        --t_SNP_number ${group_clump_threshold}
		else
Rscript ${script_dir}/D3_code_separate_groups_length_threshold_noclump.R \
        --path ${work_dir} \
        --file_noclump ${output_step4wc} \
        --file_clump ${output_step4c} \
        --out_noclump ${output_step5wc} \
        --out_clump ${output_step5c} \
        --t_SNP_number ${group_clump_threshold}		
fi               

# Step 6 - To run GCPBayes (DS)
#####################################################

# A - To run GCPBayes (DS) on groups without ld clumping (length <= group_clump_threshold)
Rscript ${script_dir}/E1_code_gcpbayes_less_extra_info_parallel.R \
        --path ${work_dir} \
        --file ${output_step5wc} \
        --t_SNP_number ${group_absolute_threshold} \
        --t_theta ${theta_exploration} \
        --pathout ${output_dir} \
        --out ${output_step6wc} \
		--cpu ${cpu_numbers}
# WARNINGS : theta_exploration is used (not theta_plot)

# B - To run GCPBayes (DS) on groups after ld clumping (group_clump_threshold < length <= group_absolute_threshold)
if $toclump ; then
Rscript ${script_dir}/E1_code_gcpbayes_less_extra_info_parallel.R \
        --path ${work_dir} \
        --file ${output_step5c} \
        --t_SNP_number ${group_absolute_threshold} \
        --t_theta ${theta_exploration} \
        --pathout ${output_dir} \
        --out ${output_step6c} \
		--cpu ${cpu_numbers}
# WARNINGS : theta_exploration is used (not theta_plot)
fi

#####################################################

dur=$(echo "$(date +%s.%N) - $start" | bc)

printf $dur




































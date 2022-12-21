# 					CODE FOR REFORMATTING THE GWAS DATA (BCAC_2020_onco_ALL)
# ================================================================================
# Summary: The code reads the GWAS Summary Statistics data for BCAC_2020_onco_all
# "icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt"
# and writes it in a proper format that could be used for our future analyses
# ================================================================================
# Written first by: Elise Lucotte
# Modified by: Yazdan Asgari
# Initial Creation Date: 12/2020
# Edited Date: 1/2022
# https://cesp.inserm.fr/en/equipe/exposome-and-heredity
# ================================================================================

# libraries used in the code
from collections import defaultdict
from scipy import stats
import datetime
# ================================================================================

# for calculation the running time of the program
begin_time = datetime.datetime.now()

# ================================================================================
# 							DEFINITION SECTION
# should be changed by a user
# ================================================================================
# directory in which your input data exist
data_dir='~/BCAC_OCAC/'

# the file names of your input data
file1='icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt'
names=['ALL']

# directory in which your output data should be written
output_dir='~/BCAC_OCAC/'
# ================================================================================

# dictionnary with the complement of nucleotides
complement=defaultdict(lambda:str)
complement['A']='T'
complement['T']='A'
complement['C']='G'
complement['G']='C'

# the file names of your output data 
fw_data_reformatted=open(f'{output_dir}BCAC_2020_onco_{names[0]}_reformatted.txt','w')
fw_snp_ambiguous=open(f'{output_dir}BCAC_2020_onco_{names[0]}_SNP_ambiguous.txt','w')
fw_snp_duplicated=open(f'{output_dir}BCAC_2020_onco_{names[0]}_SNP_dupicated.txt','w')
fw_snp_duplicated_set=open(f'{output_dir}BCAC_2020_onco_{names[0]}_SNP_duplication_set.txt','w')
fw_snp_weird_alleles=open(f'{output_dir}BCAC_2020_onco_{names[0]}_SNP_weird_alleles.txt','w')

# writing the columns names in the output file ("fw_data_reformatted")
print(' '.join(['snp','\t','chr','\t','bp_hg19','\t','Effect_A','\t','nonEffect_A','\t','beta','\t','se','\t','pval','\t','info','\t','ngt','\t','EAF','\t','MAF']),file=fw_data_reformatted)

# =================================================================================
# Reading the file for the FIRST time to create a set with the duplicated SNPs
# =================================================================================
# if two SNPs have the same rsid or the same CHR:POS, they will be removed later in the code
fr_gwas_input=open('{0}{1}'.format(data_dir, file1),'r')
header=fr_gwas_input.readline()
set_dup=set()
set_poschr=set()
set_SNP=set()

# in the "icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt" file:
# column 24: SNP_ID (rs... or chr:pos:Baseline_A:Effect_A) (Baseline_A=non-Effect Allele, Effect_A=Effect Allele)
# column 25: chr
# column 26: position
for line in fr_gwas_input: 
	line=str(line)
	line=line.split()
	ID=line[23]
	if ID[0:2]=='rs':
		ID=ID.split(':')[0]
	CHR=line[24]
	POS=line[25]
	if ID in set_SNP or CHR+':'+POS in set_poschr:
		if ID[0:2]=='rs':
			set_dup.add(ID)
			set_dup.add(CHR+':'+POS)
		else :
			set_dup.add(CHR+':'+POS)
	set_SNP.add(ID)
	set_poschr.add(CHR+':'+POS)
fr_gwas_input.close()

# End of First Reading
print('End of FIRST Reading Step')

# =================================================================================
# Reading the file for a SECOND time
# =================================================================================
fr_gwas_input=open('{0}{1}'.format(data_dir, file1),'r')
header=fr_gwas_input.readline()

# Keeping track of the different categories of SNPs:
# 1-SNPs which kept, 
# 2-SNPs with a low value for info (or r2),
# 3-SNPs which are ambiguous, 
# 4-SNPs with NULL values for beta/se/info/pval,
# 5-SNPs with zero beta&se,
# 6-SNPs which are duplicated,# SNPs include an allele that are not composed of ATGC,
# 7-SNPs with an allele that are not composed of ATGC,
# 8-SNPs with a length > 1 and with alleles that are not composed of ATGC	 
nb_SNP_kept=0
nb_SNP_info_low=0
nb_SNP_ambiguous=0
nb_NULL=0
nb_ZERO=0
nb_SNP_duplicated=0
nb_SNP_weird_alleles=0
nb_SNP_allele_size_more_than_one=0

dico_SNPs=defaultdict(list)

# in the "icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt" file:
# column 24: SNP_ID (rs... or chr:pos:Baseline_A:Effect_A) (Baseline_A=non-Effect Allele, Effect_A=Effect Allele)
# column 25: chr
# column 26: position
# column 27: Effect_Allele
# column 28: non-Effect_Allele
for line in fr_gwas_input: 
	line=str(line)
	line=line.split()
	ID=line[23]
	if ID[0:2]=='rs':
		ID=ID.split(':')[0]
	CHR=line[24]
	POS=line[25]
	A1=line[26]
	A2=line[27]

	# column 29 (line[28]): eaf = EAF_control (we considered the frequency in controls)
	# column 31 (line[30]): info (or r2)
	# column 33 (line[32]): beta
	# column 34 (line[33]): se				
	# column 38 (line[37]): P_value (P1df_risk_chi)			

	# if (eaf or info or beta or se or P_value) == NULL or NA, it removes the SNP (means it skips the line)
	if line[28]=='NULL' or line[28]=='NA' or line[30]=='NULL' or line[30]=='NA' or line[32]=='NULL' or line[32]=='NA' or line[33]=='NULL' or line[33]=='NA' or line[37]=='NULL' or line[37]=='NA':
		nb_NULL+=1
		continue

	eaf=line[28]
	info=line[30]
	beta=float(line[32])
	se=float(line[33])
	p=line[37]

	# if both beta and se == 0, it removes the SNP (means it skips the line)
	if beta==0.0 and se==0.0:
		nb_ZERO+=1
		continue	
	
	# if one of the alleles has a length not equal to 1,
	# and the SNP includes alleles that are not composed of ATGC, it will be removed 	
	if len(A1)!=1:
		alt_A1=''
		for i in range(0,len(A1)):
			try:
				alt_A1=alt_A1+complement[A1[i]]
			except TypeError:
				nb_SNP_weird_alleles+=1
				print('\t'.join([ID, CHR, POS, str(eaf), str(beta), str(se), str(p), str(info), A1, A2]), file=fw_snp_weird_alleles)
				next
	else:
		alt_A1=complement[A1]
	if len(A2)!=1:
		alt_A2=''
		for i in range(0,len(A2)):
			try:
				alt_A2=alt_A2+complement[A2[i]]
			except TypeError:
				nb_SNP_weird_alleles+=1
				print('\t'.join([ID, CHR, POS, str(eaf), str(beta), str(se), str(p), str(info), A1, A2]), file=fw_snp_weird_alleles)
				next
	else:
		alt_A2=complement[A2]
	
	# if the SNP is an ambiguous SNPs, it removes the SNP (means it skips the line) and prints the SNP in its output file
	if len(A1)==1 and len(A2)==1 and A1==complement[A2]:
		# writing the SNP information in the output file ("fw_snp_ambiguous")
		print('\t'.join([CHR, POS, A1, A2]), file=fw_snp_ambiguous)
		nb_SNP_ambiguous+=1
		continue

	# if A1 and A2 have more than one nucleotide and have the same length, 
	# it checks if they are not ambiguous
	# If so, it removes the SNP (means it skips the line)
	if len(A1)!=1 and len(A2)!=1 and len(A1)==len(A2):
		flag=0
		for i in range(0, len(A1)):
			if A1[i]==complement[A2[i]] :
				flag+=1
		if flag == len(A1):
			nb_SNP_allele_size_more_than_one+=1
			continue
			
	# if the SNP is in the duplicated set, it removes the SNP (means it skips the line) and prints the SNP in its output file
	if ID in set_dup or CHR+':'+POS in set_dup:
		# writing the duplicated SNP information in the output file ("fw_snp_duplicated")
		print('\t'.join([ID, CHR, POS, A1, A2, str(beta), str(se), str(p), info]), file=fw_snp_duplicated)
		nb_SNP_duplicated+=1
		continue

	# To get the MAF as Minor Allele Frequency based on EAF_control
	if float(eaf)<0.5:
		MAF=eaf
	else:
		MAF=str(1-float(eaf))
		
	# counting the SNPs with a value smaller than 0.9 for info (info means: r2/column 31/line[30])
	if float(info)<0.9:
		nb_SNP_info_low+=1

	# writing the line in the output file ("fw_data_reformatted")
	new_line=[ID, CHR, POS, A1, A2, str(beta), str(se), str(p), info, '0', eaf, MAF ]
	print ('\t'.join(new_line),file=fw_data_reformatted)
	nb_SNP_kept+=1
	
# End of Second Reading
print('End of SECOND Reading Step')

# =================================================================================
# calculation of some summaries of what the code did on the data
dico_SNPs['nb_SNP_kept'].append(str(nb_SNP_kept))
dico_SNPs['nb_SNP_info_low'].append(str(nb_SNP_info_low))
dico_SNPs['nb_SNP_ambiguous'].append(str(nb_SNP_ambiguous))
dico_SNPs['nb_SNP_NULL/NA'].append(str(nb_NULL))
dico_SNPs['nb_SNP_beta_se_ZERO'].append(str(nb_ZERO))
dico_SNPs['nb_SNP_duplicated'].append(str(nb_SNP_duplicated))
dico_SNPs['nb_SNP_duplicated_set'].append(str(len(set_dup)))
dico_SNPs['nb_SNP_weird_alleles'].append(str(nb_SNP_weird_alleles))
dico_SNPs['nb_SNP_allele_size_more_than_one'].append(str(nb_SNP_allele_size_more_than_one))

# create an output file for writing the summaries ("dico_SNPs")
fw_summary_SNP_All=open(f'{output_dir}BCAC_2020_onco_{names[0]}_Summary.txt','w')
print('\t'.join(['type']+names), file=fw_summary_SNP_All)
for key in dico_SNPs.keys():
	# writing the summaries in the output file ("w_summary_SNP_All")
	print('\t'.join([key]+dico_SNPs[key]), file=fw_summary_SNP_All)

# writing the duplicated SNPs set found in the GWAS data and used for SNP duplication process ("fw_snp_duplicated_set")
print("\n".join(set_dup), file=fw_snp_duplicated_set)
	
# closing all opened files 
fw_data_reformatted.close()
fw_snp_ambiguous.close()
fw_snp_duplicated.close()
fw_snp_duplicated_set.close()
fw_snp_weird_alleles.close()
fw_summary_SNP_All.close()
fr_gwas_input.close()

print('Reached to the END of program')

# ================================================================================
# printing the running time of the program 
print(datetime.datetime.now() - begin_time) 
# ================================================================================

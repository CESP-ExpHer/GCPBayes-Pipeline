# 			CODE FOR REFORMATTING AND STANDARDIZATION OF THE GWAS DATA 
#					(OCAC based on BCAC_2020_onco_ALL) 
# ================================================================================
# Summary: The code reads the GWAS Summary Statistics data for OCAC
# "extraction_OCAC.txt"
# and writes it in a proper format that could be used for our future analyses
# It also checks the "Effect" and "non-Effect" Alleles to be in accordance with 
# the reference file (here "BCAC_2020_onco_ALL_reformatted.txt" file)
# ================================================================================
# Written first by: Elise Lucotte
# Modified by: Yazdan Asgari
# Initial Creation Date: 12/2020
# Edited Date: 1/2022
# https://cesp.inserm.fr/en/equipe/exposome-and-heredity
# ================================================================================
# libraries used in the code
from collections import defaultdict
import datetime

# ================================================================================
# for calculation the running time of the program
begin_time = datetime.datetime.now()

# ================================================================================
# 							DEFINITION SECTION
# should be changed by a user
# ================================================================================
# directory in which your reference data exist
ref_dir='~/BCAC_OCAC/'

# Determining the reference data
file_ref='BCAC_2020_onco_ALL_reformatted.txt'

# directory in which your input data exist
input_dir='~/BCAC_OCAC/'

# Determining the input data
file_gwas='extraction_OCAC.txt'

# directory in which your output data should be written
output_dir='~/BCAC_OCAC/'

# Determining the names used as first and second parts of the output files names
outputfile_part1=f'OCAC'
outputfile_part2=f'_BCAC_2020_onco_ALL'

# ================================================================================
# 				            RUNNING SECTIONS
# ================================================================================
# Reading the reference file for comparing SNPs with (A1 and A2)
# here the reference file is "BCAC_2020_onco_ALL"
# --------------------------------------------------------------------------------
d={}
f_SNP=open('{0}{1}'.format(ref_dir, file_ref), 'r')
for line in f_SNP:
	line=line.split()
    # puting the SNPs information in a dictionnary as d[CHR][BP]=[SNP, A1, A2]
	# CHR : BP : [SNP, A1, A2]
	try:
		d[line[1]][line[2]]=[line[0],line[3],line[4]]
	except KeyError:
		d[line[1]]={}
		d[line[1]][line[2]]=[line[0],line[3],line[4]]
f_SNP.close()

# printing the keys of the dictionary (only for chromosome numbers)
print(d.keys())

# creation of a dictionnary with the complement of nucleotides
complement=defaultdict(lambda:str)
complement['A']='T'
complement['T']='A'
complement['C']='G'
complement['G']='C'

# --------------------------------------------------------------------------------
# defining a dictionnary to count the different categories of SNPs
dico_SNPs=defaultdict(list)

# --------------------------------------------------------------------------------
# Start of Reformatting and Standardization for each OCAC file
# --------------------------------------------------------------------------------
    
# the file names of your output data 
fw_SNP_removed_from_ref=open(f'{output_dir}{outputfile_part1}{outputfile_part2}_SNP_removed_from_ref.txt','w')
fw_SNP_removed_from_file=open(f'{output_dir}{outputfile_part1}{outputfile_part2}_SNP_removed_from_file.txt','w')
fw_snp_ambiguous=open(f'{output_dir}{outputfile_part1}{outputfile_part2}_SNP_ambiguous.txt','w')
fw_snp_duplicated=open(f'{output_dir}{outputfile_part1}{outputfile_part2}_SNP_dupicated.txt','w')
fw_snp_duplicated_set=open(f'{output_dir}{outputfile_part1}{outputfile_part2}_SNP_duplication_set.txt','w')
fw_snp_weird_alleles=open(f'{output_dir}{outputfile_part1}{outputfile_part2}_SNP_weird_alleles.txt','w')

# writing the columns names in the output file ("fw_data_reformatted")
fw_data_reformatted=open(f'{output_dir}{outputfile_part1}{outputfile_part2}_reformatted.txt','w')
print(' '.join(['snp','\t','chr','\t','bp_hg19','\t','Effect_A','\t','nonEffect_A','\t','beta','\t','se','\t','pval','\t','info','\t','ngt','\t','EAF','\t','nEAF','\t','MAF']),file=fw_data_reformatted)

# Keeping track of the different categories of SNPs:
# 1-SNPs which kept, 
# 2-SNPs which are ambiguous, 
# 3-SNPs with NULL values for freq/beta/se/pval,
# 4-SNPs with zero beta&se,
# 5-SNPs which are duplicated,
# 6-SNPs include an allele that are not composed of ATGC,
# 7-SNPs with allele size more than one, 
# 8-SNPs which removed due to mismatches of both alleles between the file and reference files, 
# 9-SNPs which removed due to not be available in the reference file
nb_SNP_kept=0
nb_SNP_ambiguous=0
nb_NULL=0
nb_ZERO=0
nb_SNP_duplicated=0
nb_SNP_weird_alleles=0
nb_SNP_allele_size_more_than_one=0
nb_SNP_removed_from_file=0
nb_SNP_removed_from_ref=0

# =================================================================================
# Reading the file for the FIRST time to create a set with the duplicated SNPs
# =================================================================================
# if two SNPs have the same CHR:POS, they will be removed later in the code
fr_gwas_input=open('{0}{1}'.format(input_dir, file_gwas), 'r')
fr_gwas_input.readline()
set_dup=set()

# in the OCAC file:
# column 2: chr
# column 3: position
# column 4: A1 (Effect Allele named A1 in OCAC Data)
# column 5: A2 (non-Effect Allele named A2 in OCAC Data)
for line in fr_gwas_input: 
	line=line.strip()
	line=line.split()
	CHR= line[1]
	POS=line[2]
	A1= line[3]
	A2= line[4]
	if CHR+':'+POS in set_dup:
		set_dup.add(CHR+':'+POS)

print('nb SNP duplicated:',len(set_dup))

# writing the duplicated SNPs set found in the GWAS data and used for SNP duplication process ("fw_snp_duplicated_set")
print("\n".join(set_dup), file=fw_snp_duplicated_set)
    
fr_gwas_input.close()

# End of First Reading
print('End of FIRST Reading Step')
    
# =================================================================================
# Reading the file for a SECOND time
# =================================================================================
fr_gwas_input=open('{0}{1}'.format(input_dir, file_gwas), 'r')
fr_gwas_input.readline()
    
for line in fr_gwas_input: 
	# in the OCAC file:
    # -----------------------------------------------------------------------------
    # column 2 (line[1]): chr
    # column 3 (line[2]): position
    # column 4: A1 (Effect Allele named A1 in OCAC Data)
    # column 5: A2 (non-Effect Allele named A2 in OCAC Data)
    # column 6 (line[5]): Effect Allele Frequency (EAF named freq1 in OCAC Data)
    # column 7 (line[6]): non-Effect Allele Frequency (nEAF named freq2 in OCAC Data) 
    # column 8 (line[7]): beta
    # column 9 (line[8]): se
    # column 10 (line[9]): p-value
    # column 11 (line[10]): info (imputation quality score) 
    # column 12 (line[11]): Number of samples
    # -----------------------------------------------------------------------------		
	line=line.strip()
	line=line.split()
	CHR= line[1]
	POS=line[2]
	A1= line[3]
	A2= line[4]
    
    # if (freq1 or freq2 or beta or se or p) == NULL or NA, it removes the SNP (means it skips the line)
	if line[5] == 'NA' or line[5]=='NULL' or line[6] == 'NA' or line[6]=='NULL' or line[7] == 'NA' or line[7]=='NULL' or line[8] == 'NA' or line[8]=='NULL' or line[9] == 'NA' or line[9]=='NULL':
		nb_NULL+=1
		continue
    
    # Effect Allele Frequency    
	eaf=float(line[5])
    # non-Effect Allele Frequency
	neaf=float(line[6])
    # Minor Allele Frequency (MAF)
	maf=min(eaf,neaf)
	beta=float(line[7])
	se=line[8]
	p=line[9]
	info=line[10]
	N=line[11]
    
#   if both beta and se == 0, we continue (means we skip the line)
	if beta==0.0 and se==0.0:
		nb_ZERO+=1
		continue

	#--------------------------------------------------------------------------------
	# if the SNP is in the duplicated set, it removes the SNP (means it skips the line) and prints the SNP in its output file
	if CHR+':'+POS in set_dup:
        # writing the duplicated SNP information in the output file ("fw_snp_duplicated")
		print('\t'.join([CHR, POS, A1, A2, str(beta), se, p, info]), file=fw_snp_duplicated)
		nb_SNP_duplicated+=1
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
    
    # -----------------------------------------------------------------------------
    # 			Standardization (Alleles Checking) Steps 
    # -----------------------------------------------------------------------------    
	try:
		# 0) checking the reference file to see if it exists (means if CHR:POS exists in the reference file)
		ref_A1=d[CHR][POS][1]
		ref_A2=d[CHR][POS][2]

        # 1) if one of the alleles has a length not equal to 1,
        # and the SNP includes alleles that are not composed of ATGC, it will be removed 
		if len(A1)!=1:
			alt_A1=''
			for i in range(0,len(A1)):
				try:
					alt_A1=alt_A1+complement[A1[i]]
				except TypeError:
					nb_SNP_weird_alleles+=1
					print('\t'.join([CHR, POS, A1, A2, str(beta), se, p, info, d[CHR][POS][0], d[CHR][POS][1], d[CHR][POS][2]]), file=fw_snp_weird_alleles)						
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
					print('\t'.join([CHR, POS, A1, A2, str(beta), se, p, info, d[CHR][POS][0], d[CHR][POS][1], d[CHR][POS][2]]), file=fw_snp_weird_alleles)
					next
		else:
			alt_A2=complement[A2]

		# 2) if the SNP is an ambiguous SNPs, it removes the SNP (means it skips the line) and prints the SNP in its output file
		if A1==alt_A2:
            # writing the SNP information in the output file ("fw_snp_ambiguous")
			print('\t'.join([CHR, POS, A1, A2, d[CHR][POS][0], d[CHR][POS][1], d[CHR][POS][2]]), file=fw_snp_ambiguous)
			nb_SNP_ambiguous+=1
			continue

		# 3) if the alleles of the file and the reference alleles are the same
        # storing the A1 and A2 from the reference file
		elif A1 == ref_A1 and A2 == ref_A2:
			new_line=[d[CHR][POS][0], CHR, POS, ref_A1, ref_A2, str(beta), se, p, info, '0', str(eaf), str(neaf), str(maf) ]
            # writing the line in the output file ("fw_data_reformatted")
			print('\t'.join(new_line),file=fw_data_reformatted)
			nb_SNP_kept+=1

		# 4) if A1 and A2 are inverse from the reference, 
        # storing the A1 and A2 from the reference file in the right order
        # it changes the sign of the beta 
		elif A1 == ref_A2 and A2 == ref_A1:
			new_line=[d[CHR][POS][0], CHR, POS, ref_A2, ref_A1, str(-beta), se, p, info, '0', str(eaf), str(neaf), str(maf) ]
            # writing the line in the output file ("fw_data_reformatted")
			print('\t'.join(new_line),file=fw_data_reformatted)
			nb_SNP_kept+=1

			
		# 5) if the alleles of the file are the complements of the reference alleles
        # storing the A1 and A2 from the reference file
		elif alt_A1==ref_A1 and alt_A2==ref_A2:
			new_line=[d[CHR][POS][0], CHR, POS, ref_A1, ref_A2, str(beta), se, p, info, '0', str(eaf), str(neaf), str(maf) ]
            # writing the line in the output file ("fw_data_reformatted")
			print ('\t'.join(new_line),file=fw_data_reformatted)
			nb_SNP_kept+=1

		# 6) if the alleles of the file are the complements of the reference alleles, BUT they are switched
        # storing the A1 and A2 from the reference file
        # it changes the sign of the beta  
		elif alt_A1==ref_A2 and alt_A2==ref_A1:
			new_line=[d[CHR][POS][0], CHR, POS, ref_A2, ref_A1, str(-beta), se, p, info, '0', str(eaf), str(neaf), str(maf) ]
            # writing the line in the output file ("fw_data_reformatted")
			print ('\t'.join(new_line),file=fw_data_reformatted)
			nb_SNP_kept+=1

		# If the SNP does not fit to any of the above categories, it removes the SNP (means it skips the line)
		else:
            # writing the line in which a SNP does not enter any of the above categories ("fw_SNP_removed_from_file")
			print('\t'.join([CHR, POS, A1, A2, d[CHR][POS][0], d[CHR][POS][1], d[CHR][POS][2]]), file=fw_SNP_removed_from_file)
			nb_SNP_removed_from_file+=1
			continue

	# if the SNP does not exist in the reference, it removes the SNP (means it skips the line)
	except KeyError: 
        # writing the SNPs that are not in the reference (fw_SNP_removed_from_ref) 
		print('\t'.join([CHR, POS, A1, A2]), file=fw_SNP_removed_from_ref)
		nb_SNP_removed_from_ref+=1
		continue

# End of Second Reading
print('End of SECOND Reading Step')
  
# =================================================================================
# calculation of some summaries of what the code did on the data
dico_SNPs['nb_SNP_kept'].append(str(nb_SNP_kept))
dico_SNPs['nb_SNP_ambiguous'].append(str(nb_SNP_ambiguous))
dico_SNPs['nb_SNP_NULL/NA'].append(str(nb_NULL))
dico_SNPs['nb_SNP_beta_se_ZERO'].append(str(nb_ZERO))
dico_SNPs['nb_SNP_duplicated'].append(str(nb_SNP_duplicated))
dico_SNPs['nb_SNP_weird_alleles'].append(str(nb_SNP_weird_alleles))
dico_SNPs['nb_SNP_allele_size_more_than_one'].append(str(nb_SNP_allele_size_more_than_one))
dico_SNPs['nb_SNP_removed_from_ref'].append(str(nb_SNP_removed_from_ref))
dico_SNPs['nb_SNP_removed_from_file'].append(str(nb_SNP_removed_from_file))	
dico_SNPs['nb_SNP_duplicated_set'].append(str(len(set_dup)))

# create an output file for writing the summaries ("dico_SNPs")
fw_summary_SNP_subtypes=open(f'{output_dir}{outputfile_part1}{outputfile_part2}_Summary_SNP.txt','w')
print('\t'.join(['type']), file=fw_summary_SNP_subtypes)
for key in dico_SNPs.keys():
	# printing the summaries in the output file ("fw_summary_SNP_subtypes")
	print('\t'.join([key]+dico_SNPs[key]), file=fw_summary_SNP_subtypes)

# closing all opened files
fw_data_reformatted.close()
fw_snp_ambiguous.close()
fw_snp_duplicated.close()
fw_snp_duplicated_set.close()
fw_snp_weird_alleles.close()
fw_SNP_removed_from_file.close()
fw_SNP_removed_from_ref.close()
fw_summary_SNP_subtypes.close()
fr_gwas_input.close()

print('Reached to the END of program')

# ================================================================================
# printing the running time of the program 
print(datetime.datetime.now() - begin_time) 
# ================================================================================

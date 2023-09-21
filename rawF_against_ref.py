#!/usr/bin/python3
from __future__ import division
import pandas as pd
import sys

"""This script is used to extract reference genome annotation (.gff3) file for region annotation 
	and for each sample, calculate:
		- calculate overall mean F  in CDS area 
		- calculate overall mean F in non coding area (outside of CDS)
		- calculate overall mean F for gene region 
		- calculate overall mean F for 5kb upstream of gene region (cis regulatory region)
		- calculate overall mean F for gene region + 5kb upstream 
		- calculate overall mean F in whole genome
		(note: overall mean F rate = sum of F rates / # of CPG sites)
	
	Input:
		annotation_handle: the path to the reference genome annotaion
		DamMet raw F files 
	Ouput:
		result text files - one file / samples 
			columns [chromosome,CDS,non_coding_gene,gene,cis,cis+gene,all_cpg_sites]

	Function:
		read_genome(gff_handle)
		read_raw_f(raw_f_handle)
		calculate_CDS_mean_f(full_annotation_df,raw_dict)
		calculate_nc_mean_f(full_annotation_df,raw_df)
		alculate_gene_and_regulatory_mean_f(full_annotation_df,raw_dict,extra_region = 5000)
"""

def read_genome(gff_handle):
	"""
	this function is used to read gff file Asian elephant genome
	input:
		gff_handle: filepath to gff reference file
	output: 
		full_annotation_df: a pandas dataframe includes chromsomes, type, start, end of annotations
	"""
	raw_annotation_df = pd.read_table(gff_handle, comment = "#", skiprows = [0], header = None)
	#print(raw_annotation_df.iloc[[1],:])
	full_annotation_df = raw_annotation_df.iloc[:,[0,2,3,4]]
	new_col_names = {0:"chromosome",2:"type",3:"start",4:"end"}
	full_annotation_df = full_annotation_df.rename(columns = new_col_names)

	return full_annotation_df

def read_raw_f(raw_f_handle):
	"""This function is used to read raw F files as dictionary 
	inpput:
		raw_f_handle: file path to raw F file 
	output:
		raw_dict: dictionary with key value pairs "center"-"f rate", and one "current_chrom"-chromsomes 
	"""
	raw_df = pd.read_table(raw_f_handle, sep = ' ', dtype = 'str', usecols = [0,5])
	
	#get the current chromosome from raw F file 
	current_chrom = raw_df['center'][1000].partition(':')[0]

	#formatting the location in chromosome into integer, mythlation rate as float
	raw_df['center'] = raw_df['center'].str.split(':').str.get(-1).astype(int)
	raw_df['f'] = raw_df['f'].astype(float)

	#turn the dataframe into a dictionary 
	raw_dict = dict(zip(raw_df['center'],raw_df['f']))

	#add a k-v pair, denoting the current chromosome 
	raw_dict['current_chrom'] = current_chrom
	return raw_dict

def calculate_mean_f(raw_dict):
	"""this function is to calculate the mean F rate across all CPG sites in one chromosome
	Input:
		raw_dict
	Ouput:
		the average F rate of this chromosomes
		the F sum of this chromsome 
		the CPG site counts of this chromosome
	"""
	sum_f_this_chrom = 0
	CPG_count_this_chrom = 0

	for site, f in raw_dict.items():
		if site != 'current_chrom':
			sum_f_this_chrom += f
			CPG_count_this_chrom += 1

	average_f_this_chrom = sum_f_this_chrom / CPG_count_this_chrom

	return sum_f_this_chrom, CPG_count_this_chrom, average_f_this_chrom

def calculate_CDS_mean_f(full_annotation_df,raw_dict):
	"""this function is to calculate the mean F rate in CDS area in one chromosome 
	Input:
		full_annotation_df
		raw_dict 
	output: 
		sum_f_this_chrom
		CPG_count_this_chrom
		meanF_in_CDS_this_chrom 
	"""
	sum_f_this_chrom = 0
	CPG_count_this_chrom = 0

	#get the current chromsome
	curr_chrom = raw_dict['current_chrom']

	#filter the annotation for CDS in this current chromosome 
	full_CDS_df = full_annotation_df[full_annotation_df['type'].str.contains('CDS')]
	this_chrom_CDS_df = full_CDS_df[full_CDS_df['chromosome'].str.contains(curr_chrom)].copy()

	#iterate through every CPG site, check if it is in CDS region
	for loc, f in raw_dict.items():
		if loc != 'current_chrom':
			is_this_in_CDS = ((this_chrom_CDS_df['start'] <= loc) & (this_chrom_CDS_df['end'] >= loc)).any()
			if is_this_in_CDS == True:
				sum_f_this_chrom += f
				CPG_count_this_chrom += 1

	if CPG_count_this_chrom != 0:
		meanF_in_CDS_this_chrom = sum_f_this_chrom/CPG_count_this_chrom
	else:
		meanF_in_CDS_this_chrom = 0
		
	return sum_f_this_chrom,CPG_count_this_chrom,meanF_in_CDS_this_chrom

def calculate_nc_mean_f(full_annotation_df,raw_dict):
	"""this function is to calculate the mean F rate outside CDS area 
	Input:
		full_annotation_df
		raw_dict 
	output: 
		sum_f_this_chrom
		CPG_count_this_chrom
		meanF_in_nc_this_chrom 

	issue!: this may be merged with the calculate_CDS_mean_f()?
	"""

	sum_f_this_chrom = 0
	CPG_count_this_chrom = 0

	#get the current chromsome
	curr_chrom = raw_dict['current_chrom']
	

	#filter the annotation for CDS in this current chromosome
	full_CDS_nc_df = full_annotation_df[full_annotation_df['type'].str.contains('CDS')]
	#full_gene_nc_df = full_annotation_df[full_annotation_df['type'].str.contains('gene')]

	this_CDS_nc_df = full_CDS_nc_df[full_CDS_nc_df['chromosome'].str.contains(curr_chrom)].copy()
	#this_gene_nc_df = full_gene_nc_df[full_gene_nc_df['chromosome'].str.contains(curr_chrom)].copy()

	#iterate through every CPG site, check if it is in CDS region
	for loc, f in raw_dict.items():
		if loc != 'current_chrom':
			nc_is_this_in_CDS = ((this_CDS_nc_df['start'] <= loc) & (this_CDS_nc_df['end'] >= loc)).any()
			#nc_is_this_in_gene = ((this_gene_nc_df['start'] <= loc) & (this_gene_nc_df['end'] >= loc)).any()

			#if this site do not fall into CDS then all consider as nc region 
			if nc_is_this_in_CDS == False:
				sum_f_this_chrom += f
				CPG_count_this_chrom += 1

	meanF_in_nc_this_chrom = sum_f_this_chrom/CPG_count_this_chrom


	return sum_f_this_chrom,CPG_count_this_chrom,meanF_in_nc_this_chrom	

def calculate_gene_and_regulatory_mean_f(full_annotation_df,raw_dict,extra_region = 5000):
	"""this function is to calculate mean f in genes and their regulatory region "extra" area 
	(default upstream 5000bp)
	Input:
		full_annotation_df
		raw_dict 
		extra_region: range of cis regulatory regions 
	output: 
		sum_f_this_chrom
		CPG_count_this_chrom
		meanF_in_CDS_this_chrom 
	"""
	sum_f_this_chrom = 0
	CPG_count_this_chrom = 0

	#get the current chromsome
	curr_chrom = raw_dict['current_chrom']
	
	#filter the annotation for gene in this current chromosome 
	full_extra_df = full_annotation_df[full_annotation_df['type'].str.contains('gene')]
	
	#adjust the start positions upstream 5000 bp 
	this_chrom_extra_df = full_extra_df[full_extra_df['chromosome'].str.contains(curr_chrom)].copy() #.copy to aviod SettingWithCopyWarning
	
	this_chrom_extra_df.loc[:,'start'] = this_chrom_extra_df.loc[:,'start'].subtract(extra_region)

	#iterate through every CPG site, check if it is in regulatory + gene regions
	for loc, f in raw_dict.items():
		if loc != 'current_chrom':
			is_this_in_genextra = ((this_chrom_extra_df['start'] <= loc) & (this_chrom_extra_df['end'] >= loc)).any()
			if is_this_in_genextra == True:
				sum_f_this_chrom += f
				CPG_count_this_chrom += 1

	meanF_in_genextra_this_chrom = sum_f_this_chrom/CPG_count_this_chrom

	
	return sum_f_this_chrom,CPG_count_this_chrom,meanF_in_genextra_this_chrom

def calculate_cis_mean_f(full_annotation_df,raw_dict,extra_region = 5000):
	"""this function is to calculate mean f in regulatory region
	(default upstream 5000bp)
	Input:
		full_annotation_df
		raw_dict 
		extra_region: range of cis regulatory regions 
	output: 
		sum_f_this_chrom
		CPG_count_this_chrom
		meanF_in_CDS_this_chrom 
	"""
	sum_f_this_chrom = 0
	CPG_count_this_chrom = 0

	#get the current chromsome
	curr_chrom = raw_dict['current_chrom']

	#filter the annotation for gene in this current chromosome 
	full_extra_df = full_annotation_df[full_annotation_df['type'].str.contains('gene')]
	this_chrom_cis_df = full_extra_df[full_extra_df['chromosome'].str.contains(curr_chrom)].copy() #.copy() to aviod SettingWithCopyWarning

	#change start position to upstream 5000 bp
	#change end position to original gene start position
	this_chrom_cis_df.loc[:,'end'] = this_chrom_cis_df.loc[:,'start']
	this_chrom_cis_df.loc[:,'start'] = this_chrom_cis_df.loc[:,'start'].subtract(extra_region)

	#iterate through every CPG site, check if it is in regulatory regions
	for loc, f in raw_dict.items():
		if loc != 'current_chrom':
			is_this_in_genextra = ((this_chrom_cis_df['start'] <= loc) & (this_chrom_cis_df['end'] >= loc)).any()
			if is_this_in_genextra == True:
				sum_f_this_chrom += f
				CPG_count_this_chrom += 1

	meanF_in_genextra_this_chrom = sum_f_this_chrom/CPG_count_this_chrom

	return sum_f_this_chrom,CPG_count_this_chrom,meanF_in_genextra_this_chrom

def calculate_gene_mean_f(full_annotation_df,raw_dict):
	"""this function is to calculate mean f in gene region
	Input:
		full_annotation_df
		raw_dict 
	output: 
		sum_f_this_chrom
		CPG_count_this_chrom
		meanF_in_CDS_this_chrom 
	"""
	sum_f_this_chrom = 0
	CPG_count_this_chrom = 0

	#get the current chromsome
	curr_chrom = raw_dict['current_chrom']

	#filter the annotation for gene in this current chromosome 
	full_gene_df = full_annotation_df[full_annotation_df['type'].str.contains('gene')]
	this_chrom_gene_df = full_gene_df[full_gene_df['chromosome'].str.contains(curr_chrom)].copy()

	#iterate through every CPG site, check if it is in gene regions
	for loc, f in raw_dict.items():
		if loc != 'current_chrom':
			is_this_in_gene = ((this_chrom_gene_df['start'] <= loc) & (this_chrom_gene_df['end'] >= loc)).any()
			if is_this_in_gene == True:
				sum_f_this_chrom += f
				CPG_count_this_chrom += 1

	meanF_in_gene_this_chrom = sum_f_this_chrom/CPG_count_this_chrom

	return sum_f_this_chrom,CPG_count_this_chrom,meanF_in_gene_this_chrom

#get annotation file
annotation_handle = sys.argv[1]
#get sample name
sample = sys.argv[2]

full_annotation_df = read_genome(annotation_handle)

#fetch chromsomes lists 
chrom_list = []
with open("/crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/chromsomes_list.txt", 'r') as f0:
	for line in f0:
		chrom_list.append(line.strip())
#specify NCPG parameters 
NCPGs = ["25","50","10","100"]
header = ["chromosome","CDS","non_coding_gene","gene","cis","cis+gene","all_cpg_sites"]

for NCPG in NCPGs:
	all_f_sum = 0
	all_f_counter = 0
	all_f_sum_auto = 0
	all_f_counter_auto = 0

	CDS_f_sum = 0
	CDS_f_counter = 0
	CDS_f_sum_auto = 0
	CDS_f_counter_auto = 0

	nc_f_sum = 0
	nc_f_counter = 0
	nc_f_sum_auto = 0
	nc_f_counter_auto = 0

	gene_f_sum = 0
	gene_f_counter = 0
	gene_f_sum_auto = 0
	gene_f_counter_auto = 0

	cis_f_sum = 0
	cis_f_counter = 0
	cis_f_sum_auto = 0
	cis_f_counter_auto = 0

	genextra_f_sum = 0
	genextra_f_counter = 0
	genextra_f_sum_auto = 0
	genextra_f_counter_auto = 0

	result_dict = {}

	for chrom in chrom_list:
		result_dict[chrom] = []
		full_f_name = 'As_'+sample +'_'+ chrom + '_' + NCPG +'.'+chrom+'.F'
		raw_F_dict = read_raw_f(full_f_name)

		#separate autosomes overall and all chromsomes counts


		#CDS region calculation
		this_CDS_f_sum, this_CDS_f_counter, this_CDS_mean = calculate_CDS_mean_f(full_annotation_df,raw_F_dict)
		CDS_f_sum += this_CDS_f_sum
		CDS_f_counter += this_CDS_f_counter
		result_dict[chrom].append(str(round(this_CDS_mean,4)))

		#w/in gene non-coding region mean F
		this_nc_f_sum, this_nc_f_counter, this_nc_mean = calculate_nc_mean_f(full_annotation_df,raw_F_dict)
		nc_f_sum += this_nc_f_sum
		nc_f_counter += this_nc_f_counter
		result_dict[chrom].append(str(round(this_nc_mean,4)))

		#gene region mean F
		this_gene_f_sum, this_gene_f_counter, this_gene_mean = calculate_gene_mean_f(full_annotation_df,raw_F_dict)
		gene_f_sum += this_gene_f_sum
		gene_f_counter += this_gene_f_counter
		result_dict[chrom].append(str(round(this_gene_mean,4)))

		#cis regulartory region mean F 
		this_cis_f_sum, this_cis_f_counter, this_cis_f_mean = calculate_cis_mean_f(full_annotation_df,raw_F_dict)
		cis_f_sum += this_cis_f_sum
		cis_f_counter += this_cis_f_counter
		result_dict[chrom].append(str(round(this_cis_f_mean,4)))

		#gene and it upstream 5kb cis regulatory region mean F
		this_genextra_f_sum, this_genextra_f_counter,this_genextra_mean = calculate_gene_and_regulatory_mean_f(full_annotation_df,raw_F_dict)
		genextra_f_sum += this_genextra_f_sum
		genextra_f_counter += this_genextra_f_counter
		result_dict[chrom].append(str(round(this_genextra_mean,4)))

		#All CPG site mean F
		this_all_f_sum, this_all_f_counter, this_all_f_mean = calculate_mean_f(raw_F_dict)
		result_dict[chrom].append(str(round(this_all_f_mean,4)))
		all_f_sum += this_all_f_sum
		all_f_counter += this_all_f_counter

		#exclude two sex chromsomes to calculate autosomes data 
		if chrom != 'NC_064846.1' and chrom != 'NC_064847.1':
			all_f_sum_auto += this_all_f_sum
			all_f_counter_auto += this_all_f_counter

			CDS_f_sum_auto += this_CDS_f_sum
			CDS_f_counter_auto += this_CDS_f_counter

			nc_f_sum_auto += this_nc_f_sum
			nc_f_counter_auto += this_nc_f_counter

			gene_f_sum_auto += this_gene_f_sum
			gene_f_counter_auto += this_gene_f_counter

			cis_f_sum_auto += this_cis_f_sum
			cis_f_counter_auto += this_cis_f_counter

			genextra_f_sum_auto += this_genextra_f_sum
			genextra_f_counter_auto += this_genextra_f_counter


	result_dict['Autosomes'] = []
	result_dict['Autosomes'].append(str(round(CDS_f_sum_auto/CDS_f_counter_auto,4)))
	result_dict['Autosomes'].append(str(round(nc_f_sum_auto/nc_f_counter_auto,4)))
	result_dict['Autosomes'].append(str(round(gene_f_sum_auto/gene_f_counter_auto,4)))
	result_dict['Autosomes'].append(str(round(cis_f_sum_auto/cis_f_counter_auto,4)))
	result_dict['Autosomes'].append(str(round(genextra_f_sum_auto/genextra_f_counter_auto,4)))
	result_dict['Autosomes'].append(str(round(all_f_sum_auto/all_f_counter_auto,4)))

	result_dict['Overall'] = []
	result_dict['Overall'].append(str(round(CDS_f_sum/CDS_f_counter,4)))
	result_dict['Overall'].append(str(round(nc_f_sum/nc_f_counter,4)))
	result_dict['Overall'].append(str(round(gene_f_sum/gene_f_counter,4)))
	result_dict['Overall'].append(str(round(cis_f_sum/cis_f_counter,4)))
	result_dict['Overall'].append(str(round(genextra_f_sum/genextra_f_counter,4)))
	result_dict['Overall'].append(str(round(all_f_sum/all_f_counter,4)))

	print(sample,NCPG)
	print(result_dict)
	#specify the result file names 
	result_f_handle ="raw_F_overview_"+ sample + '_' + NCPG 
	path = "/crex/proj/snic2022-6-144/nobackup/YUXIN/DamMet/raw_F_overview/"
	#format the dictionary into a single file 
	with open(path+result_f_handle,'w') as w1:
		w1.write('\t'.join(header) + '\n')
		for chrom, value_list in result_dict.items():
			w1.write(chrom + '\t' + '\t'.join(value_list) + '\n')


#local test
"""
full_annotation_df = read_genome('GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff')
sample = read_raw_f('As_RG.E.maximus.T18_NC_064844.1_25.NC_064844.1.F')

CDS_1, CDS_2, CDS_3 = calculate_CDS_mean_f(full_annotation_df,sample)
nc_1, nc_2, nc_3 = calculate_nc_mean_f(full_annotation_df,sample)
genextra_1, genextra_2, genextra_3 = calculate_gene_and_regulatory_mean_f(full_annotation_df,sample)

result_dict = {}
result_dict['this_sample'] = []
result_dict['this_sample'].append(str(round(CDS_3,4)))
result_dict['this_sample'].append(str(round(nc_3,4)))
result_dict['this_sample'].append(str(round(genextra_3,4)))
print(result_dict)

with open('rawf_against_ref_test','w') as w1:
	w1.write('44' + '\t' + 'mean_F_at_CDS' + '\t' + 'mean_F_at_non-coding-region' + '\t' + 'mean_F_at_gene_and_regulatory_region' + '\n')
	w1.write('this_sample' + '\t' + '\t'.join(result_dict['this_sample']) + '\n')
"""

#calculate_CDS_mean_f(full_annotation_df,raw_dict)
#calculate_nc_mean_f(full_annotation_df,raw_dict)
#calculate_gene_and_regulatory_mean_f(full_annotation_df,raw_dict)

#print(full_annotation_df['type'].unique())
"""
'region' 'gene' 'lnc_RNA' 'exon' 'pseudogene' 'mRNA' 'CDS' 'transcript'
 'snRNA' 'tRNA' 'snoRNA' 'cDNA_match' 'ncRNA' 'rRNA' 'V_gene_segment'
 'C_gene_segment'
"""
#calculate_CDS_mean_f(full_annotation_df,raw_dict)




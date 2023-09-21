#!/usr/bin/python3
from __future__ import division
import pandas as pd
import sys

"""
This script is used to get the distribution of each gene / promoter regions average F 
rate
input:
	genome annotation file (gff)
	DamMet raw F file results 
output: 
	tsv files contains target regions (start, end, chromosome, info) and their respectively average F (rawF / number of CPG sites) 
"""

def read_annotation(gff_handle):
	"""
	this function is used to read gff file Asian elephant genome
	input: 
		gff_handle: file path to reference annotation file 
	output:
		full_annotation_df
	"""
	raw_annotation_df = pd.read_table(gff_handle, comment = "#", skiprows = [0], header = None)
	#print(raw_annotation_df.iloc[[1],:])
	full_annotation_df = raw_annotation_df.iloc[:,[0,2,3,4,8]]
	new_col_names = {0:"chromosome",2:"type",3:"start",4:"end",8:"Target_Info"}
	full_annotation_df = full_annotation_df.rename(columns = new_col_names)
	return full_annotation_df

def get_full_gene(full_annotation_df):
	"""this function is to get the full annotated gene dataframe 
	input: full_annotation_df
	output: full_gene_df
	"""
	full_gene_df = full_annotation_df[full_annotation_df['type'] == 'gene'].copy()
	
	return full_gene_df

def get_full_CDS(full_annotation_df):
	"""this function is to get the full CDS region dataframe
	input: full_annotation_df
	output: full_CDS_df
	"""
	full_CDS_df = full_annotation_df[full_annotation_df['type'] == 'CDS'].copy()
	return full_CDS_df

def get_full_promoter(full_gene_df):
	"""This function is to get the full promoter dataframe modified from full gene dataframe
	default regulatory regions: upstream 5000 bp
	input: full_gene_df
	output: full_promoter_df
	"""
	full_promoter_df = full_gene_df.copy()
	full_promoter_df['end'] = full_gene_df['start']
	full_promoter_df['start'] = full_gene_df['start'].subtract(5000)
	
	return full_promoter_df

def read_raw_f(raw_f_handle):
	"""This function is to read DamMet raw F files only on center and f rate columns 
	"""
	raw_df = pd.read_table(raw_f_handle, sep = ' ', dtype = 'str', usecols = [0,5])
	raw_df['f'] = raw_df['f'].astype(float)
	raw_df['center'] = raw_df['center'].str.split(':').str.get(-1).astype(int)

	return raw_df

def get_mean_F(full_df,sample_name,NCPG,mode):
	"""This function is to caluclate all regions in full_df (CDS, promoters or genes)
	input: 
		full_df (CDS, promoter, genes)
		sample_name 
		NCPG: DamMet parameters used for raw F file to read
		mode: "CDS"or "promoter" or "gene"
	output: 
		full_df with an additional column 'Average_F_rate'	 

	"""
	n = len(full_df.index)
	meanF_list = []

	#iterate through region (row)
	for i in range(n):
		region_start = full_df['start'].iloc[i]
		region_end = full_df['end'].iloc[i]

		#avoid NW chromosomes 
		if full_df['chromosome'].iloc[i][:3] == 'NC_':
			region_chrom = full_df['chromosome'].iloc[i]

			#access the corresponding raw F file
			rawF_handle = "As_RG." + sample_name + "_" + region_chrom + '_' + NCPG + '.' + region_chrom + '.F'
			rawF_df = read_raw_f(rawF_handle)

			#select the CPG sites are within the region
			region_condition = (rawF_df['center'] > region_start) & (rawF_df['center'] < region_end)
			temp_region_df = rawF_df[region_condition]

			#check if there is any CPG sites fall within the region 
			if len(temp_region_df.index) != 0:
				region_avg_f = temp_region_df['f'].sum()/len(temp_region_df.index)
			else:
				region_avg_f = 'NA'

			meanF_list.append(region_avg_f)
		else:
			meanF_list.append('NA')

	full_df.insert(5,'Average_F_rate',meanF_list)
	result_name = "{mode}_{sample}_averageF".format(mode = mode, sample = sample_name)
	full_df.to_csv(result_name,sep = '\t')

def get_onechrom_gene(full_annotation_df,chrom_name):
	"""
	This function is to get specific chromosome annotated gene
	input:
		full_annotation_df 
		chrom_name
	output:
		partial_gene_df 
	"""
	full_gene_df = full_annotation_df[full_annotation_df['type'] == 'gene'].copy()
	partial_gene_df = full_gene_df[full_gene_df['chromosome'] == chrom_name].copy()
	return partial_gene_df

def get_onechrom_promoter(partial_gene_df):
	"""This function is to get the full promoter dataframe modified from full gene dataframe
	input: partial_gene_df
	output: partial_promoter_df
	"""
	partial_promoter_df = partial_gene_df.copy()
	partial_promoter_df['end'] = partial_gene_df['start']
	partial_promoter_df['start'] = partial_gene_df['start'].subtract(5000)

	return partial_promoter_df

#specify annotation file
gff_handle = sys.argv[1]

#specify sample name 
sample_name = sys.argv[2]
#gff_handle = "GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff"

# for doing only a subset of annotaions 
single_chrom = "NC_064830.1"
full_annotation_df = read_annotation(gff_handle)
single_gene_df = get_onechrom_gene(full_annotation_df,single_chrom)
single_promoter_df = get_onechrom_promoter(single_gene_df)

NCPGs = ['25']

for NCPG in NCPGs:
	get_mean_F(single_gene_df,sample_name,NCPG,'Gene')
	get_mean_F(single_promoter_df,sample_name,NCPG,'Promoter')

#to run whole genome annotations 
"""
full_annotation_df = read_annotation(gff_handle)
full_gene_df = get_full_gene(full_annotation_df)
full_CDS_df = get_full_CDS(full_annotation_df)
full_promoter_df = get_full_promoter(full_gene_df)

NCPGs = ['25']

for NCPG in NCPGs:
	get_mean_F(full_gene_df,sample_name,NCPG,'Gene')
	get_mean_F(full_CDS_df,sample_name,NCPG,'CDS')
	get_mean_F(full_promoter_df,sample_name,NCPG,'Promoter')

"""


#local test
"""
gff_handle = "GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff"
def search_by_id(filtered_df,gene_id):
	result = filtered_df[filtered_df['Target_Info'].str.contains(gene_id)]
	print(result)

def search_info_occurance(filtered_df,info):
	#This function is to parse the candidate df by the description in 'Target into'
	
	m = filtered_df['Target_Info'].str.contains(info).sum()
	n = len(filtered_df.index)
	print("{Info} has total {m} occurance out of {n}".format(Info = info, m = m, n = n))
	print(m/n)
	return float(m),float(n),float(m/n)
full_annotation_df = read_annotation(gff_handle)
#search_info_occurance(full_annotation_df,'olfactory')

full_gene_df = get_full_gene(full_annotation_df)
#print(full_gene_df['chromosome'].iloc[5][:3])
#print(len(full_gene_df.index))
test_promoter_df = get_full_promoter(full_gene_df).tail(10)
get_mean_F(test_promoter_df,"E.maximus.T13",'25','zPromoter')
"""





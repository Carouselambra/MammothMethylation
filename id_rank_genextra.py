#!/usr/bin/python3
from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import itertools
import sys
from scipy import stats
"""
This script is used to calculate regional averaged f rate (gene or promoter region)
of each samples

gene / cis regulatory  -> windows w/in 

input: 
	annotation gff file
	a pair of normalised F files (sample1 vs sample 2)
output: 
	a tsv file list of gene including the differences between F fold changes 
	a tsv file list of cis region ranked by  

"""

def read_annotation(annotation_handle):
	"""This function is to read the reference genome annotaion 
	this function is used to read gff file Asian elephant genome
	"""
	raw_annotation_df = pd.read_table(annotation_handle, comment = "#", skiprows = [0], header = None)
	full_annotation_df = raw_annotation_df.iloc[:,[0,2,3,4,8]].copy()
	new_col_names = {0:"chromosome",2:"type",3:"start",4:"end",8:"ID_info"}
	full_annotation_df = full_annotation_df.rename(columns = new_col_names)
	return full_annotation_df

def read_normalisedF(nF_handle):
	"""This function is to read normalised F files with columns ['window_location','average_f_rate']
	Input:
		nF_handle (str): the normalised F file handle
	Output:
		normalise F dataframe
	"""
	nF_df = pd.read_csv(nF_handle,sep = '\t', usecols = ['window_location','normalised_f_fold'])
	return nF_df

def merged_normalisedF(nF_df_list,out_handle):
	"""This function is to merge multiple normalised F files from the same spiece 
	into a new merged+averaged F files 
	input: 
		nF_df_list: a list of dataframes belong to multipe samples of the same species (should be on the same chromsomes)
		out_handle: the output file path
	output:
		out_df
		and a tsv file with one merged normalised F from multiple samples 
	"""
	sample1_df = nF_df_list[0]
	sample2_df = nF_df_list[1]
	merging_df = sample1_df.merge(sample2_df,on = 'window_location', how = 'inner')

	for i in range(2,len(nF_df_list)):
		nF_df_list[i].rename({'normalised_f_fold':'normalised_f_fold{num}'.format(num = str(i+1))},axis=1, inplace=True)
		merging_df = merging_df.merge(nF_df_list[i], on = 'window_location', how = 'inner',suffixes = (' ', ' '))
		merging_df.columns = merging_df.columns.str.strip()

	merging_df['normalised_f_fold'] = merging_df.iloc[:,1:].mean(axis = 1)
	out_df = merging_df[['window_location','normalised_f_fold']]
	out_df.to_csv(out_handle, sep = '\t')
	return out_df

def read_combine_normaliseFs(nF_handle_1,nF_handle_2):
	"""
	this function is used to read two merged normalised F files into \
	one dataframe with columns ['window_locations','normalised_f_fold_1','normalised_f_fold_2'] 
	input: 
		nF_handle_1
		nF_handle_2
	output:
		pairwise_normalise_df: a dataframe with columns ['window_locations','normalised_f_fold_1','normalised_f_fold_2'] 
	"""
	sample1_df = pd.read_csv(nF_handle_1,sep = '\t', usecols = ['window_location','normalised_f_fold'])
	sample2_df = pd.read_csv(nF_handle_2,sep = '\t', usecols = ['window_location','normalised_f_fold'])

	pairwise_normalise_df = sample1_df.merge(sample2_df,on = 'window_location',how = 'inner',suffixes = ('_1','_2'))
	return pairwise_normalise_df

def get_gene_annotation(full_annotation_df):
	full_gene_df = full_annotation_df[full_annotation_df['type'].str.contains('gene')].copy()
	return full_gene_df

def get_cis_annotation(full_annotation_df):
	full_cis_df = get_gene_annotation(full_annotation_df)
	full_cis_df['end'] = full_cis_df['start']
	full_cis_df['start'] = full_cis_df['start'].subtract(5000)
	return full_cis_df

def rank_targets(target_df,pairwise_normalise_df,test_mode = "off"):
	"""This function is to get a list of genes / cis region (key) ranked by the differences in normalised F rate
	input:
		target_df: specify regions to calculate the area f rate, extracted from reference annotation file 
		pairwise_normalise_df: a merged normalised f dataframe, columns = ['window_location','normalised_f_fold_1','normalised_f_fold_2']

	output:
		result_df: a result dataframe containing columns ['Target_info','Diff_in_f_mean']
	"""
	if test_mode != "off":
		n = int(test_mode)
	else:
		n = len(target_df.index)
	result_dict = {}
	for i in range(n):
		#choose the target - get start&end
		start = target_df['start'].iloc[i]
		end = target_df['end'].iloc[i]
		chrom = target_df['chromosome'].iloc[i]

		#set filter for windows falls into ranges of items in the target_df (either gene or cis regions)
		target_mask = pairwise_normalise_df['window_location'].apply(lambda x: int(x.split('-')[0]))\
					.between(start, end) & \
		   pairwise_normalise_df['window_location'].apply(lambda x: int(x.split('-')[1]))\
					.between(start, end)
		find_window_df = pairwise_normalise_df.loc[target_mask].copy()

		#get the genes/cis specific normalised F for two species 
		value1 = find_window_df['normalised_f_fold_1'].mean()
		value2 = find_window_df['normalised_f_fold_2'].mean()

		#get the normalised F fold differences between two species
		diff_num = value1-value2

		#get the absoulte value of the normalised F fold differences
		diff_num_abs = abs(diff_num)

		#get the additional info about the region 
		target_info = target_df['ID_info'].iloc[i]
		result_dict[target_info] = []

		#write the infomation into the dictionary 
		result_dict[target_info].append(diff_num_abs)
		result_dict[target_info].append(diff_num)
		result_dict[target_info].append(chrom)
		result_dict[target_info].append(start)
		result_dict[target_info].append(end)
	
	#convert the result dictionry into a dataframe 
	result_df = pd.DataFrame.from_dict(result_dict,orient = "index",columns = ['Diff_in_f_mean_abs','Diff_in_f_mean','chrom','start','end'])
	result_df['Target_Info'] = result_df.index
	result_df.reset_index(drop=True, inplace = True)

	#rank the regions by their absolute value of the normalised F fold
	result_df.sort_values(by = "Diff_in_f_mean_abs", ascending = False, inplace = True, ignore_index = True)
	return result_df

def sample_to_target(sample_list_1,sample_list_2):
	"""This function is to take two lists of sample names, then return the sorted list of genes / cis areas and their area binned average f rate
	input:
		sample_list_1:both should also have genus name!
		sample_list_2:
	output:
		two tsv files specify target regions (promoter cis and gene) locations, difference in f means
		(columns = ["Diff_in_f_mean_abs","Diff_in_f_mean","start","end",'chrom',"Target_Info"])
	"""
	sample1_label = sample_list_1[0]
	if len(sample_list_1) > 1:
		for i in sample_list_1[1:]:
			sample1_label = sample1_label +'.'+ i.rpartition('.')[-1]

	sample2_label = sample_list_2[0]	
	if len(sample_list_2) > 1:
		for i in sample_list_2[1:]:
			sample2_label = sample2_label +'.' +i.rpartition('.')[-1]

	#local test
	#NCPGs = ['25']
	#chrom_list = ['NC_064837.1','NC_064838.1']
	#result_filepath = ""
	#ANNOTATION = 'GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff'

	NCPGs = ["25","50","10","100"]
	#fetch the chromsome list
	chrom_list = []
	with open("/crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/chromsomes_list.txt", 'r') as f0:
		for line in f0:
			chrom_list.append(line.strip())

	#read annotation 
	ANNOTATION="/crex/proj/snic2022-6-144/nobackup/TOM/MAMMOTH/REF/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff"
	result_filepath = "/crex/proj/snic2022-6-144/nobackup/YUXIN/DamMet/id_result/"
	winsize = "1000"
	full_annotation_df = read_annotation(ANNOTATION)
	full_gene_df = get_gene_annotation(full_annotation_df)
	full_cis_df = get_cis_annotation(full_annotation_df)
	print("annotation get!")

	#iterate through all NCPG parameters
	for NCPG in NCPGs:
		#set up result dataframe
		final_result_gene_df = pd.DataFrame(columns = ["Diff_in_f_mean_abs","Diff_in_f_mean","start","end",'chrom',"Target_Info"])
		final_result_cis_df = pd.DataFrame(columns = ["Diff_in_f_mean_abs","Diff_in_f_mean","start","end",'chrom',"Target_Info"])

		#iterate through chromsomes
		for chrom in chrom_list:
			#get chromosome specific annotation
			this_chrom_gene_df = full_gene_df[full_gene_df['chromosome'] == chrom].reset_index()
			this_chrom_cis_df = full_cis_df[full_cis_df['chromosome'] == chrom].reset_index()

			#merging all the samples from the specie 1
			single_nF_list_sample1 = []
			if len(sample_list_1) > 1:
				for sample1 in sample_list_1:
					nF_handle = "normalised_binned_winsize_{window}_As_RG.".format(window = winsize) + sample1 + '_' + chrom + '_' + NCPG +'.'+chrom+'.F'
					single_nF_list_sample1.append(read_normalisedF(nF_handle))
				sample1s_df = merged_normalisedF(single_nF_list_sample1,"normlised_averaged_RG.{samples}_{chrom}_{NCPG}".format(samples = '.'.join(sample_list_1), chrom = chrom, NCPG = NCPG))
			else:
				nF_handle = "normalised_binned_winsize_{window}_As_RG.".format(window = winsize) + sample_list_1[0] + '_' + chrom + '_' + NCPG +'.'+chrom+'.F'
				sample1s_df = read_normalisedF(nF_handle)
			print(sample1_label+ '-' +chrom + " done!")

			#merging all the samples from the specie 2 
			single_nF_list_sample2 = []
			if len(sample_list_2) > 1:
				for sample2 in sample_list_2:
					nF_handle = "normalised_binned_winsize_{window}_As_RG.".format(window = winsize) + sample2 + '_' + chrom + '_' + NCPG +'.'+chrom+'.F'
					single_nF_list_sample2.append(read_normalisedF(nF_handle))
				sample2s_df = merged_normalisedF(single_nF_list_sample2,"normlised_averaged_RG.{samples}_{chrom}_{NCPG}".format(samples = '.'.join(sample_list_2), chrom = chrom, NCPG = NCPG))
			else:
				nF_handle = "normalised_binned_winsize_{window}_As_RG.".format(window = winsize) + sample_list_2[0] + '_' + chrom + '_' + NCPG +'.'+chrom+'.F'
				sample2s_df = read_normalisedF(nF_handle)
			print(sample2_label+chrom + "done!")

			#store the two merged datasets into one dataframe
			pairwise_normalise_df = sample1s_df.merge(sample2s_df,on = 'window_location',how = 'inner',suffixes = ('_1','_2'))
			print(pairwise_normalise_df.head(5))

			#calculate differences between the averaged window F scores 
			#rank all genes and promoters based on absoulte value of differences in normalised F fold changes 
			per_chrom_gene_result_df = rank_targets(this_chrom_gene_df,pairwise_normalise_df)
			per_chrom_cis_result_df = rank_targets(this_chrom_cis_df,pairwise_normalise_df)

			#combine results from all chromsomes 
			final_result_gene_df = pd.concat([final_result_gene_df,per_chrom_gene_result_df])
			final_result_cis_df = pd.concat([final_result_cis_df,per_chrom_cis_result_df])			

		#sort the all regions again (interchromsomes)
		final_result_gene_df.sort_values(by = "Diff_in_f_mean_abs",ascending = False, ignore_index = True)
		print("this is gene df preview")
		print(final_result_gene_df.head(5))

		#set outputfile name for genes
		outhandle_gene = result_filepath+"Candidates_Genes_{sample1s}_vs_{sample2s}_NCPG{NCPG}".format(sample1s = sample1_label,sample2s = sample2_label, NCPG = NCPG)
		final_result_gene_df.to_csv(outhandle_gene,sep = '\t')

		#set outputfile name for cis regulatory regions 
		final_result_cis_df.sort_values(by = "Diff_in_f_mean_abs",ascending = False, ignore_index = True)
		print("this is cis df preview")
		print(final_result_cis_df.head(5))
		outhandle_cis = result_filepath+"Candidates_Cis_{sample1s}_vs_{sample2s}_NCPG{NCPG}".format(sample1s = sample1_label,sample2s = sample2_label, NCPG = NCPG)
		final_result_cis_df.to_csv(outhandle_cis,sep = '\t')

#all sample names 
mm_samples = ['M.primigenius.E467','M.primigenius.L163','M.primigenius.L164','M.primigenius.P005','M.primigenius.M6']
ase_samples = ['E.maximus.T13','E.maximus.T16','E.maximus.T18']
afe_sample = ['L.africana.T11']
stsk_sample = ['P.antiquus.ERR2260504']
mstd_sample = ['M.americanum.ERR2260503']

#specify between which two speices the comparison needed to be made 
sample_specify = sys.argv[1]
if sample_specify == "MMs_vs_ASEs":
	#Mammoths vs. Asian Elephants
	sample_to_target(mm_samples,ase_samples)
elif sample_specify == "MMs_vs_AFE":
	#Mammoths vs. African elephant T11
	sample_to_target(mm_samples,afe_sample)
elif sample_specify == "MMs_vs_STKS":
	#Mammoths vs. straight tussk elephant sample
	sample_to_target(mm_samples,stsk_sample)
elif sample_specify == "MMs_vs_MSTD":
	#Mammoths vs. mastodon sample
	sample_to_target(mm_samples,mstd_sample)

"""
ase_samples = ["T13","T16","T18"]
mm_samples = ["E467","L163","L164","P005","M6"]

#NCPGs = ["25"]
NCPGs = ["10","25","50","100"]
#chrom_list = ['NC_064837.1','NC_064838.1']
chrom_list = []
with open("/crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/chromsomes_list.txt", 'r') as f0:
	for line in f0:
		chrom_list.append(line.strip())
winsize = "1000"
result_filepath = "/crex/proj/snic2022-6-144/nobackup/YUXIN/DamMet/id_result/"

ANNOTATION="/crex/proj/snic2022-6-144/nobackup/TOM/MAMMOTH/REF/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff"
full_annotation_df = read_annotation(ANNOTATION)
full_gene_df = get_gene_annotation(full_annotation_df)
full_cis_df = get_cis_annotation(full_annotation_df)

for NCPG in NCPGs:
	final_result_gene_df = pd.DataFrame(columns = ["Diff_in_f_mean_abs","Diff_in_f_mean","Target_Info"])
	final_result_cis_df = pd.DataFrame(columns = ["Diff_in_f_mean_abs","Diff_in_f_mean","Target_Info"])

	for chrom in chrom_list:
		single_nF_list_mm = []
		for mm in mm_samples:
			nF_handle = "normalised_binned_winsize_{window}_As_RG.M.primigenius.".format(window = winsize) + mm + '_' + chrom + '_' + NCPG +'.'+chrom+'.F'
			single_nF_list_mm.append(read_normalisedF(nF_handle))
		MMs_df = merged_normalisedF(single_nF_list_mm,"normlised_averaged_RG.M.primigenius.{samples}_{chrom}_{NCPG}".format(samples = '.'.join(mm_samples), chrom = chrom, NCPG = NCPG))
		#print(MMs_df.tail(5))

		single_nF_list_ase = []
		for ase in ase_samples:
			nF_handle = "normalised_binned_winsize_{window}_As_RG.E.maximus.".format(window = winsize) + ase + '_' + chrom + '_' + NCPG +'.'+chrom+'.F'
			single_nF_list_ase.append(read_normalisedF(nF_handle))
		ASEs_df = merged_normalisedF(single_nF_list_ase,"normlised_averaged_RG.E.maximus.{samples}_{chrom}_{NCPG}".format(samples = '.'.join(ase_samples), chrom = chrom, NCPG = NCPG))
		#print(ASEs_df.tail(5))

		pairwise_normalise_df = MMs_df.merge(ASEs_df,on = 'window_location',how = 'inner',suffixes = ('_1','_2'))
		#print(pairwise_normalise_df.tail(5))
		per_chrom_gene_result_df = rank_targets(full_gene_df,pairwise_normalise_df)
		per_chrom_cis_result_df = rank_targets(full_cis_df,pairwise_normalise_df)

		final_result_gene_df = pd.concat([final_result_gene_df,per_chrom_gene_result_df])
		final_result_cis_df = pd.concat([final_result_cis_df,per_chrom_cis_result_df])


	final_result_gene_df.sort_values(by = "Diff_in_f_mean_abs",ascending = False, ignore_index = True)
	print(final_result_gene_df.head(5))
	outhandle_gene = result_filepath+"Candidates_Genes_Mammoths_vs_AsianElephants_NCPG{NCPG}".format(NCPG = NCPG)
	final_result_gene_df.to_csv(outhandle_gene,sep = '\t')

	final_result_cis_df.sort_values(by = "Diff_in_f_mean_abs",ascending = False, ignore_index = True)
	print(final_result_cis_df.head(5))
	outhandle_cis = result_filepath+"Candidates_Cis_Mammoths_vs_AsianElephants_NCPG{NCPG}".format(NCPG = NCPG)
	final_result_cis_df.to_csv(outhandle_cis,sep = '\t')
"""



		
#local test
"""
full_annotation_df = read_annotation('GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff')
nF1_handle = 'normalised_binned_winsize_1000_As_RG.E.maximus.T13_NC_064837.1_25.NC_064837.1.F'
nF2_handle = 'normalised_binned_winsize_1000_As_RG.E.maximus.T16_NC_064837.1_25.NC_064837.1.F'
nF3_handle = 'normalised_binned_winsize_1000_As_RG.E.maximus.T18_NC_064837.1_25.NC_064837.1.F'

nFmm1_handle = 'normalised_binned_winsize_1000_As_RG.M.primigenius.M6_NC_064837.1_25.NC_064837.1.F'
nFmm2_handle = 'normalised_binned_winsize_1000_As_RG.M.primigenius.E467_NC_064837.1_25.NC_064837.1.F'
nFmm3_handle = 'normalised_binned_winsize_1000_As_RG.M.primigenius.P005_NC_064837.1_25.NC_064837.1.F'
nFmm4_handle = 'normalised_binned_winsize_1000_As_RG.M.primigenius.L163_NC_064837.1_25.NC_064837.1.F'
nFmm5_handle = 'normalised_binned_winsize_1000_As_RG.M.primigenius.L164_NC_064837.1_25.NC_064837.1.F'

combinedf = read_combine_normaliseFs(nF1_handle,nFmm1_handle)
combinedf_2 = read_combine_normaliseFs(nF2_handle,nFmm2_handle)
target_df = get_gene_annotation(full_annotation_df)


result_df_1 = rank_targets(target_df,combinedf,test_mode = "100")
result_df_2 = rank_targets(target_df,combinedf_2,test_mode = "100")
final = pd.concat([result_df_1,result_df_2]).sort_values(by = "Diff_in_f_mean",ascending = False)
print(final.head(5))
"""





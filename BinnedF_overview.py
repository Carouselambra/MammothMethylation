#!/usr/bin/python3
from __future__ import division
import pandas as pd
import numpy as np
import sys
#import matplotlib.pyplot as plt
#from matplotlib.pyplot import figure
#import itertools

"""
This script is used to calculate 
	the average F rates of an entire binned F files
	the average F rate of an whole genome (all binned F files)
Input:
	binned F files of a sample (same NCPGs and bin window size)
Output:
	a txt file: chromosome name - average binned F rate 
"""
"""
for each genome, get 
"""
def get_single_binf_mean(single_bin_df):
	single_bin_counter = len(single_bin_df.index)
	single_bin_mean = single_bin_df['average_f_rate'].mean()
	single_bin_sum = single_bin_df['average_f_rate'].sum()
	print(single_bin_mean,single_bin_counter)

	return single_bin_mean, single_bin_counter, single_bin_sum

#get "mean" bin f (Asian elephants, mammoths)

sample_1 = 'RG.E.maximus.T13.T16.T18'
#sample_2 = 'RG.M.primigenius.E467.L163.L164.P005.M6'
chrom_list = []
with open("/crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/chromsomes_list.txt", 'r') as f0:
	for line in f0:
		chrom_list.append(line.strip())
		
NCPGs = ["10","25","50","100"]
for NCPG in NCPGs:
	result_dict_1 = {}
	#result_dict_2 = {}

	all_mean_sum_1 = 0
	#all_mean_sum_2 = 0
	all_mean_sum_1_auto = 0
	#all_mean_sum_2_auto = 0


	all_counter_1 = 0
	#all_counter_2 = 0
	all_counter_1_auto = 0
	#all_counter_2_auto = 0

	for chrom in chrom_list:
		full_name_1 = "Mean_binned_winsize_1000_As_" + sample_1 + '_' +chrom+ "_NCPG" +NCPG
		#full_name_2 = "Mean_binned_winsize_1000_As_" + sample_2 + '_' +chrom+ "_NCPG" +NCPG

		df_1 = pd.read_table(full_name_1,usecols = ['merged_average_f_rate'])
		df_1.rename(columns = {'merged_average_f_rate':'average_f_rate'},inplace = True)

		#df_2 = pd.read_table(full_name_2,usecols = ['merged_average_f_rate'])
		#df_2.rename(columns = {'merged_average_f_rate':'average_f_rate'},inplace = True)

		this_mean_1, this_counter_1, this_sum_1 = get_single_binf_mean(df_1)

		#this_mean_2, this_counter_2, this_sum_2 = get_single_binf_mean(df_2)

		result_dict_1[chrom] = this_mean_1
		#result_dict_2[chrom] = this_mean_2

		all_mean_sum_1 += this_sum_1
		#all_mean_sum_2 += this_sum_2

		all_counter_1 += this_counter_1
		#all_counter_2 += this_counter_2

		if chrom != 'NC_064846.1' and chrom != 'NC_064847.1':
			all_mean_sum_1_auto += this_sum_1
			#all_mean_sum_2_auto += this_sum_2

			all_counter_1_auto += this_counter_1
			#all_counter_2_auto += this_counter_2


	result_dict_1['autosomes_average_f_rate'] = all_mean_sum_1_auto/all_counter_1_auto
	#result_dict_2['autosomes_average_f_rate'] = all_mean_sum_2_auto/all_counter_2_auto

	result_dict_1['whole_genome_average_f_rate'] = all_mean_sum_1/all_counter_1
	#result_dict_2['whole_genome_average_f_rate'] = all_mean_sum_2/all_counter_2

	print(result_dict_1['whole_genome_average_f_rate'])
	#print(result_dict_2['whole_genome_average_f_rate'])

	out_f_1 = 'Overview_mean_binnedFs_winsize_1000_averaged_' + sample_1 + "_NCPG" +NCPG
	#out_f_2 = 'Overview_mean_binnedFs_winsize_1000_averaged_' + sample_2 + "_NCPG" +NCPG
	path = '/crex/proj/snic2022-6-144/nobackup/YUXIN/DamMet/binned_F'
	with open(out_f_1,'w') as w1:
		for k,v in result_dict_1.items():
			w1.write(str(k) + '\t' + str(v) +'\n')

	#with open(out_f_2,'w') as w2:
		#for k,v in result_dict_2.items():
			#w2.write(str(k) + '\t' + str(v) +'\n')


#get overview for a single bin F 
"""
sample_name = sys.argv[1]
NCPGs = ["10","25","50","100"]
chrom_list = []
with open("/crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/chromsomes_list.txt", 'r') as f0:
	for line in f0:
		chrom_list.append(line.strip())
winsize = "1000"
for NCPG in NCPGs:
	result_dict = {}
	all_mean_sum = 0
	all_mean_sum_auto = 0
	all_counter = 0
	all_counter_auto = 0

	for chrom in chrom_list:
		binned_f_name = "binned_winsize_{window}_As_".format(window = winsize) + sample_name + '_' + chrom + '_' + NCPG +'.'+chrom+'.F'
		df = pd.read_table(binned_f_name, usecols = ['average_f_rate'])
		this_mean, this_counter, this_sum = get_single_binf_mean(df)
		result_dict[chrom] = this_mean
		
		all_mean_sum += this_sum
		all_counter += this_counter

		if chrom != 'NC_064846.1' and chrom != 'NC_064847.1':
			all_mean_sum_auto += this_sum
			all_counter_auto += this_counter

	result_dict['autosomes_average_f_rate'] = all_mean_sum_auto/all_counter_auto
	result_dict['whole_genome_average_f_rate'] = all_mean_sum/all_counter

	out_f = 'Overview_mean_binnedFs_winsize_{windownum}_'.format(windownum = winsize)+sample_name+'_' + NCPG +'.txt'
	path = '/crex/proj/snic2022-6-144/nobackup/YUXIN/DamMet/binned_F'
	with open(out_f,'w') as w:
		for k,v in result_dict.items():
			w.write(str(k) + '\t' + str(v) +'\n')
"""
#df = pd.read_table("binned_winsize_1000_As_RG.merged.E.maximus.T13_T13_18_NC_064820.1_50.NC_064820.1.F")
#get_single_binf_mean(df)
#chromosome_list = 
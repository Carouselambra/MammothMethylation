#!/usr/bin/python3
from __future__ import division
import pandas as pd
import numpy as np
import sys

"""This script is used to merge several binned F file from the same species
	Asian elephants (T13, T16, T18)
	or Mammoths (E467, L163, L164, M6, P005)
	into one single binned F file
"""
#which_one = sys.argv[1]
mm_samples = ['E467','L163','L164','P005','M6']
el_samples = ['T13','T16','T18']

chrom_list = []
with open("/crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/chromsomes_list.txt", 'r') as f0:
	for line in f0:
		chrom_list.append(line.strip())
NCPGs = ["10","25","50","100"]


def merge_dfs(input_f_list,out_f_handle):
	"""This function is used to merged several binned F files (same NCPG, same chromosome) from the same speices
	input:
		input_f_list: a list of f files handles
		out_f_handle: the name of out_put_f_files
	output: 
		a tsv file with window locatino and merged_average_f_rate
	"""
	df_list = []
	for in_f_handle in input_f_list:
		df = pd.read_csv(in_f_handle, sep = '\t')
		df_list.append(df)
	
	merging_df = df_list[0].merge(df_list[1], on = 'window_location', how = 'inner',suffixes = ('_1','_2'))
	#print(merging_df)

	for i in range(2,len(df_list)):
		df_list[i].rename({'average_f_rate':'average_f_rate_{num}'.format(num = str(i+1))},axis=1, inplace=True)
		merging_df = merging_df.merge(df_list[i], on = 'window_location', how = 'inner',suffixes = (' ', ' '))
		merging_df.columns = merging_df.columns.str.strip()
	#print(merging_df.tail(20))
	#merging_df.dropna(0)
	merging_df['merged_average_f_rate'] = merging_df.iloc[:,1:].mean(axis = 1)

	print(merging_df['merged_average_f_rate'].mean())
	out_df = merging_df[['window_location','merged_average_f_rate']]
	#print(out_df['merged_average_f_rate'].mean())
	#print(merging_df.tail(20))
	out_df.to_csv(out_f_handle,index = False, sep = '\t')
	#merging_df.to_csv(out_f_handle,index = False, sep = '\t')

"""
for chrom in chrom_list:
	for NCPG in NCPGs:
		input_f_list = []
		#if which_one == 'mammoth':
		for mm in mm_samples:
			input_f = 'binned_winsize_1000_As_RG.M.primigenius.' + mm + '_' + chrom + '_' + NCPG + '.' + chrom + '.F'
			input_f_list.append(input_f)
		out_f_handle = 'Mean_binned_winsize_1000_As_RG.M.primigenius.{mms}_{chrom}_NCPG{NCPG}'.format(mms = '.'.join(map(str,mm_samples)), chrom = chrom, NCPG = NCPG)     
		merge_dfs(input_f_list,out_f_handle)
		#elif which_one == 'elephant':
		for el in el_samples:
			input_f = 'binned_winsize_1000_As_RG.E.maximus.' + el + '_' + chrom + '_' + NCPG + '.' + chrom + '.F'
			input_f_list.append(input_f)
		out_f_handle = 'Mean_binned_winsize_1000_As_RG.E.maximus.{els}_{chrom}_NCPG{NCPG}'.format(els = '.'.join(map(str,el_samples)), chrom = chrom, NCPG = NCPG) 
		merge_dfs(input_f_list,out_f_handle)
		#else:
		#	print('Input is not right!')
"""


filespath = "/crex/proj/snic2022-6-144/nobackup/YUXIN/DamMet/binned_F/"
for chrom in chrom_list:
	for NCPG in NCPGs:
		input_f_list = []
		for mm in mm_samples:
			input_f = filespath+'binned_winsize_1000_As_RG.M.primigenius.' + mm + '_' + chrom + '_' + NCPG + '.' + chrom + '.F'
			input_f_list.append(input_f)
		out_f_handle = 'Mean_binned_winsize_1000_As_RG.M.primigenius.{mms}_{chrom}_NCPG{NCPG}'.format(mms = '.'.join(map(str,mm_samples)), chrom = chrom, NCPG = NCPG) 
		"""
		input_f_list = []
		for el in el_samples:
			input_f = filespath+'binned_winsize_1000_As_RG.E.maximus.' + el + '_' + chrom + '_' + NCPG + '.' + chrom + '.F'
			input_f_list.append(input_f)
		out_f_handle = 'Mean_binned_winsize_1000_As_RG.E.maximus.{els}_{chrom}_NCPG{NCPG}'.format(els = '.'.join(map(str,el_samples)), chrom = chrom, NCPG = NCPG) 
		"""
		merge_dfs(input_f_list,out_f_handle)

#f_handles = ['binned_winsize_1000_As_RG.E.maximus.T13_NC_064837.1_25.NC_064837.1.F','binned_winsize_1000_As_RG.E.maximus.T16_NC_064837.1_25.NC_064837.1.F','binned_winsize_1000_As_RG.E.maximus.T18_NC_064837.1_25.NC_064837.1.F']
#merge_dfs(f_handles,'MeanF_test')

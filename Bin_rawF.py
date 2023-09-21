#!/usr/bin/python3
from __future__ import division
import pandas as pd
import numpy as np
import itertools
import sys

"""This script is used to bin the DamMet raw F results by the specified window size (default 1000bp)
Input: 
	sample_name (str): to get the full names of DamMet raw F files 
	window_size: (int) default = 1000
	chromsome_list_file: filepath to a txt file containing all chromsome names 
Outout:
	tsv files: "binned_windowsize_"
"""
def bin_raw_f_file(input_f,window_size=1000):
	"""This function is to calculate f rate means for each window in a F file
	Input:
		input_f (str): The file path of the F file
		output_f (str): The file path for the output csv file
		window_size (int): The window size for taking the average f rate (default = 10000)

	Output:
		a pd dataframe with the two columns: window location and average f rate 
		(saved as a tsv file)
	"""
	#read F file into dataframe
	raw_df = pd.read_table(input_f, sep = ' ', dtype = 'str', usecols = [0,5])

	#formatting the location in chromosome into integer, mythlation rate as float
	raw_df['center'] = raw_df['center'].str.split(':').str.get(-1).astype(int)
	#print(raw_df.head(5))
	
	raw_df['f'] = raw_df['f'].astype(float)

	#get the furthest chromosome coordinates with f rate (last one in center column)
	n = raw_df['center'].iloc[-1]

	#create a dataframe to record windows and their respective average methylation rate
	win_df = pd.DataFrame(columns = ['window_location','average_f_rate'])

	#move the window with specified window size, calculate the avergae methylation rate within the current window
	#probably can start with 1
	for i in range(1,n+window_size,window_size):
		if i <= n:
			one_window = raw_df[(raw_df['center'] >=i) & (raw_df['center'] < i+window_size)]
			m_site_counter = np.count_nonzero(one_window['f'],axis = 0)
			if m_site_counter != 0:
				avg_rate = one_window.loc[:,'f'].sum()/float(m_site_counter)
				win_df.loc[len(win_df.index)] = ['{start}-{end}'.format(start = i, end = i+window_size),avg_rate]			
			else:
				avg_rate = 0
				win_df.loc[len(win_df.index)] = ['{start}-{end}'.format(start = i, end = i+window_size),avg_rate]
		else:
			one_window = raw_df[(raw_df['center'] >=i) & (raw_df['center'] <= n)]
			m_site_counter = np.count_nonzero(one_window['f'],axis = 0)
			if m_site_counter!= 0:
				avg_rate = one_window.loc[:,'f'].sum()/float(m_site_counter)
				win_df.loc[len(win_df.index)] = ['{start}-{end}'.format(start = i, end = n),avg_rate]
			else:
				avg_rate = 0
				win_df.loc[len(win_df.index)] = ['{start}-{end}'.format(start = i, end = n),avg_rate]

	#write the dataframe into the specified csv file path
	filepath = "/crex/proj/snic2022-6-144/nobackup/YUXIN/DamMet/binned_F/"
	win_df.to_csv(filepath+"binned"+ "_winsize_"+str(window_size)+'_'+input_f, sep = "\t", index = False)
	return win_df

#raw -> bin F files
#specify sample names 
sample_name = sys.argv[1]

#fetch chromsomes lists 
chrom_list = []
with open("/crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/chromsomes_list.txt", 'r') as f0:
	for line in f0:
		chrom_list.append(line.strip())

#run whole genomes, all NCPG parameters 
NCPGs = ["10","25","50","100"]

for chrom in chrom_list:
	for NCPG in NCPGs:
		full_f_name = 'As_'+sample_name +'_'+ chrom + '_' + NCPG +'.'+chrom+'.F'
		bin_raw_f_file(full_f_name,1000)

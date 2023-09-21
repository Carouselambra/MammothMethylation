#!/usr/bin/python3
from __future__ import division
import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import itertools
import sys
from scipy import stats
from matplotlib.patches import FancyArrowPatch
import matplotlib.transforms as transforms
import matplotlib.ticker as mtick


"""This script is used to analyze candidate results (output of id_rank_genextrap.py)
by providing:
	mean, standard deviation
	distribution plot (histogram)
	filtered_results (only keep the outliers, outside of mean +- 2std)

Input:
 	Tsv files with columns ("Diff_in_f_mean_abs","Diff_in_f_mean","Target_Info")

Output:
	Distribution plots (historigram) on "Diff_in_f_mean" of the input tsv file
	Filtered tsv files (contain outliers of the input files)

Functions:
	read_candidates(file_handle)
	get_statistics(candidates_df)
	plot_distr(candidates_df,plot_name)
	first_filter_result(candidates_df,out_name,manual_threshold = "off")

"""
def read_candidates(file_handle):
	"""This function is to read tsv files as pandas dataframe
	"""
	df = pd.read_csv(file_handle,sep = '\t',usecols = ["Diff_in_f_mean_abs","Diff_in_f_mean","start","end","chrom","Target_Info"])
	
	#remove sex chromosomes and no chromsome associated scaffold regions 
	mapped_mask = df['chrom'].str.contains('NC_')
	autosome_mask = (df['chrom'] != 'NC_064846.1') & (df['chrom'] != 'NC_064847.1')
	df = df[mapped_mask].copy()
	df = df[autosome_mask].copy()

	df.reset_index()
	df.sort_values(by = "Diff_in_f_mean_abs", ascending = False,inplace = True)
	return df

def get_statistics(candidates_df):
	"""This function is used to calculate the mean and the standard deviation of the result dataset
	input:
		candidate_df: a pandas dataframe contains columns ("Diff_in_f_mean_abs","Diff_in_f_mean","Target_Info")
	output:
		mean and standard deviation
	"""
	mean = candidates_df['Diff_in_f_mean'].mean()
	std = candidates_df['Diff_in_f_mean'].std()
	#print(mean, std)
	return mean, std

def plot_single_distr(candidates_df,plot_name):
	"""This function is used to plot distribution of candidate area (gene/cis) \
	differences in f rate
	input: candidates_df
	output: a histogram 
	"""
	mean, std = get_statistics(candidates_df)

	fig,ax = plt.subplots(figsize = (12,9))

	candidates_df['Diff_in_f_mean'].hist(bins = 100,color = "#274060")

	ax.set_xlabel('Difference in area f mean ')
	ax.set_ylabel('frequency counts')

	#vertical line show mean
	ax.axvline(x =mean,color = '#BF4342',linestyle = 'dashed')

	#verical lines show mean +- std, +- 2std
	ax.axvline(x = mean + std, alpha=0.9,color = '#F17300',linestyle = 'dashed')
	ax.axvline(x = mean + 2*std, alpha=0.9,color = '#F17300',linestyle = 'dashed')
	ax.axvline(x = mean - std, alpha=0.9,color = '#F17300',linestyle = 'dashed')
	ax.axvline(x = mean - 2*std, alpha=0.9,color = '#F17300',linestyle = 'dashed')

	plt.show()
	plt.savefig("Plot_"+plot_name+".jpeg")

def first_filter_result(candidates_df,out_name,manual_threshold = "off"):
	"""
	This function is to filter the results by the absolute value of differences in area F mean 
	by default filtered by 2 conditions
	1. cut off by mean +- 2*std
	2. get rid of pseduo gene
	"""
	mean, std = get_statistics(candidates_df)
	if manual_threshold == "off":
		threshold_up = mean + 2*std
		threshold_down = mean - 2*std
		candidates_df_1 = candidates_df[candidates_df["Diff_in_f_mean"] > threshold_up]
		candidates_df_2 = candidates_df[candidates_df["Diff_in_f_mean"] < threshold_down]
		print(threshold_down)
		print(threshold_up)
	#else:
	#	candidates_df = candidates_df[candidates_df["Diff_in_f_mean_abs"] > threshold]
	else:
		threshold_up = manual_threshold[0]
		threshold_down = manual_threshold[1]
		candidates_df_1 = candidates_df[candidates_df["Diff_in_f_mean"] > threshold_up]
		candidates_df_2 = candidates_df[candidates_df["Diff_in_f_mean"] < threshold_down]

	filtered_df = pd.concat([candidates_df_1,candidates_df_2])

	filtered_df = filtered_df[filtered_df['Target_Info'].str.contains("gene_biotype=pseudogene") == False]
	filtered_df.sort_values(by = "Diff_in_f_mean_abs",inplace = True,ascending = False)
	filtered_df.reset_index(drop = True)
	filtered_df.to_csv(out_name, sep = '\t',index = False)
	return filtered_df

def plot_all_distr(candidates_df_list,labels,mode):
	"""This function is used to plot distribution of candidate area (gene/cis) \
	differences in f rate
	input: 
		candidates_df_list
		labels: comparsion names
		mode: regulatory regions or genes 
	output: a figure with four subplot
	"""

	fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1,figsize = (12,9), sharex = True, sharey = True, dpi = 300)
	ax_list = [ax1,ax2,ax3,ax4]
	colors = ["#320D6D","#628395","#CF995F","#DBAD6A"]
	#print(labels[0])
	for i in range(len(candidates_df_list)):
		df = candidates_df_list[i]
		mean, std = get_statistics(df)

		ax_list[i].hist(df['Diff_in_f_mean'], weights = np.ones_like(df['Diff_in_f_mean']) / len(df['Diff_in_f_mean']), bins = 100, color = colors[i],label =labels[i])
		ax_list[i].legend()
		ax_list[i].set_ylabel('frequency counts')
		trans = ax_list[i].get_xaxis_transform()
		ax_list[i].yaxis.set_major_formatter(mtick.PercentFormatter(xmax = 1.0, decimals = 0, symbol = r'%', is_latex = True))

		#vertical line show mean
		ax_list[i].axvline(x = mean,color = '#BF4342',linestyle = 'dashed')

		#verical lines show mean +- std, +- 2std
		ax_list[i].axvline(x = mean + std, alpha=0.9,color = '#F17300',linestyle = 'dashed')
		ax_list[i].axvline(x = mean + 2*std, alpha=0.9,color = '#F17300',linestyle = 'dashed')
		ax_list[i].axvline(x = mean - std, alpha=0.9,color = '#F17300',linestyle = 'dashed')
		ax_list[i].axvline(x = mean - 2*std, alpha=0.9,color = '#F17300',linestyle = 'dashed')

	ax_list[-1].set_xlabel('Difference in window-specific f fold change')

	plt.savefig('./final_figs/id_{target}_distribution.jpeg'.format(target = mode))

	plt.show()


#local test - all sample distribution plots, promoters and genes 
path = './id_result/'
pairwise_compare = ["M.primigenius.E467.L163.L164.P005.M6_vs_E.maximus.T13.T16.T18",\
"M.primigenius.E467.L163.L164.P005.M6_vs_L.africana.T11",\
"M.primigenius.E467.L163.L164.P005.M6_vs_P.antiquus.ERR2260504",\
"M.primigenius.E467.L163.L164.P005.M6_vs_M.americanum.ERR2260503"]
candidates_df_list = []
#mode = "Genes" #either Cis or Genes
mode = "Cis"
labels = ["M.primigenius_vs_E.maximus","M.primigenius_vs_L.africana","M.primigenius_vs_P.antiquus","M.primigenius_vs_M.americanum"]
for c in pairwise_compare:
	c_handle = path + "Candidates_"+mode+'_'+c+'_NCPG25'
	out_handle = path + "Filtered_Candidates_"+mode+'_'+c+'_NCPG25'
	df = read_candidates(c_handle)
	first_filter_result(df,out_handle)
	candidates_df_list.append(df)


plot_all_distr(candidates_df_list,labels,mode)
def gaussian_test(col, values):
    stat1, p1 = stats.shapiro(values)
    stat2, p2 = stats.normaltest(values)

    print(f"Gaussian: {col}\n\t{p1:5f} (Shapiro-Wilk)\n\t{p2:5f} (D'Agostino's)")


#normality test
"""
col = 'MMvsASE'
values = candidates_df_list[1]['Diff_in_f_mean']
gaussian_test(col,values)
"""



#run on the cluster
#get filtered results (probably need to rerun this! To get rid of sex chromsomes and scafold genes/cis)
"""
comparison_check = sys.argv[1]
comparison_name = ""

if comparison_check == "MMs_vs_ASEs":
	#Mammoths vs. Asian Elephants candidate list
	comparison_name = "M.primigenius.E467.L163.L164.P005.M6_vs_E.maximus.T13.T16.T18"

elif comparison_check == "MMs_vs_AFE":
	#Mammoths vs. African elephant T11
	comparison_name = "M.primigenius.E467.L163.L164.P005.M6_vs_L.africana.T11"

elif comparison_check == "MMs_vs_STKS":
	#Mammoths vs. straight tussk elephant sample
	comparison_name = "M.primigenius.E467.L163.L164.P005.M6_vs_P.antiquus.ERR2260504"

elif comparison_check == "MMs_vs_MSTD":
	#Mammoths vs. mastodon sample
	comparison_name = "M.primigenius.E467.L163.L164.P005.M6_vs_M.americanum.ERR2260503"

NCPGs = ['25','50','100','10']
for NCPG in NCPGs:
	cis_handle = "Candidates_Cis_"+comparison_name+'_NCPG'+NCPG
	gene_handle = "Candidates_Genes_"+comparison_name+'_NCPG'+NCPG

	cis_df = read_candidates(cis_handle)
	print('cis_df converted!')
	gene_df = read_candidates(gene_handle)
	print('gene_df converted!')

	plot_distr(cis_df,cis_handle)
	print('cis distr plot drawn!')
	plot_distr(gene_df,gene_handle)
	print('gene distr plot drawn!')

	first_filter_result(cis_df,cis_handle)
	print('filtered cis!')
	first_filter_result(gene_df,gene_handle)
	print('filtered gene!')
	print(NCPG + "DONE!")
"""
#local test

"""
file_handle_cis = "Candidates_Cis_M.primigenius.E467.L163.L164.P005.M6_vs_E.maximus.T13.T16.T18_NCPG25"
file_handle_genes = "Candidates_Genes_M.primigenius.E467.L163.L164.P005.M6_vs_E.maximus.T13.T16.T18_NCPG25"

candidates_df_1 = read_candidates(file_handle_cis)
plot_single_distr(candidates_df_1,file_handle_cis)
filtered_df_1 = first_filter_result(candidates_df_1,file_handle_cis)
#search_by_id(filtered_df_1,search_id)
"""



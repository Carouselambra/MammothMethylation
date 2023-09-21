#!/usr/bin/python3
from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import itertools
import sys
from scipy import stats
from matplotlib.patches import FancyArrowPatch
from matplotlib_venn import venn2

"""This script is used to analyze the filtered candidate target (gene / cis) list
"""
def read_filtered_candidates(filtered_handle):
	filtered_df = pd.read_csv(filtered_handle,sep = '\t',usecols=["Diff_in_f_mean_abs","Diff_in_f_mean","start","end","chrom","Target_Info"])
	filtered_df = filtered_df[filtered_df['Target_Info'].str.contains('description=')].copy()
	filtered_df['description'] = filtered_df['Target_Info'].apply(lambda x: str(x.partition('description=')[-1].partition(';')[0]))
	filtered_df['GeneID'] = filtered_df['Target_Info'].apply(lambda x: str(x.partition('GeneID:')[-1].partition(';')[0]))

	return filtered_df

def analyze_description(filtered_df):
	#check how many have description 
	d = filtered_df['Target_Info'].str.contains("description").sum()
	missing_d = len(filtered_df.index) - d
	
	working_df = filtered_df[filtered_df['Target_Info'].str.contains("description") == True].copy()
	working_df['description'] = working_df['Target_Info'].apply(lambda x: str(x.partition('description=')[-1].partition(';')[0]))
	

	print(working_df['description'].unique())
	return list(working_df['description'].unique())

def reverse_filter(filtered_df,info_list):
	working_df = filtered_df[filtered_df['Target_Info'].str.contains("description") == True].copy()
	working_df['description'] = working_df['Target_Info'].apply(lambda x: str(x.partition('description=')[-1].partition(';')[0]))

def search_info_occurance(filtered_df,info, mode = "short"):
	"""This function is to parse the candidate df by the description in 'Target into'
	"""
	mask = filtered_df['Target_Info'].str.contains(info)
	target = filtered_df[mask][['Diff_in_f_mean','description','GeneID']]
	m = mask.sum()
	n = len(filtered_df.index)
	#print("{Info} has total {m} occurance out of {n}".format(Info = info, m = m, n = n))
	if mode == "short":
		return m
	elif mode == "extensive":
		return m, target

def info_pie(filtered_df):
	
	#olfactory 
	or_count, total_count = search_info_occurance(filtered_df,"olfactory receptor")
	
	#transfer RNA 
	tRNA_count = search_info_occurance(filtered_df,"transfer RNA")[0]
	
	#T cell receptor 
	tcell_count = search_info_occurance(filtered_df,"T cell receptor")[0]
	
	#histocompatibility related 
	#histoco_count = search_info_occurance(filtered_df,"histocompatibility")[0]

	#spliceosomal RNA
	#splice_count = search_info_occurance(filtered_df,"spliceosomal RNA")[0]

	#immunoglobulin
	immuglo_count = search_info_occurance(filtered_df,"immunoglobulin")[0]

	#nkanti_count = search_info_occurance(filtered_df,"natural killer cells antigen")[0]
	
	#zincf_count = search_info_occurance(filtered_df,"zinc finger")[0]
	
	#ribon_count = search_info_occurance(filtered_df,"ribonuclease")[0]

	#uncharacterized gene 
	unchar_count = search_info_occurance(filtered_df,"uncharacterized")[0]

 
	values = [or_count,tRNA_count,tcell_count,immuglo_count,unchar_count]
	other_count = total_count
	for count in values:
		#print(count)
		other_count = other_count - count
		print(other_count)
	values.append(other_count)
	print(values)
	names ="olfactory receptors","transfer RNA","T cell receptor","immunoglobulin related","uncharacterized","Others"

	plt.pie(values,labels = names, labeldistance=1.15,autopct='%1.2f%%')
	plt.show()

def search_by_info(filtered_df,info):
	filtered_df['description'] = filtered_df['Target_Info'].apply(lambda x: str(x.partition('description=')[-1].partition(';')[0]))
	result = filtered_df[filtered_df['description'].str.contains(info)]

	#print(result)

def plot_gene_on_binF(filtered_df,mode = "cis", plot_size = 100, bin_size = 1000,slide_window = 10):
	"""This function is to plot candidate area(promoter or gene) alongside the bin F graph
	input:
		candidate df columns =  ["Diff_in_f_mean_abs","Diff_in_f_mean","start","end",'chrom',"Target_Info"]
	output:
		normalised F files 
	"""
	chrom_list = ["NC_064828.1"]
	for chrom in chrom_list:
		#read normalised F files by chromosomes 
		nF_1_handle = "normlised_averaged_RG.M.primigenius.E467.M.primigenius.L163.M.primigenius.L164.M.primigenius.P005.M.primigenius.M6_NC_064828.1_25"
		nF_2_handle = "normlised_averaged_RG.E.maximus.T13.E.maximus.T16.E.maximus.T18_NC_064828.1_25"
		nF_df_1 = pd.read_csv(nF_1_handle,sep = '\t',usecols= ['window_location','normalised_f_fold'])
		nF_df_2 = pd.read_csv(nF_2_handle,sep = '\t',usecols= ['window_location','normalised_f_fold'])

		#nF_df_1['start'] = nF_df_1['window_location'].partition('-')[0]
		#nF_df_1['end'] = nF_df_1['window_location'].partition('-')[-1]

		#nF_df_2['start'] = nF_df_2['window_location'].partition('-')[0]
		#nF_df_2['end'] = nF_df_2['window_location'].partition('-')[-1]


		this_chrom_filtered_df = filtered_df[filtered_df['chrom'] == chrom].copy()
		n = len(this_chrom_filtered_df)
		print(this_chrom_filtered_df.head(5))

		#iterate through target region (promoter or gene)
		for i in range(n):
		#for i in range(n): 
			start = int(this_chrom_filtered_df['start'].iloc[i])
			end = int(this_chrom_filtered_df['end'].iloc[i])
			print(end-start)

			#if end-start >= 200*bin_size:
			if end != 0:
				gene_name = this_chrom_filtered_df['Target_Info'].iloc[i]
				#gene / promoter region annotation 
				#annotation_start = (start // bin_size) * plot_size
				#annotation_end = (end // bin_size) * plot_size

				plot_start = (start // plot_size - slide_window) * plot_size #1000 is set because the window size of binning step
				plot_end = plot_start + bin_size*plot_size


				#subset normalised F dataframes, windows has to be within plot_start and plot_end
				mask_1 = (nF_df_1['window_location'].str.partition('-')[0].astype('int') >= plot_start) & (nF_df_1['window_location'].str.partition('-')[2].astype('int') <= plot_end)
				mask_2 = (nF_df_2['window_location'].str.partition('-')[0].astype('int') >= plot_start) & (nF_df_2['window_location'].str.partition('-')[2].astype('int') <= plot_end)

				
				plot_df1 = nF_df_1[mask_1].copy().reset_index(drop = True)
				plot_df2 = nF_df_2[mask_2].copy().reset_index(drop = True)
				#print(len(plot_df1.index),len(plot_df2.index))
				

				#obtain the relative X location 


				#get the differences between the two normalised F dataframes 
				highlight_df = plot_df1.merge(plot_df2, how = "inner", on = "window_location", suffixes = ("_1","_2"))
		
				if highlight_df.empty == False:
					highlight_df_1 = highlight_df[highlight_df['normalised_f_fold_1'] > highlight_df['normalised_f_fold_2']].copy()
					highlight_df_2 = highlight_df[highlight_df['normalised_f_fold_2'] > highlight_df['normalised_f_fold_1']].copy()

					#convert them into a list 
					x_highlights_1 = highlight_df_1.index.values.tolist()
					x_highlights_2 = highlight_df_2.index.values.tolist()

				#obtain the relative X location of the target area
				annotated_condition = (highlight_df['window_location'].str.partition('-')[0].astype('int') >= start) & (highlight_df['window_location'].str.partition('-')[2].astype('int') <= end)
				annotated_arrow_xs = highlight_df[annotated_condition].index.values.tolist()
				print(annotated_arrow_xs)

				
				#set up x axis ticks, labels
				x_ = np.arange(0,plot_size,1).tolist()
				xticks = [1,int(plot_size*0.2),int(plot_size*0.4),int(plot_size*0.6),int(plot_size*0.8),plot_size]
				xticks_label = []
				for j in xticks:
						if j <= len(plot_df1.index):
							x_label = "{:.2E}".format(int(plot_df1.iloc[j-1]['window_location'].partition('-')[0]))
							xticks_label.append(x_label)
						else:
							xticks_label.append(" ")

			
				fig,(ax1,ax2) = plt.subplots(2,1,figsize = (16,9))

				#plot the normalised F 
				plt.subplot(2,1,1)
				highlight_df.plot(y = 'normalised_f_fold_1', kind = 'bar', width = 1, ax = ax1, color = '#F4AC45')
				#ax1.legend([label[0]])
				ax1.set_ylim(0,6)
				ax1.set_ylabel('F fold change')
				#ax1.set_xlabel(label[2],loc='right')
				ax1.set_xticks(xticks,xticks_label,rotation = 0)

				plt.subplot(2,1,2)
				highlight_df.plot(y = 'normalised_f_fold_2', kind = 'bar', width = 1, ax = ax2, color = '#4281A4')
				#ax2.legend(label[1])
				ax2.set_ylim(0,6)
				ax2.set_ylabel('F fold change')
				#ax2.set_xlabel(label[2],loc='right')
				ax2.set_xticks(xticks,xticks_label,rotation = 0)

				#plot the highlighted differences
				
				#for highlight in x_highlights_1:
				#	ax1.axvline(x=highlight,alpha=0.5,color = '#F56476')
				#	#ax1.text(highlight+0.1,0.1,'threshold = {num}'.format(num = threshold))
				#	ax2.axvline(x=highlight,alpha=0.5,color = '#F56476')
					#ax2.text(highlight+0.1,0.1,'threshold = {num}'.format(num = threshold))

				#for highlight in x_highlights_2:
				#	ax1.axvline(x=highlight,alpha=0.5,color = '#2A7F62')
				#	#ax1.text(highlight+0.1,0.1,'threshold = {num}'.format(num = threshold))
				#	ax2.axvline(x=highlight,alpha=0.5,color = '#2A7F62')
				
				#plot annotated gene or promoter area
				#ax1.annotate('', xy=(0, -0.1), xycoords=ax1, xytext=(1, -0.1),
				#arrowprops=dict(arrowstyle="<->", color='b'))

				#specify start and end point of target area
				#draw highlighted area instead of the arrows
				arrow_centers = ((annotated_arrow_xs[0]-1,2),(annotated_arrow_xs[-1]+1,2))
				target_area = FancyArrowPatch(*arrow_centers)
				target_area_2 = FancyArrowPatch(*arrow_centers)
				ax1.add_patch(target_area)
				ax1.annotate(gene_name, xy=(0,3), xycoords=ax1,xytext=(1, 1))
				
				ax2.add_patch(target_area_2)
				#ax2.annotate(gene_name, xy=(0,3), xycoords=ax1, xytext=(1, -0.1))
				plt.show()
			
def further_outliers(filtered_df,threshold):

	outlier_df = filtered_df[filtered_df['Diff_in_f_mean_abs'] > threshold].copy()
	outlier_df['description'] = outlier_df['Target_Info'].apply(lambda x: str(x.partition('description=')[-1].partition(';')[0]))
	print(outlier_df['description'].unique())
	print(outlier_df)
	#outlier_df.hist(bins = 10,color = "#274060")
	#plt.show()
	return outlier_df

def parse_for_panther(filtered_df,outhandle):
	"""This function is used to converted filtered df to a list of genes 
	for GO enrichment analysis 
	"""
	working_df = filtered_df[filtered_df['Target_Info'].str.contains("Dbxref=GeneID:") == True].copy()
	working_df['Entrez_ID'] = working_df['Target_Info'].apply(lambda x: str(x.partition('Dbxref=')[-1].partition(';')[0]))
	working_df.to_csv(outhandle,sep = '\t',index = False,header = False,columns = ['Entrez_ID'])
	print(working_df.head(5))

def filter_descriptor_counts(filtered_df,search_terms,mode):
	"""This function is used to find the frequency counts of major descriptors  
	input: 
		filtered_df 
		mode: ASE or AFE
	output:
		descriptor_count_df
	"""
	descriptor_count_dict = {}
	n = len(filtered_df.index)
	#print(filtered_df[filtered_df['description'].str.contains("small proline-rich protein")])
	occurance_df = filtered_df.copy()
	for term in search_terms:
		count = search_info_occurance(occurance_df,term)
		print(count, term)

		#descriptor_count_dict[term]
		filtered_df = filtered_df[filtered_df['description'].str.contains(term) == False]

		#result_df = filtered_df[filtered_df['description'].str.contains(term)]
		#print(result_df)
	filtered_df.to_csv('./id_result/z_further_check_{mode}_'.format(mode = mode),sep = '\t',index = False,columns = ["GeneID","description","Diff_in_f_mean_abs","Diff_in_f_mean"])
	
	#print(filtered_df['description'].unique)
	#descriptor_count_dict = {}
	#for term in search_terms:
	#	count = search_info_occurance(filtered_df,term)
	#	descriptor_count_dict[term] = int(count)
	#print(descriptor_count_dict)

def polish_category(filtered_df,category_index,suffix):
	"""This function is to filtered through description columns 
	to find ones were not in category_index
	"""
	others_count = len(filtered_df.index)
	for big_term, sub_terms_dict in category_index.items():
		for sub_term, count in sub_terms_dict.items():
			count = search_info_occurance(filtered_df,sub_term)
			filtered_df = filtered_df[filtered_df['description'].str.contains(sub_term) == False]
			sub_terms_dict[sub_term] = count
			others_count -= count
	
	category_index['others']['others'] = others_count
	filtered_df.to_csv('./id_result/z_further_check_{suffix}'.format(suffix = suffix),sep = '\t',index = False,columns = ["GeneID","description","Diff_in_f_mean_abs","Diff_in_f_mean"])

	for big_term, sub_terms_dict in category_index.items():
		for sub_term, count in list(sub_terms_dict.items()):
			if count == 0:
				del sub_terms_dict[sub_term]
	for big_term, sub_terms_dict in category_index.items():
		for sub_term, count in sub_terms_dict.items():
			print(sub_term,count)
	return category_index

def count_filtered_by_category_index(filtered_df,category_index):
	"""This function is to filtered through the filitered_df 
	and count the frequency of the corresponding values of each description 
	input:
		filitered_df
		category_index a nested dictionary
	output:
		category_index: a nested dictionary with updated frequency count
	"""
	total_count = len(filtered_df.index)
	for big_term, sub_terms_dict in category_index.items():
		for sub_term, count in sub_terms_dict.items():
			count, details = search_info_occurance(filtered_df,sub_term, mode = "extensive")
			total_count -= count
			sub_terms_dict[sub_term] = count

	category_index['others']['others'] = total_count
	#for big_term, sub_terms_dict in category_index.items():
		#for sub_term, count in sub_terms_dict.items():
			#print(big_term,sub_term,count)
	return category_index

def plot_lollipop(filtered_df,category_index,mode = "simple"):
	"""This function is to group the filtered candidate targets (genes, cis)
	then plot them into a lollipop plot 
	input:
		filtered_df
		mode: "board"
	"""
	plot_dict = {}
	if mode == "simple":
		for big_term, sub_terms_dict in category_index.items():
			plot_dict[big_term] = 0
			big_term_count = 0
			for sub_term, count in sub_terms_dict.items():
				big_term_count += count
			plot_dict[big_term] = big_term_count
		plot_df = pd.DataFrame.from_dict(plot_dict, orient = "index",columns = ["count"])
		plot_df = plot_df.rename_axis('category').reset_index()
		print(plot_df)
		fig,ax = plt.subplots(1,1)
		ax.stem(plot_df['count'])
		ax.set_xticks(np.arange(len(plot_df['category'].index)))
		ax.set_xticklabels(plot_df['category'],rotation = 45)
		plt.show()

		#missing first (0) label???
		#ax.set_xticklabels(plot_df['category'].to_list())
		#plt.show()

def plot_single_target(target_dict,species2):
	"""this plot should take one gene / regulatory region coordinates 
	then plot their background 
	input:
		target_dict: contains keys "chromosomes", "start", "end","GeneID"

	output:
		a plot with two subplot 
	"""
	chrom = target_dict['chrom']
	target_start = int(target_dict['start'])
	target_end = int(target_dict['end'])
	id_text = target_dict['GeneID']

	if target_end - target_start == 5000:
		plot_mode = "promoter"
		#plot_size = 100
	elif target_end - target_start > 5000:
		plot_mode = "large_gene"
	elif target_end - target_start < 5000:
		plot_mode = "small_gene"
		#plot_size = 1000

	#find the first and the last window of the highlighted regions
	target_first_window_start = target_start // 1000 * 1000 +1
	target_last_window_start = target_end // 1000 * 1000+1

	target_first_window = "{start}-{end}".format(start = str(target_first_window_start), end = str(target_first_window_start+1000))
	target_last_window = "{start}-{end}".format(start = str(target_last_window_start), end = str(target_last_window_start+1000))

	#read the normalised F files 
	path = './normalised_Fs/'
	nF_1_handle = path + "normlised_averaged_RG.M.primigenius.E467.M.primigenius.L163.M.primigenius.L164.M.primigenius.P005.M.primigenius.M6_{chrom_name}_25".format(chrom_name = chrom)
	
	if species2 == "ASE":
		label_2 = "E.maximus"
		nF_2_handle = path + "normlised_averaged_RG.E.maximus.T13.E.maximus.T16.E.maximus.T18_{chrom_name}_25".format(chrom_name = chrom)
	elif species2 == "AFE":
		label_2 = "L.africana"
		nF_2_handle = path + "normalised_binned_winsize_1000_As_RG.L.africana.T11_{chrom_name}_25.{chrom_name}.F".format(chrom_name = chrom)
	elif species2 == "MSTD":
		label_2 = "M.americanum"
		nF_2_handle = path + "normalised_binned_winsize_1000_As_RG.M.americanum.ERR2260503_{chrom_name}_25_{chrom_name}.F".format(chrom_name = chrom)
	elif species2 == "STSK":
		label_2 = "P.antiquus"
		nF_2_handle = path + "normalised_binned_winsize_1000_As_RG.P.antiquus.ERR2260504_{chrom_name}_25_{chrom_name}.F".format(chrom_name = chrom)

	nF_df_1 = pd.read_csv(nF_1_handle,sep = '\t',usecols= ['window_location','normalised_f_fold'])
	nF_df_2 = pd.read_csv(nF_2_handle,sep = '\t',usecols= ['window_location','normalised_f_fold'])
	
	target_index_start_1 = nF_df_1[nF_df_1['window_location'] == target_first_window].index.values.astype(int)[0]
	target_index_start_2 = nF_df_2[nF_df_2['window_location'] == target_first_window].index.values.astype(int)[0]
	#print(target_index_start_1,target_index_start_2)


	target_index_end_1 = nF_df_1[nF_df_1['window_location'] == target_last_window].index.values.astype(int)[0]
	target_index_end_2 = nF_df_2[nF_df_2['window_location'] == target_last_window].index.values.astype(int)[0]

	#check if the two normalised F is parallel or not 
	if target_index_start_1 != target_index_start_2:
		print("your index doesn't match!!!")

	if plot_mode == "promoter":
		subplot1_index_start = target_index_start_1 - 25
		subplot1_index_end = target_index_end_1 + 70
	
		subplot2_index_start = target_index_start_2 - 25
		subplot2_index_end = target_index_end_2 + 70

	elif plot_mode == "large_gene":
		subplot1_index_start = target_index_start_1 - 50
		subplot1_index_end = target_index_end_1 + 50
	
		subplot2_index_start = target_index_start_2 - 50
		subplot2_index_end = target_index_end_2 + 50
	
	elif plot_mode == "small_gene":
		subplot1_index_start = target_index_start_1 - 25
		subplot1_index_end = subplot1_index_start + 99

		subplot2_index_start = target_index_start_2 - 25
		subplot2_index_end = target_index_end_2 + 99


	else:	
		print("specify the mode, either 'promoter' or 'gene'")
	print(plot_mode)
	#elif in gene mode (-200 +300)


	plot_df_1 = nF_df_1.iloc[subplot1_index_start:subplot1_index_end].copy().reset_index()
	plot_df_2 = nF_df_2.iloc[subplot2_index_start:subplot2_index_end].copy().reset_index()

	highlight1_start = plot_df_1[plot_df_1['window_location'] == target_first_window].index.values.astype(int)[0]
	highlight1_end = plot_df_1[plot_df_1['window_location'] == target_last_window].index.values.astype(int)[0]
	#print(highlight1_start,highlight1_end)
	highlight2_start = plot_df_2[plot_df_2['window_location'] == target_first_window].index.values.astype(int)[0]
	highlight2_end = plot_df_2[plot_df_2['window_location'] == target_last_window].index.values.astype(int)[0]

	plot_df_1['window_start'] = plot_df_1['window_location'].str.partition('-')[0]
	plot_df_2['window_start'] = plot_df_2['window_location'].str.partition('-')[0]
	#print(len(plot_df_1.index))
	#label_list = plot_df_2['window_start']

		#xticks_label1 = plot_df_1['window_location'].iloc[0,25,50,75,100]
		#print(xticks_label1)
	#elif plot_mode == "gene":


	fig, (ax1,ax2) = plt.subplots(2,1, figsize = (12,4),sharex = True, sharey = True)
	fig.canvas.draw()
	#plot selected background 
	ax1.bar(plot_df_1.index,plot_df_1['normalised_f_fold'],label = 'M.primigenius')
	ax1.legend()

	#plot target region at subplot 1
	ax1.axvspan(highlight1_start,highlight1_end,alpha = 0.5)

	#plot the normalised average line 
	ax1.axhline(y = 1, linestyle =  '--', color = "#000000")

	#plot the background at subplot 2
	ax2.bar(plot_df_2.index,plot_df_2['normalised_f_fold'],color = "#8EA604",label = label_2)
	ax2.legend()
	#plot target region at subplot 2 
	ax2.axvspan(highlight2_start,highlight2_end,alpha = 0.5)

	#plot the normalised average line 
	ax2.axhline(y = 1, linestyle =  '--', color = "#000000")

	# labeled as promoter?
	if plot_mode == "promoter" or plot_mode == "small_gene":
		x_ticks_list = [0,19,39,59,79,99]

	elif plot_mode == "large_gene":
		x_ticks_total = len(plot_df_1.index)
		x_ticks_list = [0,int(x_ticks_total*0.2)-1,int(x_ticks_total*0.4)-1,int(x_ticks_total*0.6)-1,int(x_ticks_total*0.8)-1,int(x_ticks_total-1)]
	ax2.set_xticks(x_ticks_list)
	ax2.set_xticklabels(list(plot_df_2['window_start'][x_ticks_list]))


	#ax2.set_xticklabels(plot_df_2['window_start'])
	ax2.xaxis.set_tick_params(labelsize=9)
	
	#adding chromosome to bottom left side 
	ax2.set_xlabel(chrom + ':' + '\n' + 'coordinates')
	ax2.xaxis.set_label_coords(-0.07, -0.06)

	plt.show()

def plot_grouped_bar(filtered_df_list,category_index):
	"""
	"""
	series_df_list = []
	for filtered_df in filtered_df_list:
		big_category_dict = {"olfactory receptor":{"olfactory receptor":0,"trace amine-associated receptor 8":0,"gustatory receptor clone PTE03":0},\
		"immune system related":{"butyrophilin subfamily 1 member A1":0,"MHC class I":0,"protein FAM90A27P":0,"T cell receptor":0,"natural killer cells antigen":0,"immunoglobulin":0,"histocompatibility antigen":0,\
		"small proline-rich protein":0,"defensin":0,"interferon":0,"guanylate-binding protein":0, "E3 ubiquitin-protein ligase":0, "schlafen-like protein":0,\
		"NKG2-A/NKG2-B type II integral membrane protein":0,"copper transport protein ATOX1":0,"tripartite motif":0,\
		"NLR family pyrin domain containing":0,"toll-like receptor":0,"protein S100-A7":0,"TNF receptor superfamily member 17":0,\
		"late cornified envelope-like proline-rich protein":0,"basic salivary proline-rich protein 1":0,"RHO family interacting cell polarization regulator 2":0,\
		"T-cell receptor-associated transmembrane adapter":0,"interleukin":0},\
		"non coding RNA related":{'transfer RNA':0,"small nucleolar":0,"spliceosomal":0,"ribonuclease":0,"transcription initiation factor TFIID subunit 4":0,"small nuclear RNA":0}, \
		"Metabolism":{"aldehyde dehydrogenase 1A1":0,"cytochrome":0,"quinone oxidoreductase":0, "chymotrypsin-like elastase family member 2A":0,\
		"carboxypeptidase Q":0,"glutamine amidotransferase class 1 domain containing 1":0},\
		"uncharacterized/function unknown":{"uncharacterized":0,"protein FRG2":0,"NUT family member":0},
		"others":{"asporin":0,"collagen":0,"coagulation":0, "pancreatic polypeptide":0}}
		category_index = count_filtered_by_category_index(filtered_df,big_category_dict)
		plot_dict = {}

		for big_term, sub_terms_dict in category_index.items():
			plot_dict[big_term] = 0
			big_term_count = 0
			for sub_term, count in sub_terms_dict.items():
				big_term_count += count
			plot_dict[big_term] = big_term_count
		single_series_plot_df = pd.DataFrame.from_dict(plot_dict, orient = "index",columns = ["count"])
		single_series_plot_df = single_series_plot_df.rename_axis('category').reset_index()
		series_df_list.append(single_series_plot_df)
	#plot_df = pd.merge(series_df_list[0], series_df_list[1],how = 'inner',on ='category')
	#print(plot_df)

	fig,ax = plt.subplots(1,1,figsize = (16,4))
	plot_df = pd.merge(series_df_list[0], series_df_list[1],how = 'inner',on ='category')

	#plot_df.rename(columns = {'count_x':'Woolly mammoths vs Asian elephants','count_y':"Woolly mammoths vs African elephants"})
	ax = plot_df.plot.bar()
	ax.set_xticks(np.arange(len(plot_df['category'].index)))
	ax.set_xticks(np.arange(len(plot_df['category'].index)))
	ax.set_xticklabels(plot_df['category'], rotation = 0)
	ax.legend(['Woolly mammoths vs Asian elephants',"Woolly mammoths vs African elephants"])
	plt.show()


#for result section
comparisons = ["M.primigenius.E467.L163.L164.P005.M6_vs_E.maximus.T13.T16.T18",
"M.primigenius.E467.L163.L164.P005.M6_vs_L.africana.T11",
"M.primigenius.E467.L163.L164.P005.M6_vs_M.americanum.ERR2260503",
"M.primigenius.E467.L163.L164.P005.M6_vs_P.antiquus.ERR2260504"]
targets = ["Genes", "Cis"]
cis_big_category_index = {"olfactory receptor":{"olfactory receptor":0,"trace amine-associated receptor 8":0,"gustatory receptor clone PTE03":0},\
"immune system related":{"T cell receptor":0,"natural killer cells antigen":0,"immunoglobulin":0,"histocompatibility antigen":0,\
"small proline-rich protein":0,"defensin":0,"interferon":0,"guanylate-binding protein":0, "E3 ubiquitin-protein ligase":0, "schlafen-like protein":0,\
"NKG2-A/NKG2-B type II integral membrane protein":0,"copper transport protein ATOX1":0,"tripartite motif":0,\
"NLR family pyrin domain containing":0,"toll-like receptor":0,"protein S100-A7":0,"TNF receptor superfamily member 17":0,\
"late cornified envelope-like proline-rich protein":0,"basic salivary proline-rich protein 1":0,"RHO family interacting cell polarization regulator 2":0,\
"T-cell receptor-associated transmembrane adapter":0},\
"DNA and RNA related":{'transfer RNA':0,"small nucleolar":0,"spliceosomal":0,"ribonuclease":0,"transcription initiation factor TFIID subunit 4":0,"small nuclear RNA":0}, \
"hair and body growth, metabolism":{"collagen":0,"asporin":0,"keratin-associated protein":0,"chymotrypsin-like elastase family member 2A":0,"endothelin receptor type A":0,"butyrophilin subfamily 1 member A1":0,\
"aldehyde dehydrogenase 1A1":0,"cytochrome P450":0,"adhesion G protein-coupled receptor L4":0,"quinone oxidoreductase":0, \
"carboxypeptidase Q":0,"glutamine amidotransferase class 1 domain containing 1":0,"statherin":0,"crystallin beta A1":0," leucine-rich repeat-containing protein 37A2":0},
"uncharacterized/function unknown":{"uncharacterized":0,"protein FRG2":0,"NUT family member":0},
"movement":{"myosin":0}, \
"transport":{"solute carrier":0,"aquaporin":0,"UDP-glucuronosyltransferase":0,"heparan sulfate 2-O-sulfotransferase":0},\
"cell development and cell fate":{"yippe":0,"ransforming protein RhoA-like":0,"NUMB like endocytic adaptor protein":0,"serine/threonine-protein kinase MARK2-like":0,"dual serine/threonine and tyrosine protein kinase":0},\
"others":{}}

all_samples = {}
for target in targets:
	all_samples[target] = []
	for sample in comparisons:
		handle = "./id_result/Filtered_Candidates_"+ target +'_'+sample+"_NCPG25"
		#print(sample, target)

		df = read_filtered_candidates(handle)
		#print(df.iloc[-1])
		#print(count_filtered_by_category_index(df,cis_big_category_index)['immune system related'])
		#search by single phrase
		#print(search_info)
		#print(search_info_occurance(df,search_info,mode = "short"))
		#print(len(df.index))
		all_samples[target].append(df)

#plot one single gene or cis regions 

def plot_target_by_id(geneID,species2,region):
	if species2 == "ASE":
		target_df = all_samples[region][0]
	elif species2 == "AFE":
		target_df = all_samples[region][1]

	target_dict = target_df[target_df['GeneID'] == geneID].copy() #.to_dict(orient = 'records')[0]
	print(target_dict)
	target_dict = target_dict.to_dict(orient = 'records')[0]

	plot_single_target(target_dict,species2)

geneID = "126068250"
species2 = "ASE"
region = "Genes"
plot_target_by_id(geneID,species2,region)


#plot categories bar plot
"""
gene_filtered_df_list = all_samples['Genes'][:2]

gene_big_category_dict = {"olfactory receptor":{"olfactory receptor":0,"trace amine-associated receptor 8":0,"gustatory receptor clone PTE03":0},\
"immune system related":{"butyrophilin subfamily 1 member A1":0,"MHC class I":0,"protein FAM90A27P":0,"T cell receptor":0,"natural killer cells antigen":0,"immunoglobulin":0,"histocompatibility antigen":0,\
"small proline-rich protein":0,"defensin":0,"interferon":0,"guanylate-binding protein":0, "E3 ubiquitin-protein ligase":0, "schlafen-like protein":0,\
"NKG2-A/NKG2-B type II integral membrane protein":0,"copper transport protein ATOX1":0,"tripartite motif":0,\
"NLR family pyrin domain containing":0,"toll-like receptor":0,"protein S100-A7":0,"TNF receptor superfamily member 17":0,\
"late cornified envelope-like proline-rich protein":0,"basic salivary proline-rich protein 1":0,"RHO family interacting cell polarization regulator 2":0,\
"T-cell receptor-associated transmembrane adapter":0,"interleukin":0},\
"non coding RNA":{'transfer RNA':0,"small nucleolar":0,"spliceosomal":0,"ribonuclease":0,"transcription initiation factor TFIID subunit 4":0,"small nuclear RNA":0}, \
"Metabolism":{"aldehyde dehydrogenase 1A1":0,"cytochrome":0,"quinone oxidoreductase":0, "chymotrypsin-like elastase family member 2A":0,\
"carboxypeptidase Q":0,"glutamine amidotransferase class 1 domain containing 1":0},\
"uncharacterized/function unknown":{"uncharacterized":0,"protein FRG2":0,"NUT family member":0},
"others":{"asporin":0,"collagen":0,"coagulation":0, "pancreatic polypeptide":0}}

plot_grouped_bar(gene_filtered_df_list,gene_big_category_dict)
"""
#polish the gene big category list
"""
ASE_gene_big_category_dict = polish_category(all_samples['Genes'][0],gene_big_category_dict, suffix = "ASE_genes")
AFE_gene_big_category_dict = polish_category(all_samples['Genes'][1],gene_big_category_dict, suffix = "AFE_genes")

count_filtered_by_category_index(all_samples['Genes'][0],ASE_gene_big_category_dict)
count_filtered_by_category_index(all_samples['Genes'][1],AFE_gene_big_category_dict)
"""

#search by single key word
"""
#immune related 
#search_word = "histocompatibility antigen"
#search_word = "T cell receptor"
#search_word = "natural killer"
#search_word = "NKG"
#search_word = "interferon"
#search_word = "guanylate-binding protein"
#search_word = "E3 ubiquitin-protein ligase TRIM22"
#search_word = "TNF receptor superfamily member 17"
#search_word = "defensin"
#search_word = "pancreatic polypeptide"
#search_word = "tripartite motif-containing protein"
#search_word = "chymotrypsin-like elastase family member 2A"
#search_word = "coagulation factor X"
#search_word = "transfer RNA"
#search_word = "immunoglobulin"
#search_word = "butyrophilin subfamily 1 member A1-like"
#search_word = "CD4"
#search_word = "cytochrome P450"

#search_word = "keratin-associated protein" #keep!
#search_word = "collagen"
#search_word = "asporin"

#search_word = "cytochrome"
#search_word = "aldehyde dehydrogenase 1A1" #meh 
#search_word = "adhesion G protein-coupled receptor L4" #meh angiogenesis, only one hit  
#search_word = "quinone oxidoreductase" #keep in cat, but no extensive known function
#search_word = "carboxypeptidase Q" #keep in cat, but no extensive known function
#search_word = "crystallin beta A1" #move to "others", eye lens related, may 

#search_word = "solute carrier" #may or may not include

#search_word = "olfactory receptor"

#search_word = "transfer RNA"
#search_word = "small nucleolar"
#search_word = "spliceosomal"
#search_word = "ribonuclease" #protein!!! may worth to bring it to "others section"
#search_word = "transcription initiation factor TFIID subunit 4" #this is protein, need to exclude 
#search_word = "small nuclear"
#search_word = "RNA"
#search_word = "cytochrome c"
#search_word = "interleukin"
#search_word = "CD94"
search_word = "guanylate-binding protein"
gene_compare_df_list = []
for i in all_samples['Genes'][:2]:
	count_g, info_g = search_info_occurance(i,search_word,mode = "extensive")
	print(len(i.index))

	info_g['Diff_in_f_mean'] = info_g['Diff_in_f_mean'].astype(float)
	
	print(count_g)
	print(info_g[['description','GeneID','Diff_in_f_mean']])
	gene_compare_df_list.append(info_g)
gene_check_df = pd.merge(gene_compare_df_list[0],gene_compare_df_list[1],how = 'inner',on = ['GeneID'])
print("shared:")
print(len(gene_check_df.index))
#print(gene_check_df[['description_x','GeneID']])

print("Below are cis regulatory results")

cis_compare_df_list = []
for j in all_samples['Cis'][:2]:
	print(len(j.index))
	#revised_gene_big_category = 

	count_p, info_p = search_info_occurance(j,search_word,mode = "extensive")
	info_p['Diff_in_f_mean'] = info_p['Diff_in_f_mean'].astype(float)
	print(count_p)
	print(info_p[['description','GeneID','Diff_in_f_mean']])
	cis_compare_df_list.append(info_p)
cis_check_df = pd.merge(cis_compare_df_list[0],cis_compare_df_list[1],how = 'inner', on = ['GeneID'])
print("shared:")
print(len(cis_check_df))
"""

"""
ASE_cis_filtered_handle = "./id_result/Filtered_Candidates_Cis_M.primigenius.E467.L163.L164.P005.M6_vs_E.maximus.T13.T16.T18_NCPG25"
ASE_gene_filtered_handle = "./id_result/Filtered_Candidates_Genes_M.primigenius.E467.L163.L164.P005.M6_vs_E.maximus.T13.T16.T18_NCPG25"

AFE_cis_filtered_handle = './id_result/Filtered_Candidates_Cis_M.primigenius.E467.L163.L164.P005.M6_vs_L.africana.T11_NCPG25'

ASE_cis_filtered_df = read_filtered_candidates(ASE_cis_filtered_handle)
ASE_gene_filtered_df = read_filtered_candidates(ASE_gene_filtered_handle)

AFE_cis_filtered_df = read_filtered_candidates(AFE_cis_filtered_handle)

ASE_cis_size = len(ASE_cis_filtered_df.index)
AFE_cis_size = len(AFE_cis_filtered_df.index)
"""
#plot one single gene or cis regions 
"""
geneID = "126064590"
ASE_cis_df = all_samples['Cis'][0]
ASE_gene_df = all_samples['Genes'][0]

promoter_dict = ASE_cis_df[ASE_cis_df['GeneID'] == geneID].copy() #.to_dict(orient = 'records')[0]
#print(promoter_dict)
promoter_dict = promoter_dict.to_dict(orient = 'records')[0]

gene_dict = ASE_gene_df[ASE_gene_df['GeneID'] == geneID].copy()
print(gene_dict)
gene_dict = gene_dict.to_dict(orient = 'records')[0]

print(gene_dict)
plot_single_target(gene_dict,"ASE")
plot_single_target(promoter_dict,"ASE")
#cis_dict = all_samples['Cis'][0].to_dict(orient = 'records')[13]
#gene_dict = all_samples['Genes'][0].to_dict(orient = 'records')[1]
#gene_dict = ASE_gene_filtered_df[ASE_gene_filtered_df['chrom'] == 'NC_064828.1'].to_dict(orient = 'records')[19]
#print(gene_dict)
#plot_single_target(gene_dict,"ASE")
#print(cis_dict)
#plot_single_target(cis_dict,"ASE")


#merged_cis_filered_df = pd.concat([ASE_cis_filtered_df,AFE_cis_filtered_df]).drop_duplicates().reset_index(drop=True)
#intersec_df = pd.merge(ASE_cis_filtered_df,AFE_cis_filtered_df,on = 'GeneID')
#shared_size = len(merged_df.index)
"""
#!extremely important 

cis_AFE_big_category_index = {"olfactory receptor":{"olfactory receptor":0,"trace amine-associated receptor 8":0,"gustatory receptor clone PTE03":0},\
"immune system related":{"T cell receptor":0,"natural killer cells antigen":0,"immunoglobulin":0,"histocompatibility antigen":0,\
"small proline-rich protein":0,"defensin":0,"interferon":0,"guanylate-binding protein":0, "E3 ubiquitin-protein ligase":0, "schlafen-like protein":0,\
"NKG2-A/NKG2-B type II integral membrane protein":0,"copper transport protein ATOX1":0,"tripartite motif":0,\
"NLR family pyrin domain containing":0,"toll-like receptor":0,"protein S100-A7":0,"TNF receptor superfamily member 17":0,\
"late cornified envelope-like proline-rich protein":0,"basic salivary proline-rich protein 1":0,"RHO family interacting cell polarization regulator 2":0,\
"T-cell receptor-associated transmembrane adapter":0},\
"DNA and RNA related":{'transfer RNA':0,"small nucleolar":0,"spliceosomal":0,"ribonuclease":0,"transcription initiation factor TFIID subunit 4":0,"small nuclear RNA":0}, \
"hair and body growth, metabolism":{"collagen":0,"asporin":0,"keratin-associated protein":0,"chymotrypsin-like elastase family member 2A":0,"endothelin receptor type A":0,"butyrophilin subfamily 1 member A1":0,\
"aldehyde dehydrogenase 1A1":0,"cytochrome P450":0,"adhesion G protein-coupled receptor L4":0,"quinone oxidoreductase":0, \
"carboxypeptidase Q":0,"glutamine amidotransferase class 1 domain containing 1":0,"statherin":0,"crystallin beta A1":0," leucine-rich repeat-containing protein 37A2":0},
"uncharacterized/function unknown":{"uncharacterized":0,"protein FRG2":0,"NUT family member":0},
"movement":{"myosin":0}, \
"transport":{"solute carrier":0,"aquaporin":0,"UDP-glucuronosyltransferase":0,"heparan sulfate 2-O-sulfotransferase":0},\
"cell development and cell fate":{"yippe":0,"ransforming protein RhoA-like":0,"NUMB like endocytic adaptor protein":0,"serine/threonine-protein kinase MARK2-like":0,"dual serine/threonine and tyrosine protein kinase":0},\
"others":{}}

cis_AFE_big_category_index = {"olfactory receptor":{"olfactory receptor":0,"trace amine-associated receptor 8":0,"gustatory receptor clone PTE03":0},\
"immune system related":{"T cell receptor":0,"natural killer cells antigen":0,"immunoglobulin":0,"histocompatibility antigen":0,\
"small proline-rich protein":0,"defensin":0,"interferon":0,"guanylate-binding protein":0, "E3 ubiquitin-protein ligase":0, "schlafen-like protein":0,\
"NKG2-A/NKG2-B type II integral membrane protein":0,"copper transport protein ATOX1":0,"tripartite motif":0,\
"NLR family pyrin domain containing":0,"toll-like receptor":0,"protein S100-A7":0,"TNF receptor superfamily member 17":0,\
"late cornified envelope-like proline-rich protein":0,"basic salivary proline-rich protein 1":0,"RHO family interacting cell polarization regulator 2":0,\
"T-cell receptor-associated transmembrane adapter":0},\
"DNA and RNA related":{'transfer RNA':0,"small nucleolar":0,"spliceosomal":0,"ribonuclease":0,"transcription initiation factor TFIID subunit 4":0,"small nuclear RNA":0}, \
"hair and body growth, metabolism":{"collagen":0,"asporin":0,"keratin-associated protein":0,"chymotrypsin-like elastase family member 2A":0,"endothelin receptor type A":0,"butyrophilin subfamily 1 member A1":0,\
"aldehyde dehydrogenase 1A1":0,"cytochrome P450":0,"adhesion G protein-coupled receptor L4":0,"quinone oxidoreductase":0, \
"carboxypeptidase Q":0,"glutamine amidotransferase class 1 domain containing 1":0,"statherin":0,"crystallin beta A1":0," leucine-rich repeat-containing protein 37A2":0},
"uncharacterized/function unknown":{"uncharacterized":0,"protein FRG2":0,"NUT family member":0},
"movement":{"myosin":0}, \
"transport":{"solute carrier":0,"aquaporin":0,"UDP-glucuronosyltransferase":0,"heparan sulfate 2-O-sulfotransferase":0},\
"cell development and cell fate":{"yippe":0,"ransforming protein RhoA-like":0,"NUMB like endocytic adaptor protein":0,"serine/threonine-protein kinase MARK2-like":0,"dual serine/threonine and tyrosine protein kinase":0},\
"others":{}}


#polished_AFE_category_index = polish_category(AFE_cis_filtered_df,cis_AFE_big_category_index,mode = 'cis_AFE')
#plot_lollipop(AFE_cis_filtered_df,polished_AFE_category_index)











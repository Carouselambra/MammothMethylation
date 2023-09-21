#!/usr/bin/python3

#this script is to plot single gene/regulatory region against their background

from __future__ import division
import pandas as pd
import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure
from scipy import stats
from matplotlib.patches import Rectangle
from matplotlib.ticker import NullLocator


def read_annotation(gff_handle):
	"""
	this function is used to read gff file Asian elephant genome
	input: 
		gff_handle: file path to reference annotation file 
	output:
		full_annotation_df
	"""
	raw_annotation_df = pd.read_table(gff_handle, comment = "#", skiprows = [0], header = None)
	#print(raw_annotation_df.columns)
	#print(raw_annotation_df.iloc[[1],:])
	full_annotation_df = raw_annotation_df.iloc[:,[0,2,3,4,8]]
	new_col_names = {0:"chromosome",2:"type",3:"start",4:"end",8:"Target_Info"}
	full_annotation_df = full_annotation_df.rename(columns = new_col_names)

	return full_annotation_df

def convert_chrom_name(chrom_number):
	chrom_list = ["NC_064819.1","NC_064820.1","NC_064821.1","NC_064822.1","NC_064823.1","NC_064824.1", \
	"NC_064825.1","NC_064826.1","NC_064827.1","NC_064828.1","NC_064829.1","NC_064830.1","NC_064831.1", \
	"NC_064832.1","NC_064833.1","NC_064834.1","NC_064835.1","NC_064836.1","NC_064837.1","NC_064838.1", \
	"NC_064839.1","NC_064840.1","NC_064841.1","NC_064842.1","NC_064843.1","NC_064844.1","NC_064845.1", \
	"NC_064846.1","NC_064847.1"]
	out_put = ""
	for i in range(len(chrom_list)-1):
		if chrom_list[i] == chrom_number:
			#shouldn't include sex chromosome in the candidate gene list 
			if chrom_number == "NC_064846.1":
				out_put = "chr{number} ({NC_num})".format(number = "X", NC_num = chrom_number)
			elif chrom_number == "NC_064847.1":
				out_put = "chr{number} ({NC_num})".format(number = "Y", NC_num = chrom_number)
			else:
				out_put = "chr{number} ({NC_num})".format(number = str(i+1), NC_num = chrom_number) 
	return out_put


def read_filtered_candidates(filtered_handle):
	filtered_df = pd.read_csv(filtered_handle,sep = '\t',usecols=["Diff_in_f_mean_abs","Diff_in_f_mean","start","end","chrom","Target_Info"])
	filtered_df = filtered_df[filtered_df['Target_Info'].str.contains('description=')].copy()
	filtered_df['description'] = filtered_df['Target_Info'].apply(lambda x: str(x.partition('description=')[-1].partition(';')[0]))
	filtered_df['GeneID'] = filtered_df['Target_Info'].apply(lambda x: str(x.partition('GeneID:')[-1].partition(';')[0]))
	return filtered_df

def plot_single_gene_ver2(target_dict,species2):
	"""this plot should take one gene / regulatory region coordinates 
	then plot their background 
	input:
		target_dict: contains keys "chrom", "gene_start", "gene_end","GeneID","IRR_start","IRR_end"
	output:
	"""
	chrom = target_dict['chrom']
	id_text = target_dict['GeneID']

	IRR_start = int(target_dict['IRR_start'])
	IRR_end = int(target_dict['IRR_end'])
	gene_start = int(target_dict['gene_start'])
	gene_end = int(target_dict['gene_end'])

	

	#determine the nature of the region
	#if target_end - target_start <= 5000:
	#	plot_mode = "small_gene"
		#plot_size = 100
	#else: 
	#	plot_mode = "large_gene"

	#find the first and the last window of the highlighted regions
	#target_first_window_start = target_start // 1000 * 1000 +1
	#target_last_window_start = target_end // 1000 * 1000+1

	#target_first_window = "{start}-{end}".format(start = str(target_first_window_start), end = str(target_first_window_start+1000))
	#target_last_window = "{start}-{end}".format(start = str(target_last_window_start), end = str(target_last_window_start+1000))

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
	
	#determine how much area will be covered in one plot 
	if plot_mode == "small_gene":
		subplot1_index_start = target_index_start_1 - 25
		subplot1_index_end = target_index_end_1 + 70
	
		subplot2_index_start = target_index_start_2 - 25
		subplot2_index_end = target_index_end_2 + 70

	elif plot_mode == "large_gene":
		subplot1_index_start = target_index_start_1 - 50
		subplot1_index_end = target_index_end_1 + 50
	
		subplot2_index_start = target_index_start_2 - 50
		subplot2_index_end = target_index_end_2 + 50
	print(plot_mode)

	#select the background region 
	plot_df_1 = nF_df_1.iloc[subplot1_index_start:subplot1_index_end].copy().reset_index()
	plot_df_2 = nF_df_2.iloc[subplot2_index_start:subplot2_index_end].copy().reset_index()
	
	highlight1_start = plot_df_1[plot_df_1['window_location'] == target_first_window].index.values.astype(int)[0]
	highlight1_end = plot_df_1[plot_df_1['window_location'] == target_last_window].index.values.astype(int)[0]
	#print(highlight1_start,highlight1_end)
	highlight2_start = plot_df_2[plot_df_2['window_location'] == target_first_window].index.values.astype(int)[0]
	highlight2_end = plot_df_2[plot_df_2['window_location'] == target_last_window].index.values.astype(int)[0]

	plot_df_1['window_start'] = plot_df_1['window_location'].str.partition('-')[0]
	plot_df_2['window_start'] = plot_df_2['window_location'].str.partition('-')[0]

	#generate the figure 
	fig, (ax1,ax2) = plt.subplots(2,1, figsize = (12,4),sharex = True, sharey = True)
	fig.canvas.draw()

	#viridis color palette three categories
	clrs_1 = ['#fde725' if x < 0.8  else ('#440154' if x > 1.2 else '#21918c') for x in plot_df_1['normalised_f_fold']]
	ax1.bar(plot_df_1.index,plot_df_1['normalised_f_fold'],label = 'M.primigenius',color = clrs_1)
	ax1.legend()

	clrs_2 = ['#fde725' if x < 0.8  else ('#440154' if x > 1.2 else '#21918c') for x in plot_df_2['normalised_f_fold']]
	ax2.bar(plot_df_2.index,plot_df_2['normalised_f_fold'],color = clrs_2,label = label_2)
	ax2.legend()

	#draw the gene region annotation 
	

	#plot target region at subplot 2 
	#ax2.axvspan(highlight2_start,highlight2_end,alpha = 0.5)

	#plot the normalised average line 
	ax2.axhline(y = 1, linestyle =  '--', color = "#000000")

	# labeled the x ticks (chromosome choordination)
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
	#ax2.set_xlabel(chrom + ':' + '\n' + 'coordinates')
	#ax2.xaxis.set_label_coords(-0.07, -0.06)

	plt.show()

def get_full_target_by_id(geneID,species2):
	#all_samples = {}
	comparisons = ["M.primigenius.E467.L163.L164.P005.M6_vs_E.maximus.T13.T16.T18",
	"M.primigenius.E467.L163.L164.P005.M6_vs_L.africana.T11",
	"M.primigenius.E467.L163.L164.P005.M6_vs_M.americanum.ERR2260503",
	"M.primigenius.E467.L163.L164.P005.M6_vs_P.antiquus.ERR2260504"]
	species_2  = ["ASE","AFE","MSTD","STKS"]
	for i in range(len(species_2)):
		if species2 == species_2[i]:
			#all_samples[species2] = []
			cis_handle = "./id_result/Filtered_Candidates_Cis_"+comparisons[i] +"_NCPG25"
			cis_df = read_filtered_candidates(cis_handle)
			#all_samples[species2].append(cis_df)
			gene_handle = "./id_result/Filtered_Candidates_Genes_"+comparisons[i] +"_NCPG25"
			gene_df = read_filtered_candidates(gene_handle)
			#all_samples[species2].append(gene_df)

	target_dict = {}
	target_dict['GeneID'] = geneID
	target_dict['gene_start'] = 0 
	target_dict['gene_end'] = 0
	target_dict['IRR_start'] = 0
	target_dict['IRR_end'] = 0
	target_dict['description'] = ""


	#search id through cis candidate lists, get cis start and end coordinates  
	if not cis_df[cis_df['GeneID'] == geneID].empty:
		target_dict_cis = cis_df[cis_df['GeneID'] == geneID].copy().to_dict(orient = 'records')[0]
		target_dict['IRR_start'] = target_dict_cis['start']
		target_dict['IRR_end'] = target_dict_cis['end']
		target_dict['chrom'] = target_dict_cis['chrom']
		print("in cis df")

	#serarch id through gene body candidate list, get gene body start and end coordinates 
	if not gene_df[gene_df['GeneID'] == geneID].empty:
		target_dict_gene = gene_df[gene_df['GeneID'] == geneID].copy().to_dict(orient = 'records')[0]
		target_dict['gene_start'] = target_dict_gene['start']
		target_dict['gene_end'] = target_dict_gene['end']
		target_dict['chrom'] = target_dict_gene['chrom']
		print("in gene df")

	gff_handle = "GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff"
	full_annotation_df = read_annotation(gff_handle)

	for k, v in target_dict.items():
		if v == 0:
			#cis empty, then calculate IRR by start and end
			if 'IRR' in k:
				if target_dict['gene_start'] != 0:
					target_dict['IRR_start'] = target_dict['gene_start'] - 5000
					target_dict['IRR_end'] = target_dict['gene_start']

			#gene body empty, then search this gene start and end in full reference file
			if 'gene' in k:
				target_dict_gene = full_annotation_df[full_annotation_df['Target_Info'].str.contains(geneID)].copy().to_dict(orient = 'records')[0]
				target_dict['gene_start'] = target_dict_gene['start']
				target_dict['gene_end'] = target_dict_gene['end']
		if k == "description":
			target_dict_gene = full_annotation_df[full_annotation_df['Target_Info'].str.contains(geneID)].copy().to_dict(orient = 'records')[0]
			target_dict['description'] = target_dict_gene['Target_Info'].partition('description=')[-1].partition(';')[0]
	print(target_dict)
	return target_dict

def plot_full_target(target_dict,species2,full_plot_range = 150,context_ahead_bar_count = 50):
	"""this function is to plot the normalised F rate 
	with gene IRR and body regions marked below as shape 
	(highlight first as well?)
	"""
	#full_plot_range = 150 #one plot 
	#context_ahead_bar_count = 50

	chrom = target_dict['chrom']
	chrom_label = convert_chrom_name(chrom)
	id_text = target_dict['GeneID']
	description = target_dict['description']

	IRR_start = int(target_dict['IRR_start'])
	IRR_end = int(target_dict['IRR_end'])
	gene_start = int(target_dict['gene_start'])
	gene_end = int(target_dict['gene_end'])

	if gene_end - gene_start > 100000:
		print("This gene needs to be plot separately!")
		print(target_dict)
	else:


		first_bar_start = (IRR_start - context_ahead_bar_count*1000) // 1000 * 1000 +1
		first_bar_x = "{start}-{end}".format(start = str(first_bar_start), end = str(first_bar_start+1000))

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

		subplot1_index_start = nF_df_1[nF_df_1['window_location'] == first_bar_x].index.values.astype(int)[0]
		subplot2_index_start = nF_df_2[nF_df_2['window_location'] == first_bar_x].index.values.astype(int)[0]

		if subplot1_index_start != subplot2_index_start:
			print("your index doesn't match!!!")
		else:
			#select the background region 
			plot_df_1 = nF_df_1.iloc[subplot1_index_start:subplot1_index_start+full_plot_range].copy().reset_index()
			plot_df_2 = nF_df_2.iloc[subplot2_index_start:subplot2_index_start+full_plot_range].copy().reset_index()

		#get the indexes of IRR and gene body 
		#get IRR region binned window start and end
		IRR_window_start = IRR_start // 1000 * 1000 + 1
		IRR_first_window = "{start}-{end}".format(start = str(IRR_window_start), end = str(IRR_window_start+1000))
		
		IRR_window_end = IRR_end // 1000 * 1000 + 1
		IRR_last_window = "{start}-{end}".format(start = str(IRR_window_end), end = str(IRR_window_end + 1000))

		#get the index of IRR start and end within plot region
		IRR_start_index = plot_df_1[plot_df_1['window_location'] == IRR_first_window].index.values.astype(int)[0]
		IRR_end_index = plot_df_1[plot_df_1['window_location'] == IRR_last_window].index.values.astype(int)[0]
		print(IRR_start_index, IRR_end_index)


		#get gene region binned window start and end
		gene_window_start = gene_start // 1000 * 1000 + 1
		gene_first_window = "{start}-{end}".format(start = str(gene_window_start), end = str(gene_window_start+1000))

		gene_window_end = gene_end // 1000 * 1000 + 1
		gene_last_window = "{start}-{end}".format(start = str(gene_window_end), end = str(gene_window_end + 1000))

		#get the index of gene body start and end within plot region
		gene_start_index = plot_df_1[plot_df_1['window_location'] == gene_first_window].index.values.astype(int)[0]
		gene_end_index = plot_df_1[plot_df_1['window_location'] == gene_last_window].index.values.astype(int)[0]
		print(gene_start_index, gene_end_index)

		#generate the figure 
		fig, (ax1,ax2,ax3) = plt.subplots(3,1, figsize = (12,6), gridspec_kw={'height_ratios': [6,6,1]}, sharex = True)
		#bar graph share y axis
		ax2.sharey(ax1)

		fig.canvas.draw()

		#viridis color palette three categories
		clrs_1 = ['#66c2a5' if x < 0.4  else ('#abdda4' if x < 0.8 else ('#e6f598' if x < 1.2 else ("#fdae61" if x < 1.6 else '#f46d43'))) for x in plot_df_1['normalised_f_fold']]
		ax1.bar(plot_df_1.index,plot_df_1['normalised_f_fold'],label = 'M.primigenius',color = clrs_1)
		ax1.set_ylim(0,7)
		ax1.axhline(y = 1, xmin = 0.045,linestyle =  '--', color = "#000000")
		ax1.text(-0.1,0.5,'M.primigenius',transform = ax1.transAxes)
		ax1.text(-0.1,0.35,'F fold change',transform = ax1.transAxes)

		#add a scale line (5000bp)
		ax1.axhline(y = 6.4, xmin = 0.1, xmax = 0.1 + 4/full_plot_range, color = "#000000")
		ax1.text(0.09, 0.95, '5000 bp',transform = ax1.transAxes)

		clrs_2 = ['#66c2a5' if x < 0.4  else ('#abdda4' if x < 0.8 else ('#e6f598' if x < 1.2 else ("#fdae61" if x < 1.6 else '#f46d43'))) for x in plot_df_2['normalised_f_fold']]
		ax2.bar(plot_df_2.index,plot_df_2['normalised_f_fold'],color = clrs_2,label = label_2)
		ax2.axhline(y = 1, xmin = 0.045,linestyle =  '--', color = "#000000")
		ax2.text(-0.1,0.5,label_2, transform = ax2.transAxes)
		ax2.text(-0.1,0.35,'F fold change',transform = ax2.transAxes)

		#set ticks and label
		x_ticks_list = [0,len(plot_df_1.index)]
		ax2.set_xticks(x_ticks_list)

		last_bar_end = plot_df_2['window_location'].iloc[-1].partition('-')[-1]
		ax3.set_xticklabels([first_bar_start,last_bar_end],ha = 'left')

		ax3.axhline(y = 0.25, xmin = 0.045,zorder=1,color = "#000000")
		ax3.text(0,1,chrom_label)

		#IRR_rect = mpl.patches.Rectangle((IRR_start_index, 0), width=IRR_end_index - IRR_start_index, height=0.5, color="#2DD881", clip_on=False, fill = True, linewidth = 1.5, zorder=5)
		IRR_arrow = mpl.patches.Arrow(IRR_start_index,0.25,IRR_end_index - IRR_start_index,0,width = 1.5,color="#2DD881", clip_on=False, fill = True)
		ax3.add_patch(IRR_arrow)
		gene_rect = mpl.patches.Rectangle((gene_start_index, 0), width=gene_end_index - gene_start_index, height=0.5, color="#a7a7a7", clip_on=False, fill = True,zorder=5)
		ax3.add_patch(gene_rect)

		ax3.text(gene_start_index, -1, '(Entrez ID:{geneID})'.format(geneID = geneID),zorder = 5)
		ax3.text(gene_start_index, -0.5,description,zorder = 5)
		ax3.tick_params(left = False, bottom = False,labelleft = False)


		ax1.spines['bottom'].set_position('zero')
		ax2.spines['bottom'].set_position('zero')
		ax1.spines['left'].set_position(('data',-0.5))
		ax2.spines['left'].set_position(('data',-0.5))


		ax1.spines['top'].set_visible(False)
		ax2.spines['top'].set_visible(False)
		ax3.spines['top'].set_visible(False)
		ax1.spines['right'].set_visible(False)
		ax2.spines['right'].set_visible(False)
		ax3.spines['right'].set_visible(False)
		ax3.spines['left'].set_visible(False)
		ax3.spines['bottom'].set_visible(False)

		if "/" in description:
			final_description = description.replace("/"," or ")
		else:
			final_description = description

		plt.savefig('/Users/arcquinn/Documents/MammothMethylation/SCRIPTS/manuscript_fig/gene_plots/gene_plot_{species2}_{id}_{name}.pdf'.format(species2 = species2,id = geneID, name = final_description))
		#plt.show()
		plt.close()


def plot_full_target_unfilter(target_dict,species2,full_plot_range = 150,context_ahead_bar_count = 50):
	"""this function is to plot the normalised F rate 
	with gene IRR and body regions marked below as shape 
	(highlight first as well?)
	"""
	#full_plot_range = 150 #one plot 
	#context_ahead_bar_count = 50

	chrom = target_dict['chrom']
	chrom_label = convert_chrom_name(chrom)
	id_text = target_dict['GeneID']
	description = target_dict['description']

	IRR_start = int(target_dict['IRR_start'])
	IRR_end = int(target_dict['IRR_end'])
	gene_start = int(target_dict['gene_start'])
	gene_end = int(target_dict['gene_end'])


	first_bar_start = (IRR_start - context_ahead_bar_count*1000) // 1000 * 1000 +1
	first_bar_x = "{start}-{end}".format(start = str(first_bar_start), end = str(first_bar_start+1000))

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

	subplot1_index_start = nF_df_1[nF_df_1['window_location'] == first_bar_x].index.values.astype(int)[0]
	subplot2_index_start = nF_df_2[nF_df_2['window_location'] == first_bar_x].index.values.astype(int)[0]

	if subplot1_index_start != subplot2_index_start:
		print("your index doesn't match!!!")
	else:
		#select the background region 
		plot_df_1 = nF_df_1.iloc[subplot1_index_start:subplot1_index_start+full_plot_range].copy().reset_index()
		plot_df_2 = nF_df_2.iloc[subplot2_index_start:subplot2_index_start+full_plot_range].copy().reset_index()

	#get the indexes of IRR and gene body 
	#get IRR region binned window start and end
	IRR_window_start = IRR_start // 1000 * 1000 + 1
	IRR_first_window = "{start}-{end}".format(start = str(IRR_window_start), end = str(IRR_window_start+1000))
	
	IRR_window_end = IRR_end // 1000 * 1000 + 1
	IRR_last_window = "{start}-{end}".format(start = str(IRR_window_end), end = str(IRR_window_end + 1000))

	#get the index of IRR start and end within plot region
	IRR_start_index = plot_df_1[plot_df_1['window_location'] == IRR_first_window].index.values.astype(int)[0]
	IRR_end_index = plot_df_1[plot_df_1['window_location'] == IRR_last_window].index.values.astype(int)[0]
	print(IRR_start_index, IRR_end_index)


	#get gene region binned window start and end
	gene_window_start = gene_start // 1000 * 1000 + 1
	gene_first_window = "{start}-{end}".format(start = str(gene_window_start), end = str(gene_window_start+1000))

	gene_window_end = gene_end // 1000 * 1000 + 1
	gene_last_window = "{start}-{end}".format(start = str(gene_window_end), end = str(gene_window_end + 1000))

	#get the index of gene body start and end within plot region
	gene_start_index = plot_df_1[plot_df_1['window_location'] == gene_first_window].index.values.astype(int)[0]
	gene_end_index = plot_df_1[plot_df_1['window_location'] == gene_last_window].index.values.astype(int)[0]
	print(gene_start_index, gene_end_index)

	#generate the figure 
	fig, (ax1,ax2,ax3) = plt.subplots(3,1, figsize = (12,6), gridspec_kw={'height_ratios': [6,6,1]}, sharex = True)
	#bar graph share y axis
	ax2.sharey(ax1)

	fig.canvas.draw()

	#viridis color palette three categories
	clrs_1 = ['#66c2a5' if x < 0.4  else ('#abdda4' if x < 0.8 else ('#e6f598' if x < 1.2 else ("#fdae61" if x < 1.6 else '#f46d43'))) for x in plot_df_1['normalised_f_fold']]
	ax1.bar(plot_df_1.index,plot_df_1['normalised_f_fold'],label = 'M.primigenius',color = clrs_1)
	ax1.set_ylim(0,7)
	ax1.axhline(y = 1, xmin = 0.045,linestyle =  '--', color = "#000000")
	ax1.text(-0.1,0.5,'M.primigenius',transform = ax1.transAxes)
	ax1.text(-0.1,0.35,'F fold change',transform = ax1.transAxes)

	#add a scale line (5000bp)
	ax1.axhline(y = 6.4, xmin = 0.1, xmax = 0.1 + 4/full_plot_range, color = "#000000")
	ax1.text(0.09, 0.95, '5000 bp',transform = ax1.transAxes)

	clrs_2 = ['#66c2a5' if x < 0.4  else ('#abdda4' if x < 0.8 else ('#e6f598' if x < 1.2 else ("#fdae61" if x < 1.6 else '#f46d43'))) for x in plot_df_2['normalised_f_fold']]
	ax2.bar(plot_df_2.index,plot_df_2['normalised_f_fold'],color = clrs_2,label = label_2)
	ax2.axhline(y = 1, xmin = 0.045,linestyle =  '--', color = "#000000")
	ax2.text(-0.1,0.5,label_2, transform = ax2.transAxes)
	ax2.text(-0.1,0.35,'F fold change',transform = ax2.transAxes)

	#set ticks and label
	x_ticks_list = [0,len(plot_df_1.index)]
	ax2.set_xticks(x_ticks_list)

	last_bar_end = plot_df_2['window_location'].iloc[-1].partition('-')[-1]
	ax3.set_xticklabels([first_bar_start,last_bar_end],ha = 'left')

	ax3.axhline(y = 0.25, xmin = 0.045,zorder=1,color = "#000000")
	ax3.text(0,1,chrom_label)

	#IRR_rect = mpl.patches.Rectangle((IRR_start_index, 0), width=IRR_end_index - IRR_start_index, height=0.5, color="#2DD881", clip_on=False, fill = True, linewidth = 1.5, zorder=5)
	IRR_arrow = mpl.patches.Arrow(IRR_start_index,0.25,IRR_end_index - IRR_start_index,0,width = 1.5,color="#2DD881", clip_on=False, fill = True)
	ax3.add_patch(IRR_arrow)
	gene_rect = mpl.patches.Rectangle((gene_start_index, 0), width=gene_end_index - gene_start_index, height=0.5, color="#a7a7a7", clip_on=False, fill = True,zorder=5)
	ax3.add_patch(gene_rect)

	ax3.text(gene_start_index, -1, '(Entrez ID:{geneID})'.format(geneID = geneID),zorder = 5)
	ax3.text(gene_start_index, -0.5,description,zorder = 5)
	ax3.tick_params(left = False, bottom = False,labelleft = False)


	ax1.spines['bottom'].set_position('zero')
	ax2.spines['bottom'].set_position('zero')
	ax1.spines['left'].set_position(('data',-0.5))
	ax2.spines['left'].set_position(('data',-0.5))


	ax1.spines['top'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax1.spines['right'].set_visible(False)
	ax2.spines['right'].set_visible(False)
	ax3.spines['right'].set_visible(False)
	ax3.spines['left'].set_visible(False)
	ax3.spines['bottom'].set_visible(False)

	if "/" in description:
		final_description = description.replace("/"," or ")
	else:
		final_description = description

	plt.savefig('/Users/arcquinn/Documents/MammothMethylation/SCRIPTS/manuscript_fig/gene_plots/gene_plot_{species2}_{id}_{name}.pdf'.format(species2 = species2,id = geneID, name = final_description))
	#plt.show()
	plt.close()

#single gene plot

geneID = "126081159"
species2 = "ASE"
target_dict = get_full_target_by_id(geneID,species2)
plot_full_target_unfilter(target_dict,species2,full_plot_range = 500,context_ahead_bar_count = 100)


#loop: plot the top 50 for mm vs. ASE, and mm vs. AFE
"""
top_gene_num = 50
species2_list = ["ASE"]

for species2 in species2_list:
	if species2 == "ASE":
		filtered_df_path_cis = "/Users/arcquinn/Documents/MammothMethylation/SCRIPTS/id_result/Filtered_Candidates_Cis_M.primigenius.E467.L163.L164.P005.M6_vs_E.maximus.T13.T16.T18_NCPG25"
		filtered_df_path_gene = "/Users/arcquinn/Documents/MammothMethylation/SCRIPTS/id_result/Filtered_Candidates_Genes_M.primigenius.E467.L163.L164.P005.M6_vs_E.maximus.T13.T16.T18_NCPG25"
	elif species2 == "AFE":
		filtered_df_path_cis = "/Users/arcquinn/Documents/MammothMethylation/SCRIPTS/id_result/Filtered_Candidates_Cis_M.primigenius.E467.L163.L164.P005.M6_vs_L.africana.T11_NCPG25"
		filtered_df_path_gene = "/Users/arcquinn/Documents/MammothMethylation/SCRIPTS/id_result/Filtered_Candidates_Genes_M.primigenius.E467.L163.L164.P005.M6_vs_L.africana.T11_NCPG25"	
	
	filtered_df_cis = read_filtered_candidates(filtered_df_path_cis)
	filtered_df_gene = read_filtered_candidates(filtered_df_path_gene)

	cis_id_list = filtered_df_cis['GeneID'].to_list()[:top_gene_num-1]
	gene_id_list = filtered_df_gene['GeneID'].to_list()[:top_gene_num-1]

	final_geneID_list = list(set(cis_id_list + gene_id_list))
	for geneID in final_geneID_list:
		target_dict = get_full_target_by_id(geneID,species2)
		plot_full_target(target_dict,species2,full_plot_range = 150)
"""


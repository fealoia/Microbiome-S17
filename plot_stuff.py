import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import check_adjacents

'''
	Plotting code for each of the different plots
'''

path = '/Users/williampascucci/Documents/Research/MIDAS/data_out/pool_Apair_res/genes/'

species_reads,species_info = check_adjacents.check_adjacents()
#print(species_info.keys())

def get_reads_for_order(path):
	count_dict = {}
	count_list = []
	with open(path+'summary.txt','r') as fn:
		for line in fn.readlines():
			stuff_list = line.split('\t')
			if(stuff_list[0] == 'species_id'):
				continue
			#print(int(stuff_list[-1].strip()))
			count_dict[int(stuff_list[-1].strip())] = stuff_list[0]
			count_list.append(int(stuff_list[-1].strip()))
		return(count_dict,count_list)

def plot_percent_rough_reads(count_dict,count_list):
	plt.clf()
	spec_x = []
	perc_y = []
	for species in reversed(sorted(count_dict.items())):
		s_index = species[1]+'.genes.gz'
		spec_x.append(species[1])
		perc_y.append(species_info[s_index]['percent_rough']*100)

	#print(perc_y)

	fig, ax = plt.subplots()
	fig.suptitle('Pool A LowAcc Percent of rough mapped reads', fontsize=12)
	x = range(len(spec_x))
	plt.ylabel('Percentage of Rough Reads')
	plt.xticks(x,spec_x,rotation=40,ha='right')

	fmt = '%.1f%%' # Format you want the ticks, e.g. '40%'
	yticks = FormatStrFormatter(fmt)
	ax.yaxis.set_major_formatter(yticks)

	plt.tight_layout()
	plt.bar(x,perc_y)
	#plt.show()
	fig.savefig('Pool_A_LowAcc_rough_reads.png')

	return

def plot_percent_smooth_reads(count_dict,count_list):
	plt.clf()
	spec_x = []
	perc_y = []
	for species in reversed(sorted(count_dict.items())):
		s_index = species[1]+'.genes.gz'
		spec_x.append(species[1])
		perc_y.append(species_info[s_index]['percent_smooth']*100)

	#print(perc_y)

	fig, ax = plt.subplots()
	fig.suptitle('Pool A LowAcc Percent of smooth mapped reads', fontsize=12)
	x = range(len(spec_x))
	plt.ylabel('Perc of Smooth Reads')
	plt.xticks(x,spec_x,rotation=40,ha='right')

	fmt = '%.1f%%' # Format you want the ticks, e.g. '40%'
	yticks = FormatStrFormatter(fmt)
	ax.yaxis.set_major_formatter(yticks)

	plt.tight_layout()
	plt.bar(x,perc_y)
	#plt.show()
	fig.savefig('Pool_A_LowAcc_smooth_reads.png')

	return

def plot_num_families(count_dict,count_list):
	plt.clf()
	spec_x = []
	perc_y = []
	for species in reversed(sorted(count_dict.items())):
		s_index = species[1]+'.genes.gz'
		spec_x.append(species[1])
		perc_y.append(species_info[s_index]['families_count'])

	#print(perc_y)

	fig, ax = plt.subplots()
	fig.suptitle('Pool A LowAcc number of families mapped to per species', fontsize=12)
	x = range(len(spec_x))
	plt.ylabel('# Families')
	plt.xticks(x,spec_x,rotation=40,ha='right')

	plt.tight_layout()
	plt.bar(x,perc_y)
	#plt.show()
	fig.savefig('Pool_A_LowAcc_num_families.png')

	return

def plot_total_genes(count_dict,count_list):
	plt.clf()
	spec_x = []
	perc_y = []
	for species in reversed(sorted(count_dict.items())):
		s_index = species[1]+'.genes.gz'
		spec_x.append(species[1])
		perc_y.append(species_info[s_index]['genes_count'])

	#print(perc_y)

	fig, ax = plt.subplots()
	fig.suptitle('Pool A LowAcc number of total genes included per species', fontsize=12)
	x = range(len(spec_x))
	plt.ylabel('# Total Genes')
	plt.xticks(x,spec_x,rotation=40,ha='right')

	plt.bar(x,perc_y)
	plt.tight_layout()
	#plt.show()
	fig.savefig('Pool_A_LowAcc_total_genes.png')

	return

def plot_sig_mapped_genes(count_dict,count_list):
	plt.clf()
	spec_x = []
	perc_y = []
	for species in reversed(sorted(count_dict.items())):
		s_index = species[1]+'.genes.gz'
		spec_x.append(species[1])
		perc_y.append(species_info[s_index]['smooth_mappings']+species_info[s_index]['rough_mappings'])

	#print(perc_y)

	fig, ax = plt.subplots()
	fig.suptitle('Pool A LowAcc number of genes with significant mapping to per species', fontsize=12)
	x = range(len(spec_x))
	plt.ylabel('# Significant Mappings to Genes')
	plt.xticks(x,spec_x,rotation=40,ha='right')

	plt.bar(x,perc_y)
	plt.tight_layout()
	#plt.show()
	fig.savefig('Pool_A_LowAcc_#_sig_mapped_genes.png')

	return

def plot_num_mappings(count_dict,count_list):
	plt.clf()
	spec_x = []
	perc_y = []
	for species in reversed(sorted(count_dict.items())):
		s_index = species[1]+'.genes.gz'
		spec_x.append(species[1])
		perc_y.append(species_info[s_index]['smooth_mappings']+species_info[s_index]['rough_mappings']+species_info[s_index]['low_read_mappings'])

	#print(perc_y)

	fig, ax = plt.subplots()
	fig.suptitle('Pool A LowAcc number of genes with a mapping to per species', fontsize=12)
	x = range(len(spec_x))
	plt.ylabel('#Mappings to Genes')
	plt.xticks(x,spec_x,rotation=40,ha='right')

	plt.bar(x,perc_y)
	plt.tight_layout()
	#plt.show()
	fig.savefig('Pool_A_LowAcc_#_mapped_genes.png')

	return


count_dict,count_list = get_reads_for_order(path)
plot_percent_rough_reads(count_dict,count_list)
plot_percent_smooth_reads(count_dict,count_list)
plot_num_families(count_dict,count_list)
plot_total_genes(count_dict,count_list)
plot_sig_mapped_genes(count_dict,count_list)
plot_num_mappings(count_dict,count_list)



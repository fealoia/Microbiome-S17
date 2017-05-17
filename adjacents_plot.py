import sys
import gzip
import collections
import os
import matplotlib.pyplot as plt

start_path = '../data_out/pool_Apair_res/genes/output/'
species_reads = {}
species_info = {}

def build_gene_dict():
	for i in range(4):
		load_file(kleb_list[i],kleb_dict[i])
	#for bact in os.listdir(start_path)
	return	
		

def load_file(filename,bact_dict):
	with gzip.open(filename,'rt') as fn:
		for line in fn.readlines():
			#print(line.strip())
			info_list = line.split('\t')
			#print(info_list)
			if(info_list[0]=='gene_id'):
				continue
			gene_id = info_list[0].split('.')
			#print(gene_id)
			family_id = gene_id[0]+'.'+gene_id[1]
			if(family_id in bact_dict):
				bact_dict[family_id][int(gene_id[3])] = int(info_list[1])
			else:
				bact_dict[family_id] = collections.OrderedDict()
				bact_dict[family_id][int(gene_id[3])] = int(info_list[1])

def check_broken(bact_dict,species_info,spec_name):
	species_info['smooth_mappings'] = 0
	species_info['rough_mappings'] = 0
	species_info['total_mapped_reads'] = 0
	species_info['low_read_mappings'] = 0
	species_info['genes_count'] = 0
	species_info['families_count'] = 0
	y_rough = []
	y_smooth = []
	y_low_read = []
	for family in bact_dict.keys():
		species_info['families_count'] += 1
		for pos in bact_dict[family].keys():
			species_info['genes_count'] += 1
			species_info['total_mapped_reads']+=bact_dict[family][pos]
			#y.append(bact_dict[family][pos])
			if((bact_dict[family][pos]>100)and(pos-1 in bact_dict[family])and bact_dict[family][pos-1]<100):
				#print('break at '+family+'.'+str(pos))
				species_info['rough_mappings']+=1
				y_rough.append(bact_dict[family][pos])
				y_smooth.append(0)
				y_low_read.append(0)
			elif((bact_dict[family][pos]>100)and(pos+1 in bact_dict[family])and bact_dict[family][pos+1]<100):
				#print('break at '+family+'.'+str(pos))
				species_info['rough_mappings']+=1
				y_rough.append(bact_dict[family][pos])
				y_smooth.append(0)
				y_low_read.append(0)
			elif(bact_dict[family][pos]>100):
				species_info['smooth_mappings']+=1
				y_rough.append(0)
				y_smooth.append(bact_dict[family][pos])
				y_low_read.append(0)
			elif(bact_dict[family][pos] >= 1 and bact_dict[family][pos]<100):
				species_info['low_read_mappings']+=1
				y_rough.append(0)
				y_smooth.append(0)
				y_low_read.append(bact_dict[family][pos])
	x = range(species_info['rough_mappings']+species_info['smooth_mappings']+species_info['low_read_mappings'])
	plt.clf()
	plt.bar(x,y_rough,color='r',label='rough reads',edgecolor='none')
	plt.bar(x,y_smooth,color='g',label='smooth reads',edgecolor='none')
	plt.bar(x,y_low_read,color='b',label = 'low # reads',edgecolor='none')
	plt.legend()
	#plt.show()
	plt.ylabel('# of reads mapped')
	plt.xlabel('Position on pangenome')
	name = spec_name.split('.')
	plt.savefig('Pool_A95perc_all_reads_'+name[0]+'.png')

def check_unique_families(start_path):
	family_dict = {}
	for bact in os.listdir(start_path):
		with gzip.open(start_path+bact,'rt') as fn:
			for line in fn.readlines():
				info_list = line.split('\t')
				if(info_list[0]=='gene_id'):
					continue
				gene_id = info_list[0].split('.')
				#different than load_file's procedure
				family_id = gene_id[0]
				if(not(family_id) in family_dict):
					family_dict[family_id] = [bact]
				elif(not(bact in family_dict[family_id])):
					print('found!')
					family_dict[family_id].append(bact)

def load_all_bacteria_infile_intodict(start_path,species_dict,species_info):
	#print(os.listdir(start_path))
	for bact in os.listdir(start_path):
		species_info[bact] = {}
		species_dict[bact] = {}
		#print(bact)
		load_file(start_path+bact,species_dict[bact])

#check_unique_families(start_path)
#build_gene_dict()
#print(kleb_dict[1])
#check_broken(kleb_dict[1])

def check_adjacents():
	start_path = '../data_out/pool_Apair_res/genes/output/'
	species_reads = {}
	species_info = {}
	load_all_bacteria_infile_intodict(start_path,species_reads,species_info)
	for spec in species_reads:
		check_broken(species_reads[spec],species_info[spec],spec)
	return species_reads,species_info

check_adjacents()











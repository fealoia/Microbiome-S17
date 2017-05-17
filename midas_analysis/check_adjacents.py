import sys
import gzip
import collections
import os
import matplotlib.pyplot as plt

start_path = '../data_out/pool_Apair_res/main_genes/output/'
species_reads = {}
species_info = {}

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
	#y = []
	for family in bact_dict.keys():
		species_info['families_count'] += 1
		for pos in bact_dict[family].keys():
			species_info['genes_count'] += 1
			species_info['total_mapped_reads']+=bact_dict[family][pos]
			#y.append(bact_dict[family][pos])
			'''
			10-fold comparison
			if((pos-1 in bact_dict[family])and bact_dict[family][pos-1]>10*bact_dict[family][pos] and bact_dict[family][pos-1]>10):
				#print('break at '+family+'.'+str(pos))
				species_info['rough_mappings']+=1
			elif((pos+1 in bact_dict[family])and bact_dict[family][pos+1]>10*bact_dict[family][pos] and bact_dict[family][pos+1]>10):
				#print('break at '+family+'.'+str(pos))
				species_info['rough_mappings']+=1
			'''
			Raw threshold comparison
			if((bact_dict[family][pos]>100)and(pos-1 in bact_dict[family])and bact_dict[family][pos-1]<100):
				#print('break at '+family+'.'+str(pos))
				species_info['rough_mappings']+=1
			elif((bact_dict[family][pos]>100)and(pos+1 in bact_dict[family])and bact_dict[family][pos+1]<100):
				#print('break at '+family+'.'+str(pos))
				species_info['rough_mappings']+=1
			elif(bact_dict[family][pos]>100):
				species_info['smooth_mappings']+=1
			elif(bact_dict[family][pos] >= 1 and bact_dict[family][pos]<10):
				species_info['low_read_mappings']+=1

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


def check_adjacents():
	start_path = '../data_out/pool_Apair_res/genes/output/'
	species_reads = {}
	species_info = {}
	load_all_bacteria_infile_intodict(start_path,species_reads,species_info)
	for spec in species_reads:
		check_broken(species_reads[spec],species_info[spec],spec)
		#print("num families: "+str(species_info[spec]['families_count']))
		#print("num genes: " + str(species_info[spec]['genes_count']))
		#print("smooth: " + str(species_info[spec]['smooth_mappings']))
		#print("rough: " + str(species_info[spec]['rough_mappings']))
		#print("low count: " + str(species_info[spec]['low_read_mappings']))
		#print("total: " + str(species_info[spec]['total_mapped_reads']))
		species_info[spec]['percent_rough'] = species_info[spec]['rough_mappings']/(species_info[spec]['smooth_mappings']+species_info[spec]['rough_mappings'])
		species_info[spec]['percent_smooth'] = species_info[spec]['smooth_mappings']/(species_info[spec]['smooth_mappings']+species_info[spec]['rough_mappings'])
		#print("percent rough: " + str(species_info[spec]['percent_rough']))
	return species_reads,species_info


check_adjacents()

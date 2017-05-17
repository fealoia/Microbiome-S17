from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import sys
import gzip
import collections
import os
import matplotlib.pyplot as plt
import pickle

res_path = '/Users/williampascucci/Documents/Research/MIDAS/data_out/pool_Apair_res/main_genes/output'

target_genera = []
target_species = []

num_valid_swaps = 0

#builds target genera for matching potentials to see if a BLAST output is relevant
def build_genera():
	for name in os.listdir(res_path):
		name_arr = name.split('_')
		if not name_arr[0] in target_genera:
			target_genera.append(name_arr[0])

#builds target species for matching blast alignments for relevant swaps
def build_species():
	for name in os.listdir(res_path):
		name_arr = name.split('_')
		if not name_arr[0] in target_species:
			target_species.append(name_arr[0]+name_arr[1])

relation_dict = {}

#Build a simple dict in memory between species and the PATRIC gene ids
def build_relations():
	for name in os.listdir(res_path):
		name_arr = name.split('_')
		with gzip.open(res_path+'/'+name,'rt') as fn:
			for line in fn.readlines():
				if(not line.startswith('gene_id')):
					gid = line.split('.')[0]
					if(name_arr[0]+name_arr[1] in relation_dict):
						if(not gid in relation_dict[name_arr[0]+name_arr[1]]):
							relation_dict[name_arr[0]+name_arr[1]].append(gid)
					else:
						relation_dict[name_arr[0]+name_arr[1]] = [gid]

#reads the XML blast output using biopython
def read_blast(fn):

	res_handle = open(fn,"r")

	blast_record = NCBIXML.parse(res_handle)
	#print(blast_record)
	#normal_bast_res = 0
	#unexpected_blast_res = 0

	alt_dict = {}
	first = True
	fam = ''
	pickle_name = ''
	
	for record in blast_record:
		if(first):
			arr = record.query.split('.')
			fam = arr[0]
			pickle_name = "alts_"+fam+".pkl"
			if("alts_"+fam+".pkl" in os.listdir()):
				alt_dict = pickle.load(open(pickle_name,"rb"))
			first = False

		alt_dict[record.query] = []
		cter = 0
		
		for al in record.alignments:
			cter+=1
			should_print = False
			hsp_list = []
			for hsp in al.hsps:
				length = hsp.align_length
				score = hsp.score
				if(hsp.align_length and record.query_length and hsp.align_length<.9*record.query_length):
					break
				else:
					should_print = True
					alt_dict[record.query].append((al.title,hsp.query))

	pickle.dump(alt_dict,open(pickle_name,"wb"))

	res_handle.close()

	return

#iterates over all the xml files in the directory (a lot)
def web_resp_into_dict():
	for fn in os.listdir():
		if(fn.startswith('blastpkl_res_rough_reads_latest_')):
			read_blast(fn)
	return

#checks the valid swaps from the blast outputs, returns the total number of valid swaps
def check_valid_swaps(num_valid_swaps):
	for fn in os.listdir():
		if(fn.startswith('alts_')):
			alts_dict = pickle.load(open(fn,'rb'))
			for source_id in alts_dict.keys():
				for k in relation_dict.keys():
					if(source_id.split('.')[0] in relation_dict[k]):
						source_bact = k
				if(alts_dict[source_id]):
					for pos in alts_dict[source_id]:
						arr = pos[0].split(' ')
						targ_bact = arr[1]+arr[2]

						if(source_bact!=targ_bact and targ_bact in target_species):
							num_valid_swaps += 1
	return

#build_genera()
build_species()
#web_resp_into_dict()
build_relations()
check_valid_swaps(num_valid_swaps)

print(num_valid_swaps)
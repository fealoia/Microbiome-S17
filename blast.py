from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import sys
import gzip
import collections
import os
import matplotlib.pyplot as plt
import pickle

blast_in = {}

#bl_str = '>640131.3.peg.3667\nATGTCTAAGATTAAAGGTAACGTTAAGTGGTTTAATGAGTCCAAAGGATTCGGTTTCATTACTCCGGAAGATGGCAGCAAAGATGTATTCGTACACTTCTCTGCAATCCAGTCCAACGGTTTCAAAACTCTGGCTGAAGGTCAGCGTGTAGAGTTCGAAATCACTAACGGTGCCAAAGGCCCTTCTGCTGCAAACGTAATGGCTATCTAA\n>1328374.3.peg.1050\nATGAAATATGCAGCTTTAGTGGCAGGTGCAGCATTGTTTTTAACGGCGAACGCCTTCTCCGCTGAGATTTTAACGAAAGAAGCCTTTAATAAGGTCCATACGCAATATACCAAAATCGGCACCATTTCCTCTACCGGCCAGACTGCGCCGAGCGATGCCCGTGAAGAGCTGATTAAAAAGGCCGATGAAAAAGGCGCGGACGTTATCGTGCTCTCTTCCGGTAATACCGATAAAAAACTGCACGTCTCCGCTAACATCTATAAGAAGAAGTAA'

def find_alts(id_string,spec_name):
	alt_path = '../../midas_db_v1.2/pan_genomes/'+spec_name+'/centroids.ffn.gz'
	targ_seq = ''
	save = False
	with gzip.open(alt_path,'rt') as fn:
		for line in fn.readlines():
			if(save and line[0]=='>'):
				break
			if(save):
				targ_seq += line.strip()
			if(id_string in line):
				save = True
	blast_in[id_string] = targ_seq
	return targ_seq


def blast_n(blast_str,name,i):
	res_handle = NCBIWWW.qblast('blastn','nt',blast_str)

	#blast_record = NCBIXML.read(res_handle)
	with open("blastpkl_res_rough_reads_kleb_latest_"+name+'_'+str(i)+"_genes.xml","w") as out_file:
		out_file.write(res_handle.read())

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
				bact_dict[family_id][int(gene_id[3])] = {}
				bact_dict[family_id][int(gene_id[3])]['num'] = int(info_list[1])
				bact_dict[family_id][int(gene_id[3])]['alts'] = []
			else:
				bact_dict[family_id] = collections.OrderedDict()
				bact_dict[family_id][int(gene_id[3])] = {}
				bact_dict[family_id][int(gene_id[3])]['num'] = int(info_list[1])
				bact_dict[family_id][int(gene_id[3])]['alts'] = []

def check_broken(bact_dict,species_info,spec_name,fix):
	#y = []
	for family in bact_dict.keys():
		for pos in bact_dict[family].keys():
			if((bact_dict[family][pos]['num']>100)and(pos-1 in bact_dict[family])and bact_dict[family][pos-1]['num']<100):
				#print('break at '+family+'.'+str(pos))
				if(fix):
					find_alts(str(family)+'.peg.'+str(pos),spec_name)
			elif((bact_dict[family][pos]['num']>100)and(pos+1 in bact_dict[family])and bact_dict[family][pos+1]['num']<100):
				#print('break at '+family+'.'+str(pos))
				if(fix):
					find_alts(str(family)+'.peg.'+str(pos),spec_name)
			elif(bact_dict[family][pos]['num']>100):
				continue
			elif(bact_dict[family][pos]['num'] >= 1 and bact_dict[family][pos]['num']<100):
				continue

def load_all_bacteria_infile_intodict(start_path,species_dict,species_info):
	#print(os.listdir(start_path))
	for bact in os.listdir(start_path):
		species_info[bact] = {}
		species_dict[bact] = {}
		load_file(start_path+bact,species_dict[bact])

def check_adjacents():
	start_path = '../../data_out/pool_Apair_res/main_genes/output/'
	species_reads = {}
	species_info = {}
	load_all_bacteria_infile_intodict(start_path,species_reads,species_info)
	for spec in species_reads:
		spec_name = spec.split('.')[0]
		check_broken(species_reads[spec],species_info[spec],spec_name,True)
		pickle.dump(blast_in,open('rough_reads_seqs_' + spec_name + '.pkl','wb'))
	return species_reads,species_info

def write_dict():
	with open("blast_these_genes.txt","w") as fn:
		for key in blast_in.keys():
			fn.write(">"+key)
			fn.write(blast_in[key])

#pickle.dump(blast_in,open('rough_reads_seqs.pkl','wb'))
#print(blast_in)
#blast_in = pickle.load(open('rough_reads_seqs.pkl','rb'))
#blast_str = ''


def main_blast(start_ct):
	blast_str=''
	cter = 0
	for f in os.listdir():
		if(f.startswith('rough')):
			if(f.startswith('rough_reads_seqs_Bacteroides')):
				continue
			name = f[17:len(f)-4] 
			blast_in = pickle.load(open(f,'rb'))
			for key in sorted(blast_in.keys()):
				print(cter)
				if(cter<start_ct):
					cter+=1
					#print('notyet')
				else:
					arr = key.split('.')
					cter+=1
					#if(arr[0] in kleb_ref):
					blast_str+='>'+key+'\n'
					blast_str+=blast_in[key]+'\n'
					if(cter%50==0):
						blast_n(blast_str,name,cter)
						blast_str=''

main_blast(34150)
#check_adjacents()
#blast_n(blast_str)


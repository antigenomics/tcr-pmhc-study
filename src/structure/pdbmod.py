from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB import Select
from Bio.PDB import PPBuilder
from Bio.PDB import Polypeptide
import pandas as pd
import numpy as np
import subprocess
from math import ceil

PEPTIDE_BORDER = 30

class Interaction:
	def __init__(self, table, pdb):
		table = table.reset_index(drop=True)
		struct = PDBParser().get_structure(table['pdb_id'][0], pdb)
			
		table.insert(table.columns.get_loc('tcr_chain'), 'tcr_chain_name', '')
		table.loc[table['tcr_v_allele'].str.startswith('TRA'), 'tcr_chain_name'] = 'alpha'
		table.loc[table['tcr_v_allele'].str.startswith('TRB'), 'tcr_chain_name'] = 'beta'
		if sum(table['tcr_chain_name'].notnull()) != table.shape[0]:
			self.printerr('SOMETHING IS WRONG WITH V REGION ANNOTATION')
			
		table = table.fillna('')
		
		self.__table = table
		self.__name = str(struct.get_id())		
		self.__struct = struct
		self.__chains = [chain.get_id() for chain in struct[0]]
		self.__regions = table.groupby(['tcr_chain_name', 'tcr_region'])
	
		# Dictionary of regions residues;
		# looks like : { ('alpha', 'CDR1') : [residue list],
		# ('alpha', 'CDR2') : [residue list], ... } 
		self.__regions_res = self.__regions.groups 	
		for key in self.__regions_res.keys():	
			self.__regions_res[key] = []
	
		# Pepdide residue list	
		self.__peptide = []				
	
		# Dictionaries with pairwise region matrices;
		# look like : { (('alpha', 'CDR1'), ('alpha', 'CDR2')) : dataframe, 
		# (('alpha', 'CDR1'), ('alpha', 'CDR3')) : dataframe, ... }
		self.__d_matrices = {}	
		self.__e_matrices = {}	
		self.verbose = True
		if not self.getRegionsResidues():
			print 'SOME REGION WAS NOT FOUND IN PDB'
		if not self.definePeptideChain():
			print 'PEPTIDE WAS NOT FOUND IN PDB'
		
	def __str__(self):
		s = '\nprotein name: ' + self.__name + \
			'\npeptide: ' + self.getClearPeptideSeq() + \
			'\ncdr3 alpha: ' + self.__cdr3_a_seq + \
			'\ncdr3 beta: ' + self.__cdr3_b_seq + '\n'
		return s
		
	def printerr(self, err):
		if self.verbose:
			print err
		
	def getSeqLocation(self, seq):		# return sequence position and chain id
		ppb=PPBuilder()
		bltpep = ppb.build_peptides(self.__struct[0])
		for pp in bltpep: 
			beg = 0
			end = 0
			s = str(pp.get_sequence())
			ind = s.find(seq, 0, len(s))
			if (ind != -1):
				beg = beg + ind
				end = beg + len(seq) - 1
				chain = pp[0].get_parent().get_id()
				break	
		if beg == end == 0:	
			line = '\n' + seq + ' not found in '+str(self.__struct.get_id()) + '!\n'
			self.printerr(line)
			return None, None, None
		return beg, end, chain
		
	def getRegionsResidues(self):		# fill self.__regions_res dictionary with list of residues
		ppb=PPBuilder()			# for every region contained in self.__regions_res
		res = []
		bltpep = ppb.build_peptides(self.__struct[0])
		for key in self.__regions_res:
			for pp in bltpep: 
				s = str(pp.get_sequence())
				reg_seq = list(self.__regions.get_group(key)['tcr_region_seq'])[0]
				ind = s.find(reg_seq, 0, len(s))
				if (ind != -1):
					for i in range(ind, ind + len(reg_seq)):
						res.append(pp[i])
					self.__regions_res[key] = res
					break
			if not res:	
				line = '\n' + reg_seq + ' not found in '+ self.__name + '!\n'
				self.printerr('getRegionResidues(): ' + line)
				return 0
			res = []
		return 1
				
	def definePeptideChain(self):		# find peptide chain if not stated in self.__table and fill self.__peptide
		l = 'INFINITY'			# with list of peptide residues if length is less than 30
		if not self.__table['chain_antigen'][0]:
			for i in self.__chains:
				buf = len(PPBuilder().build_peptides(self.__struct[0][i])[0])
				if (buf <= l):
					l = buf
					chid = i
			self.__table.loc['chain_antigen', :] = chid
		else:
			chid = self.__table['chain_antigen'][0]
					
		pp = list(self.__struct[0][chid])
		
		if (len(pp) > PEPTIDE_BORDER):
			line = self.__name + '\t;TOO MANY AMINO ACIDS (' + str(len(pp)) + ') TO BE A PEPTIDE :(\n'
			self.printerr('definePeptideChain(): ' + line)
			return 0
				
		pep_res = []
		for r in pp:
			if (Polypeptide.is_aa(r.get_resname(), standard=True)):
				pep_res.append(r)
		self.__peptide = pep_res
		self.__regions_res.update({'peptide':pep_res})
		return 1
						
	def calcDistMatrices(self, key1, key2): 	# calculate and store distance matrix to self.__d_matrices for a pair of regions
		res_list_1 = self.getRegion(key1)	# key1 and key2 refer to keys of self.__regions_res dictionary
		res_list_2 = self.getRegion(key2)
		
		if not res_list_1 or not res_list_2:
			self.printerr('calcDistMatrices(): RESIDUE LIST IS EMPTY\n')
			return 0
			
		values = []
		for res1 in res_list_1:
			values.append([])	
			for res2 in res_list_2:	
				values[len(values)-1].append(residuesMinDist(res1, res2))
				
		rows = [Polypeptide.three_to_one(x.get_resname()) for x in res_list_1]	
		cols = [Polypeptide.three_to_one(x.get_resname()) for x in res_list_2]
		mat = pd.DataFrame(values, index = rows, columns = cols)
		self.__d_matrices.update({(key1, key2): mat})
		return 1
			
	def calcNrgMatrices(self, gromacs_dir, xpm_dir, *keys): 	# calculate and store energy matrices to self.__e_matrices for pairs of 
		res_list = [self.getRegion(key) for key in keys]	# regions; keys is a list of keys from self.__regions_res dictionary
		if not all(res_list):
			self.printerr('calcNrgMatrices: RESIDUE LIST IS EMPTY\n')
			return 0
	
		# Create pdb containing only desirable regions
		self.pushToPDB(gromacs_dir, *keys)
		
		# Execute Gromacs script
		list_of_seqs = [''.join(map(lambda res: Polypeptide.three_to_one(res.get_resname()), x)) for x in res_list]
		bashCommand = " ".join(["./grom_script.sh", gromacs_dir, xpm_dir, self.__name] + list_of_seqs)
		exitcode = subprocess.call(bashCommand, shell=True)
		if not exitcode:
			self.printerr('calcNrgMatrices: GROMACS FAILURE\n')
			return 0
		
		# Extract DataFrames from xpm files
		bad_mat, annotation = extractXPM(''.join([xpm_dir, '/', 'total', self.__name, '.xpm']))
		if len(annotation) > 2: 
			aminos = reduce(lambda x, y: list(x) + list(y), annotation[1:])
		else:
			aminos = list(annotation[-1])
		mat = pd.DataFrame(bad_mat, columns = aminos, index = aminos)
		
		ind_list = [0]
		[ind_list.append(len(x) + ind_list[-1]) for x in annotation[1:]]
		
		sub_mats = []
		sub_mats_keys = []
		for i in range(len(res_list)):
			for j in range(i, len(res_list)):
				sub_mats_keys.append((keys[i], keys[j]))
				sub_mats.append(mat.iloc[ind_list[i]:ind_list[i + 1], ind_list[j]:ind_list[j + 1]])
				
		self.__e_matrices.update(dict(zip(sub_mats_keys, sub_mats)))
		return 1
		
	def getClearPeptideSeq(self):
		if not self.__peptide:
			self.printerr('getClearPeptideSeq(): PEPTIDE (' + self.__name +') IS EMPTY\n')
			return 0
		s = ''
		for r in list(self.__peptide):
			s = s + Polypeptide.three_to_one(r.get_resname())
		return s
		
	def getCDR3AlphaSeq(self):
		return list(self.__regions.get_group(('alpha', 'CDR3'))['tcr_region_seq'])[0]
		
	def getCDR3BetaSeq(self):
		return list(self.__regions.get_group(('beta', 'CDR3'))['tcr_region_seq'])[0]
	
	def getRegion(self, key):
		try:
			res = self.__regions_res[key]
		except KeyError:
			self.printerr('NO ENERGY MATRIX FOR SUCH KEY: ' + key)
			return []
		else:
			return res
	
	def getDistMatrix(self, key):
		try:
			res = self.__d_matrices[key]
		except KeyError:
			self.printerr('NO DISTANCE MATRIX FOR SUCH KEY: ' + key)
			return pd.DataFrame()
		else:
			return res
		
	def getNrgMatrix(self, key):
		try:
			res = self.__e_matrices[key]
		except KeyError:
			self.printerr('NO ENERGY MATRIX FOR SUCH KEY: ' + key)
			return pd.DataFrame()
		else:
			return res
	
	def writeInFile_CDR3_Pept_Dist(self, path = '.'):
		if self.getDistMatrix(('peptide', ('alpha', 'CDR3'))).empty or self.getDistMatrix(('peptide', ('beta', 'CDR3'))).empty:
			self.printerr('writeInFile_CDR3_Pept_Dist: DISTANCE MATRIX IS EMPTY\n')
			return 0
			
		f = open(path + '/' + self.__name+'(0).txt', 'w')
		self.getDistMatrix(('peptide', ('alpha', 'CDR3'))).to_csv(f, sep = '\t')
		f.close()
		f = open(path + '/' + self.__name+'(1).txt', 'w')
		self.getDistMatrix(('peptide', ('beta', 'CDR3'))).to_csv(f, sep = '\t')
		f.close()
		return 1
			
	def writeInFile_CDR3_Pept_Nrg(self, path = '.'):
		if self.getNrgMatrix(('peptide', ('alpha', 'CDR3'))).empty or self.getNrgMatrix(('peptide', ('beta', 'CDR3'))).empty:
			self.printerr('writeInFile_CDR3_Pept_Nrg: ENERGY MATRIX IS EMPTY\n')
			return 0
			
		f = open(path + '/' + self.__name+'(0).txt', 'w')
		self.getNrgMatrix(('peptide', ('alpha', 'CDR3'))).to_csv(f, sep = '\t')
		f.close()
		f = open(path + '/' + self.__name+'(1).txt', 'w')
		self.getNrgMatrix(('peptide', ('beta', 'CDR3'))).to_csv(f, sep = '\t')
		f.close()
		return 1
			
	def ResiduesToStr(self):
		if not self.__peptide or not self.__regions_res[('alpha', 'CDR3')] or not self.__regions_res[('beta', 'CDR3')]:
			self.printerr('ResiduesToStr: RESIDUE LIST IS EMPTY\n')
			return 0
			
		s = ''	
		for r in self.__peptide:
			s = s + str(r.get_id())+' '
		s = s + '\n'
		for r in self.__regions_res[('alpha', 'CDR3')]:
			s = s + str(r.get_id())+' '
		s = s + '\n'
		for r in self.__regions_res[('beta', 'CDR3')]:
			s = s + str(r.get_id())+' '
		s = s + '\n\n'
		return s
		
	def printIndices(self, f):
		f.write(self.indicesToStr())
		
	class _ResSelect(Select):
	    def __init__(self, *res_lists):
	    	self.__res_list = reduce(lambda x, y: x + y, res_lists)
	    def accept_residue(self, residue):
		if residue in self.__res_list:
		    return 1
		else:
		    return 0
	
	def pushToPDB(self, path, *keys):
		seq_list = [self.getRegion(key) for key in keys]
		if not all(seq_list):
			self.printerr('pushToPDB: RESIDUE LIST IS EMPTY\n')
			return 0

		io = PDBIO()
		io.set_structure(self.__struct)
		io.save(path + '/' +self.__name + ".pdb", self._ResSelect(*seq_list))
		return 1
	
	def overallTable(self, main_key, *keys):
		seq_list = [self.getRegion(key) for key in keys] + self.getRegion(main_key)
		if not all(seq_list):
			self.printerr('overallTable: RESIDUE LIST IS EMPTY\n')
			return 0
		
		table = self.__table.copy()	
		table.drop(['chain_mhc_a', 
		    'chain_mhc_b', 
		    'chain_mhc_b', 
		    'tcr_region_start', 
		    'tcr_region_end', 
		    'chain_antigen'], axis=1, inplace=True)
		
		table.insert(table.shape[1], 'len_tcr_region', '')
		table.insert(table.shape[1], 'pos_tcr_region', '')
		table.insert(table.shape[1],  'aa_tcr_region', '')
		
		region_name = str(main_key)
		table.insert(table.shape[1], 'len_' + region_name, '')
		table.insert(table.shape[1], 'pos_' + region_name, '')
		table.insert(table.shape[1],  'aa_' + region_name, '')
		table.insert(table.shape[1], 'nrg_' + region_name, '')
		table.insert(table.shape[1], 'dst_' + region_name, '')
		
		for key in keys:
			cond = (table['tcr_chain_name'] == key[0]) & (table['tcr_region'] == key[1])
			cut_index = table[cond].index.tolist()[0]
			
			nrg_table = self.getNrgMatrix((main_key, key))
			dst_table = self.getDistMatrix((main_key, key))
			if nrg_table.empty or dst_table.empty:
				continue
		
			nrg_table_mod 		= nrg_table.copy()
			nrg_table_mod.index 	= range(len(nrg_table.index))
			nrg_table_mod.columns 	= range(len(nrg_table.columns))
			
			dst_table_mod 		= dst_table.copy()
			dst_table_mod.index 	= range(len(dst_table.index))
			dst_table_mod.columns 	= range(len(dst_table.columns))
			
			nrg_melt_table 		= pd.melt(nrg_table.reset_index(), id_vars = ['index'])
			dst_melt_table		= pd.melt(dst_table.reset_index(), id_vars = ['index'])
			nrg_melt_table_mod 	= pd.melt(nrg_table_mod.reset_index(), id_vars = ['index'])
			dst_melt_table_mod 	= pd.melt(dst_table_mod.reset_index(), id_vars = ['index'])
			
			len_tcr_region = len(self.getRegion(key))
			len_region_name =  len(self.getRegion(main_key))
			table.loc[cond, 'len_tcr_region'] = len_tcr_region
			
			small = pd.concat(nrg_melt_table.shape[0]*[table[cond]])
			small.loc[:, 'pos_tcr_region'] 		= [x - ceil(len_tcr_region / 2) for x in nrg_melt_table_mod['variable'].tolist()]
			small.loc[:, 'pos_' + region_name] 	= [x - ceil(len_region_name / 2) for x in nrg_melt_table_mod['index'].tolist()]
			small.loc[:, 'aa_tcr_region'] 		= nrg_melt_table['variable'].tolist()
			small.loc[:, 'aa_' + region_name] 	= nrg_melt_table['index'].tolist()
			small.loc[:, 'len_' + region_name] 	= len_region_name
			small.loc[:, 'dst_' + region_name] 	= ["{:.2f}".format(x) for x in dst_melt_table['value'].tolist()]
			small.loc[:, 'nrg_' + region_name] 	= nrg_melt_table['value'].tolist()
			#small = small.reset_index(drop=True)
			#small.to_csv('sad.dat', sep='\t')

			table = pd.concat([table.iloc[:cut_index, :], small, table.iloc[cut_index + 1:, :]])
			table = table.reset_index(drop=True)
		self.summary_table = table
		return 1
		
		
	'''def convertTable(self):	
		if self.isEmpty():
			self.printerr(self.__info)
			return 0
			
		acdr3 = self.getSeqLocation(self.__cdr3_a_seq)
		bcdr3 = self.getSeqLocation(self.__cdr3_b_seq)
		#avreg = self.getSeqLocation(self.__vreg_a)
		#bvreg = self.getSeqLocation(self.__vreg_b)
		achain = acdr3[2]
		bchain = bcdr3[2]
		regnum = 4
		tb = pd.DataFrame(columns = ['complex',
			'chain',
			'type',
			'name',
			'id',
			'region',
			'start',
			'end'])
		chains = self.__chains
		achainpos = chains.index(achain)
		chains = chains[:achainpos] + regnum * [achain] + chains[achainpos + 1:]
		bchainpos = chains.index(bchain)
		chains = chains[:bchainpos] + regnum * [bchain] + chains[bchainpos + 1:]
		
		tb.loc[:, 'chain'] = chains
		tb.loc[:, 'complex'] = self.__name
		tb.loc[tb['chain'].isin(self.__mhc.keys()), 'type'] = 'MHC'
		tb.loc[tb['chain'].isin([achain, bchain]), 'type'] = 'TCR'
		tb.loc[tb['chain'] == self.__pepchain, 'type'] = 'Antigen'
		tb.loc[tb['chain'] == achain, 'name'] = 'alpha'
		tb.loc[tb['chain'] == bchain, 'name'] = 'beta'
		tb.loc[tb['chain'].isin(self.__mhc.keys()), 'name'] = list(self.__mhc.values())
		tb.loc[tb['name'] == 'alpha', 'region'] = ['cdr1', 'cdr2', 'cdr3', 'vreg']
		tb.loc[tb['name'] == 'beta', 'region'] = ['cdr1', 'cdr2', 'cdr3', 'vreg']
		tb.loc[(tb['name'] == 'alpha') & (tb['region'] == 'cdr3'), ['start', 'end']] = acdr3[:2]
		tb.loc[(tb['name'] == 'beta') & (tb['region'] == 'cdr3'), ['start', 'end']] = bcdr3[:2]
		#tb.loc[(tb['name'] == 'alpha') & (tb['region'] == 'vreg'), ['start', 'end']] = avreg[:2]
		#tb.loc[(tb['name'] == 'beta') & (tb['region'] == 'vreg'), ['start', 'end']] = bvreg[:2]
		return tb'''
		
	
def residuesMeanDist(res1, res2):
	dist = 0
	num = 0
	if (res1 == res2):
		return 0;
	else:
		for at1 in res1:
			for at2 in res2:
				dist += (at1-at2)
				num += 1
		return dist/num
	
def residuesMinDist(res1, res2):
	dist = 'INFINITY'
	if (res1 == res2):
		return 0;
	else:
		for at1 in res1:
			for at2 in res2:
				buf = at1-at2
				if (buf < dist):
					dist = buf		
		return dist
		
def residuesCaDist(res2, res1):
	if (res1 == res2):
		return 0;
	else:
		return res1['CA']-res2['CA']
		

'''def parseDataFile(path):
	f = open(path, 'r')
	l = list(f)
	d = {}
	for i in range (1, len(l)):
		s = l[i].split("\t")
		if (len(s) == 20):
			pname = s[0].lower()
			if (d.get(pname) == None):
				d.update({pname : [[s[1].strip()]]})
			else:
				d[pname][0].append(s[1].strip())
			if (s[17].isalpha()):
				d[pname].append(s[17].strip())
			else
				d[pname].append('Nan')
			if
	f.close()
	return d'''

'''	
def parseDataFile(path):
	fr = pd.DataFrame(pd.read_table(path, sep='\t'))
	fr['PDB ID'] = map(lambda x: x.lower(), fr['PDB ID']) 
	
	return fr
'''
	
'''
def getProtein(frame, name):
	smframe = frame[frame['PDB ID'] == name]
	if smframe.shape[0] == 0:
		print "Protein " + name + " is not listed in the table! Closing...\n"
		return 0
	
	# get the names of the chains listed in table
	pepchains = map(lambda x: x.strip(), smframe['Chain ID'].dropna())

	#get the dict with MHC types for each chain listed in table
	buframe = smframe[frame['HLA Type'].notnull()]
	mhc = dict(zip(buframe['Chain ID'], buframe['HLA Type']))

	#get CDR3 and V Region
	buframe = smframe[smframe['cdr3'].notnull()]

	if buframe.shape[0] != 2:
		print "Number CDR3 listed does not equal 2! Closing...\n"
		return 0

	it = buframe.iterrows()
	index, row = next(it)
	cdr3_a = [row['cdr3'].strip(), row['v'].strip()]
	index, row = next(it)
	cdr3_b = [row['cdr3'].strip(), row['v'].strip()]

	tra_crit = buframe['v'].map(lambda x: x.startswith('TRA'))
	trb_crit = buframe['v'].map(lambda x: x.startswith('TRB'))
	if any(tra_crit) and any(trb_crit):
		cdr3_a = [buframe.loc[tra_crit, 'cdr3'].values[0].strip(), buframe.loc[tra_crit, 'v'].values[0].strip()]
		cdr3_b = [buframe.loc[trb_crit, 'cdr3'].values[0].strip(), buframe.loc[trb_crit, 'v'].values[0].strip()]

	return pepchains, mhc, cdr3_a, cdr3_b
'''
	
def extractXPM(path):
	xpm = open(path, 'r')
	lines = list(xpm)
	xpm.close()
	ind = 0
	for line in lines:
		pline = line.split()
		if pline[0] == 'static':
		    break
		ind += 1

	ind += 1 
	print lines[ind].split('"')[0]
	snum = int(lines[ind].split('"')[1].split()[3])
	ind += 1 

	alphabet = {}
	while(True):
		if lines[ind][0] == '/':
		    break
		pline = lines[ind].split('"')
		code = pline[1].split()[0]
		alphabet.update([(code, pline[3])])
		ind += 1
	ind += 2

	mat = []
	while(True):
		if lines[ind][0] == '/':
		    break
		line = lines[ind].split('"')[1]
		l = []
		ins = ''
		for i in range(0, len(line)):
		    ins += line[i]
		    if (i+1)%snum == 0:
			l.append(float(alphabet[ins]))
			ins = ''
		mat.append(l)
		ind += 1
		
	lline = lines[ind]
	annot = lline.replace('/', ' ').replace('*', ' ').split()
	mat = list(reversed(mat))
		
	'''a_mat = []
	a_bs = [len(annot[1]), len(annot[1]) + len(annot[2])]
	for line in mat[:len(annot[1])]:
		a_mat.append(line[a_bs[0]:a_bs[1]])
		
	b_mat = []
	b_bs = [a_bs[1], a_bs[1] + len(annot[3])]
	for line in mat[:len(annot[1])]:
		b_mat.append(line[b_bs[0]:b_bs[1]])'''
		
	return mat, annot

'''def writeInFile_CDR3_CDR3(parser, datadist, item, f):
	structure = parser.get_structure(item, '../pdbs/'+item+'.pdb')
	cdr3seqlist = datadist[item]
	cdr3reslist = []
	for i in range(1, len(cdr3seqlist)):
		cdr3reslist.append(getCDR3Indices(structure, cdr3seqlist[i]))
	i = 0
	for cdr3res in cdr3reslist:
		f2 = open('generated/pdbcdr3/dist_mats/cdr3+cdr3/'+item+'('+str(i)+').txt', 'w')
		m = calcDistMatrix(item, cdr3res, cdr3res, True)
		fwriteMatrix(f2, m)
		f2.write('\n')
		f2.close()
		f.write(item+'\n')
		fwriteMatrix(f, m)
		f.write('\n')
		i += 1'''

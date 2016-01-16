import pdbmod
import pandas as pd
import sys
import argparse
from os.path import isfile
from Bio.PDB import PDBParser

#==================================================================================
# For calculating energies grom_script.sh must be in the same folder with this file
#
#			    
# Usage example: pdbstat.py table -n name1.pdb name2.pdb name3.pdb
# Use --help for more info
#
# this will create summary table in the folder specified by result_table_dir
# with the regions listed in the variables main_region and regions
#==================================================================================
# variants for final flag -a: nrg, dist, comp, table
#==================================================================================

pdb_dir = '../../../../pdbs'				# where .pdb files lie
data_path = '../../../../pdbs/final.annotations.txt'	# annotation file
nrg_dir = '../../../generated/pdbcdr3/energy_mats' 	# where to write nrg matrices
dist_dir = '../../../generated/pdbcdr3/dist_mats'	# where to write dist matrices
gromacs_dir = '../gromacs'				# sandbox for gromacs actions
xpm_dir = '../gromacs/xpm'				# folder with .xmp files (result of gromacs activity)
result_table_dir = '.'					# where to write summary table

# specify regions you want to calculate intercation between (for Table(name))
main_region = 'peptide' # other variants of region names : ('alpha', 'CDR3') or ('alpha', 'CDR1')
			# or ('beta', 'CDR1') etc. according to the names in the final.annotations table
			
	    # first region      second region     enter as much region as you wish and launch
regions = [('alpha', 'CDR3'), ('beta', 'CDR3')] 

#==================================================================================

def Nrg(name): # store a couple of energy matrices for peptide|CDR3alpha (0) and peptide|CDR3beta (1) to nrg_dir
	print "\nProtein name: " + name + "\n"
	info = data[data['pdb_id'] == name]
	pdb_file = pdb_dir + '/' + name + '.pdb'
	c = pdbmod.Interaction(info, pdb_file)
		
	if not c.calcNrgMatrices(gromacs_dir, xpm_dir, 'peptide', ('alpha', 'CDR3'), ('beta', 'CDR3')):
		return 0
	if not c.writeInFile_CDR3_Pept_Nrg(nrg_dir):
		return 0
	return c

def Dist(name):	# store a couple of distance matrices for peptide|CDR3alpha (0) and peptide|CDR3beta (1) to dist_dir
	print "\nProtein name: " + name + "\n"
	info = data[data['pdb_id'] == name]
	pdb_file = pdb_dir + '/' + name + '.pdb'
	c = pdbmod.Interaction(info, pdb_file)
		
	if not c.calcDistMatrices('peptide', ('alpha', 'CDR3')):
		return 0
	if not c.calcDistMatrices('peptide', ('beta', 'CDR3')):
		return 0
	if not c.writeInFile_CDR3_Pept_Dist(dist_dir):
		return 0
	return c
	
def CompressPDB(name): # store truncated pdbs contaning only specified regions to ../truncpdbs
	print "\nProtein name: " + name + "\n"
	info = data[data['pdb_id'] == name]
	pdb_file = pdb_dir + '/' + name + '.pdb'
	c = pdbmod.Interaction(info, pdb_file)
		
	if not c.pushToPDB('../../../../truncpdbs', 'peptide', ('alpha', 'CDR3'), ('beta', 'CDR3')):
		return 0
	return c
	
def Table(name): # create big table with pairwise energy and distance
	print "\nProtein name: " + name + "\n"
	info = data[data['pdb_id'] == name]
	pdb_file = pdb_dir + '/' + name + '.pdb'
	c = pdbmod.Interaction(info, pdb_file)
		
	for region in regions:
		if not c.calcDistMatrices(main_region, region):
			return 0
	if not c.calcNrgMatrices(gromacs_dir, xpm_dir, main_region, *regions):
		return 0
	if not c.overallTable(main_region, *regions):
		return 0
	return c
	
def Iterate(items, FUN):	# run one of upper functions for several proteins
	output = [FUN(i) for i in items]
	return filter(lambda x: x != 0, output)
		
	
def Err():
	print 'Wrong flag (nrg/dist/comp/table)'
	exit(1)
	
def CheckInput(data, items):
	table_names_set = set(data['pdb_id'].tolist())
	input_names_set = set(items)

	to_download_list = [item for item in items if not isfile(pdb_dir + '/' + item + '.pdb')]
	input_minus_table = input_names_set.difference(table_names_set)
	f = open('pdbstat.log', 'w')
	f.write('These should be downloaded:' + str(to_download_list) + '\n\n')
	f.write('These should are not present in the table:' + str(input_minus_table) + '\n')
	f.close()
	items = list(table_names_set.intersection(input_names_set).difference(set(to_download_list)))
	return items

#============================ Main =============================================

parser = argparse.ArgumentParser(description='Use --help for more info')

parser.add_argument('--names', '-n', nargs='+', dest='names', type=str, default=['all'], \
	help='type a list of pdb files (ex: 1ao7.pdb 1awx.pdb etc.) or "all" to process all listed in the table at once (default: all)')
parser.add_argument('action', nargs=1, type=str,
	choices=['nrg', 'dist', 'comp', 'table'], help=str('what to do with the pdbs;\n \
		nrg: puts pairwise energy matrices into folder specified by nrg_dir;\n \
		dist: puts pairwise distance matrices into folder specified by dist_dir;\n \
		comp: puts compressed pdb files into folder specified by pdb_dist;\n \
		table: calculates overall table\n'))
args = parser.parse_args()

switch = {'nrg': lambda x: Iterate(x, Nrg),
	'dist': lambda x: Iterate(x, Dist),
	'comp': lambda x: Iterate(x, CompressPDB), 
	'table': lambda x: Iterate(x, Table)}
exception = lambda x : Err();

flag = args.action[0]
names = args.names

FUN = switch.get(flag, exception)
data = pd.DataFrame(pd.read_table(data_path, sep='\t'))

if len(names) == 1 and names[0] == 'all':
	names = list(set(data['pdb_id'].tolist()))
else:
	names = map(lambda x: x[-8:-4], names)

names = CheckInput(data, names)

cpx = FUN(names)
if not cpx:
	print 'ERROR OCCURED!'
	
try:
	if flag == 'table':
		result_table = pd.concat([x.summary_table for x in cpx])
		result_table = result_table.reset_index(drop=True)
		result_table.to_html(result_table_dir + '/summary_table.html')
		result_table.to_csv(result_table_dir + '/summary_table.txt', sep = '\t')
except ValueError:
	print 'NOTHING TO RETURN'





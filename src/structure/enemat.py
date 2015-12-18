
# coding: utf-8

# In[237]:

from Bio.PDB import PDBParser
from Bio.PDB import Polypeptide
from gromacs import cbook
import pandas as pd
import gromacs.cbook
import pdbmod
import sys


# In[238]:

def parseGro(f):
    
    fp = open(f, 'r')
    pp = list(fp)
    res = ''
    atlist = []
    aas = []
    prevnum = '0'
    for i in range(2, len(pp) - 1):
        line = pp[i].split()
        num = line[0][:-3]
        aa = line[0][-3:]
            
        if (num != prevnum):
            if not Polypeptide.is_aa(aa):
            	shortaa = 'X'
            else:
            	shortaa = Polypeptide.three_to_one(aa)
            atlist.append(i - 1)
            aas.append(aa + '_' + num)
            prevnum = num
            res = res + shortaa
            
    fp.close()
    atlist.append(len(pp) - 3)
    return res, atlist, aas


# In[239]:

def findSeqsInGro(protein_seq, seq_list):
	begin = []
	end = []
	for seq in seq_list:
		ind = protein_seq.find(seq)
		if (ind == -1):
			sys.stderr.write('SEQUENCE ' + seq + ' WAS NOT FOUND IN .gro')
			exit(1)
		begin.append(ind + 1)
		end.append(begin[-1] + len(seq) - 1)
    	return zip(begin, end)


# In[240]:

def appendSeqs(path, borders_list, atlist, aas):
    fl = open(path, 'a')
    groups = gromacs.cbook.get_ndx_groups(path)
    names = []
    for pair in borders_list:
        for aaid in range(pair[0] - 1, pair[1]):
            name = 'r_' + str(pair[0]) + '-' + str(pair[1]) + '_' + aas[aaid]
            names.append(name)
            if any(group['name'] == name for group in groups):
                continue
            fl.write('[ ' + name + ' ]\n')
            for item in range(atlist[aaid], atlist[aaid + 1]):
                fl.write(str(item) + ' ')
            fl.write('\n')
    fl.close()
    return names


# In[241]:

def appendToMDP(path, grnames):
    replaceText = 'energygrps\t='
    for name in grnames:
        replaceText = replaceText + ' ' + name
    replaceText = replaceText + '\n'
    
    fl = open(path, 'r')
    text = fl.read()
    fl.close()
    
    fl = open(path, 'r')
    s = ''
    for line in fl:
        if line.split()[0] == 'energygrps':
            s = line
            break
    fl.close()
    
    fl = open(path, 'w')
    if not s:
        fl.write(text + replaceText)
    else:
        fl.write(text.replace(s, replaceText))
    fl.close()


# In[242]:

def appendToGroupsDat(path, names):
    fl = open(path, 'w')
    fl.write(str(len(names)) + '\n')
    for name in names:
        fl.write(name + '\n')


# In[243]:

def appendSeqs(path, idlist, atlist, aas):
	fl = open(path, 'a')
	groups = gromacs.cbook.get_ndx_groups(path)
	names = []
	for pair in idlist:
		for aa in range(pair[0] - 1, pair[1]):
		    name = 'r_' + str(pair[0]) + '-' + str(pair[1]) + '_' + aas[aa]
		    names.append(name)
		    if any(group['name'] == name for group in groups):
			continue
		    fl.write('[ ' + name + ' ]\n')
		    for atom in range(atlist[aa], atlist[aa + 1]):
			fl.write(str(atom) + ' ')
		    fl.write('\n')
	fl.close()
	return names

if len(sys.argv) < 5:
	sys.stderr.write('NEED MORE PARAMETERS TO EXECUTE')		
	exit(1)

item = sys.argv[3]
	
index_file_dir = sys.argv[1]
mdp_file_dir = index_file_dir + '/params'
G = gromacs.cbook.IndexBuilder(index_file_dir + '/' + item + '.gro')
G.cat(index_file_dir + '/index.ndx')
st, atlist, aas = parseGro(index_file_dir + '/' +item + '.gro')
seq_list = sys.argv[4:]

borders_list = findSeqsInGro(st, seq_list)
names = appendSeqs(index_file_dir + '/index.ndx', borders_list, atlist, aas)
appendToMDP(mdp_file_dir + '/minim.mdp', names)
appendToGroupsDat(index_file_dir + '/groups' + item + '.dat', names)

print '/*' + '*/ /*'.join([item] + seq_list) + '*/'




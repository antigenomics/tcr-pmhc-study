
# coding: utf-8

# In[21]:

#get_ipython().magic('matplotlib inline')

import pandas as pd
import numpy as np
from functools import *

# Fetch structural data, refine it a bit, and assign Myazava energy value to each pair
df = pd.DataFrame(pd.read_table('structure.txt'))

df['pos_tcr'] = (df['pos_tcr'] - list(map(np.floor, (df['len_tcr'] / 2)))).apply(int)
df['pos_antigen'] = (df['pos_antigen'] - list(map(np.floor, (df['len_antigen'] / 2)))).apply(int)
df.dropna(inplace=True)
df.index = range(df.shape[0])

#Some defective structures, that shoul not be taken into account
#
#* 1u3h - to del
#* 2xn9 - tcr too far
#* 4eup - tcr shifted 
#* 4ozg - to del
#* 5d2l - tcr too far
#* 5d2n - to del
#* 4c56 - article: "..demonstrates absence of Tcr-Peptide contacts"

bad_ids = ['1u3h', '2xn9', '4eup', '4ozg', '5d2l', '5d2n', '4c56']
groups = df.groupby('pdb_id')
for name in bad_ids:
    df = df.drop(groups.get_group(name).index)

#Bad structures deleted
#3d39's antigen LLFGPVYV should be changed to LLFGPVY, as V is missed in the distance/energy columns

df.loc[df.pdb_id == '3d39', 'antigen_seq'] = 'LLFGPVY'
df.loc[df.pdb_id == '3d39', 'len_antigen'] = 7

#We should also eliminate structures having indentical both cdr3s, mhc and peptide at once

temp = df.groupby('pdb_id').apply(lambda x: [x[x.tcr_v_allele.str.startswith('TRA')]['tcr_region_seq'].iloc[0], 
                                             x[x.tcr_v_allele.str.startswith('TRB')]['tcr_region_seq'].iloc[0],
                                             x['mhc_a_allele'].iloc[0], x['mhc_b_allele'].iloc[0], x['antigen_seq'].iloc[0]])
temp = pd.DataFrame(temp.apply(str), index=temp.index)
temp.groupby(0).groups.values()
bad_ids = reduce(lambda x, y: x + y, [x[1:] for x in temp.groupby(0).groups.values() if len(x) > 1])

#So we delete these: '3mv8', '3mv9', '2ckb', '3kxf', '4jrx', '2z31', '2f53', '2p5e', '2p5w', '2pye', '2vlj', '2vlk', '2vlr', '3utt'

groups = df.groupby('pdb_id')
for name in bad_ids:
    df = df.drop(groups.get_group(name).index)
df.to_csv('structure_refined.txt', sep='\t', index=False)


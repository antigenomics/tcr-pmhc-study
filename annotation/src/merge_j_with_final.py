import pandas as pd
import numpy as np

jreg = pd.read_csv('../../result/tcr.jreg.annotations.txt', sep='\t')
['pdb_id', 'tcr_chain', 'species']

final = pd.read_csv('../../result/final.annotations.txt', sep='\t')
ind = int(np.where(final.columns == 'tcr_v_allele')[0][0])
if (not 'tcr_j_allele' in final.columns):
    final.insert(ind+1, 'tcr_j_allele', np.nan)

for ind in range(final.shape[0]):
    jname = jreg[(jreg.pdb_id == final.ix[ind, 'pdb_id']) & (jreg.tcr_chain == final.ix[ind, 'tcr_v_allele'][:3]) & (jreg.species == final.ix[ind, 'species'])].j_allele
    final.ix[ind, 'tcr_j_allele'] = '' if jname.empty else jname.iloc[0]
final.to_csv('../../result/final.annotations.txt', sep='\t', index=False)

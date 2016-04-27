import pandas as pd

df = pd.read_csv('../../result/structure.txt', sep='\t')
dfold = pd.read_csv('../../result/structure_old.txt', sep='\t')
df['energy'] = dfold['energy']
df.to_csv('../../result/new_structure.txt', sep='\t', index=False)



import pandas as pd
import numpy as np

alphabet = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', \
            'G', 'H', 'I', 'L', 'K', 'M', 'F', \
            'P', 'S', 'T', 'W', 'Y', 'V']

def empty_aa_mat():
    return pd.DataFrame(np.nan, columns=alphabet, index=alphabet)
    
def get_matrix_from_bryant(bryant):
    cols = []
    for distance in bryant.columns:
        col = bryant[distance]
        mat = empty_aa_mat()
        for i in alphabet:
            for j in alphabet:
                temp = 0
                try:
                    temp = col[i, j]
                except KeyError:
                    temp = col[j, i]
                mat.loc[i, j] = temp
        cols.append(mat.stack().to_frame(distance))
    return reduce(lambda x, y: x.join(y), cols)
    
tl = pd.read_table('bryant_table_pairwise.txt', sep='\s+')
vals = pd.concat([tl.iloc[x:x+6, :].transpose() for x in np.arange(0, 123, 7)])
vals = vals.iloc[:-3, :]
pairs = [list(tl.columns)] + [[tl.index[x]] + list(tl.iloc[x, :-1]) for x in np.arange(6, 123, 7)]
pairs = reduce(lambda x, y: x + y, pairs)[:-3]
pairs = map(lambda x: (x[0], x[2]), pairs)
pairs = pd.MultiIndex.from_tuples(pairs)
vals.index = pairs
vals = pd.DataFrame(vals.applymap(float))

get_matrix_from_bryant(vals).to_csv('bryant_table_pairwise_nice.txt')

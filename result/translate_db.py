
# coding: utf-8

# In[20]:

import pandas as pd
import numpy as np

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


# In[21]:

table = pd.read_csv('segments.all.txt', sep='\s+')
translate = lambda seq: str(Seq(seq, generic_dna).translate())
#open('segments.all.fasta', 'w').close()
#f = open('segments.all.fasta', 'a')
#table.apply(lambda x: f.write('>'+'|'.join([x.id, x.gene, x['#species']])+'\n'+translate(str(x.sequence))+'\n'), axis=1);
#f.close()

table.sequence = table.sequence.apply(translate)
table.to_csv('segments.all.prot.txt', sep='\t', index=False)


# In[ ]:




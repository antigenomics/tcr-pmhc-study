
# coding: utf-8

# In[142]:

import numpy as np
import pandas as pd
import urllib.request
import bs4
import re
import time

from google import search


# In[143]:

df = pd.read_csv('../../result/final.annotations.txt', sep='\t')
grouped = df.groupby('pdb_id')


# In[144]:

df.head()


# In[145]:

columns = ['chunk.id', 'cdr3.alpha', 'v.alpha', 'j.alpha', 'cdr3.beta', 'v.beta',
       'd.beta', 'j.beta', 'species', 'mhc.a', 'mhc.b', 'mhc.class',
       'antigen.epitope', 'antigen.gene', 'antigen.species', 'reference.id',
       'method.identification', 'method.frequency', 'method.singlecell',
       'method.sequencing', 'method.verification', 'meta.study.id',
       'meta.cell.subset', 'meta.subset.frequency', 'meta.subject.cohort',
       'meta.subject.id', 'meta.replica.id', 'meta.clone.id',
       'meta.epitope.id', 'meta.tissue', 'meta.donor.MHC',
       'meta.donor.MHC.method', 'meta.structure.id']


# In[146]:

# in grouped.groups.keys()


# In[150]:

print('Create database...')

rows = []
for key, group in grouped:
    row = pd.Series([np.nan]*len(columns), columns)
    row['chunk.id'] = key
    row['cdr3.alpha'] = group.loc[(group.tcr_v_allele.str.startswith('TRA')) & 
                                  (group.tcr_region == 'CDR3'), 'tcr_region_seq'].iloc[0]
    row['v.alpha'] = group.loc[(group.tcr_v_allele.str.startswith('TRA')), 'tcr_v_allele'].iloc[0]
    row['j.alpha'] = group.loc[(group.tcr_v_allele.str.startswith('TRA')), 'tcr_j_allele'].iloc[0]
    row['cdr3.beta'] = group.loc[(group.tcr_v_allele.str.startswith('TRB')) & 
                                 (group.tcr_region == 'CDR3'), 'tcr_region_seq'].iloc[0]
    row['v.beta'] = group.loc[(group.tcr_v_allele.str.startswith('TRB')), 'tcr_v_allele'].iloc[0]
    row['j.beta'] = group.loc[(group.tcr_v_allele.str.startswith('TRB')), 'tcr_j_allele'].iloc[0]
    row['species'] = group['species'].iloc[0]
    row['mhc.a'] = group['mhc_a_allele'].iloc[0]
    row['mhc.b'] = group['mhc_b_allele'].iloc[0]
    row['mhc.class'] = group['mhc_type'].iloc[0]
    row['antigen.epitope'] = group['antigen_seq'].iloc[0]
    row['mhc.class'] = group['mhc_type'].iloc[0]
    rows.append(row)    
    
db = pd.concat(rows, axis=1).T

# Somehow id '4e41' transformed into float 4e+41
ind = db['chunk.id'].apply(type) == float
if any(ind):
    db.loc[ind, 'chunk.id'] = '4e41'


# In[88]:

db.head()


# In[89]:

def get_pubmed_id(pdb_id):
    try:
        url = 'http://www.rcsb.org/pdb/explore.do?structureId='+pdb_id.lower()
        source = urllib.request.urlopen(url).read()
        soup = bs4.BeautifulSoup(source, 'html5lib')
        string = str(soup.find('meta', {'name':'description'}))
        pattern = re.compile('<meta content="[0-9A-Za-z]{4}:\s+(.+)" name')
        article = pattern.match(string).group(1)
    except BaseException:
        print(pdb_id.lower()+': '+"Something's wrong")
        return np.nan
    for counter in range(1):
        try:
            if counter > 0:
                time.sleep(np.random.random_integers(10, 30))
            links = [url for url in search(article, stop=40)]
            pmids = [split[-1] if split[-2] == 'pubmed' else '' for split in [link.split('/') for link in links]]
            global pm
            pm = pmids[np.where(np.array(pmids) != '')[0][0]]
        except BaseException as e:
            if str(e).find('503'):
                raise Exception('Seems google blocked you')
            continue
        else:
            break
    print(pdb_id.lower()+': '+pm)
    return pm


# In[90]:

#refs = pd.Series([get_pubmed_id(pdb_id) for pdb_id in db['chunk.id']], list(db['chunk.id']))


# In[91]:

#get_pubmed_id('1ao7')


# In[92]:

#refs = refs.apply(lambda x: 'PMID'+x)
#db['reference.id'] = refs
#db.to_csv('database.txt', sep='\t', index=False)


# In[93]:

def get_pubmed_id2(pdb_id):
    try:
        url = 'http://www.rcsb.org/pdb/explore.do?structureId='+pdb_id
        source = urllib.request.urlopen(url).read()
        soup = bs4.BeautifulSoup(source, 'html5lib')
        string = str(soup.find('meta', {'name':'description'}))
        pattern = re.compile('<meta content="[0-9A-Za-z]{4}:\s+(.+)" name')
        article = pattern.match(string).group(1)

        mkquery = lambda x: 'http://www.ncbi.nlm.nih.gov/pubmed?term=%28'+'+'.join(x.split(' '))+'[Title]%29'
        url = mkquery(article)
        print('Article url: '+url)
        source = urllib.request.urlopen(url).read()
        soup = bs4.BeautifulSoup(source, 'html.parser')
        if any(np.array([x.text for x in soup.findAll('h3')]) == 'Abstract'):
            find_pmid = lambda x: x.findAll('div', id='maincontent')[0].findAll('div', 'resc')[0].find('dd').text
            res = find_pmid(soup)
            print('PMID: '+res+'\n')
            return res
        if any(np.array([x.text for x in soup.findAll('h2')]) == 'Search results'):
            find_pmid = lambda x: x.findAll('div', id='maincontent')[0].findAll('dl', 'rprtid')[0].find('dd').text
            res = find_pmid(soup)+'?'
            print('PMID: '+res+'\n')
            return res
        print('PMID: not found\n')
        return ''
    except BaseException:
        return ''

'''
# In[94]:

print('Fetch PubMed IDs...')

refs = pd.Series([get_pubmed_id2(pdb_id) for pdb_id in db['chunk.id']], list(db['chunk.id']))

# In[ ]:

print('Filling missed values')

# In[166]:

refs[refs == '']


# In[154]:

refs['1bd2'] = '9586631'
refs['1fo0'] = '11017099'
refs['1fyt'] = '11060013'
refs['2f53'] = '16600963'
refs['2f54'] = '16600963'
refs['3d39'] = '19698083'
refs['4eup'] = ''
refs['4ftv'] = ''
refs['4e41'] = '17334368'
refs['5d2l'] = '26429912'
refs['5d2n'] = '26429912'


# In[165]:

refs[refs.apply(lambda x: x[-1] == '?' if x else False)]


# In[164]:

refs['1ao7'] = '8906788'
refs['1oga'] = '12796775'

# In[169]:

pretty_refs = refs.apply(lambda x: 'PMID:'+str(x))
'''

pretty_refs = pd.read_csv('../tmp/references.txt', sep='\t', index_col=0)
pretty_refs = pretty_refs['reference.id']

# In[172]:

db.index = db['chunk.id']
db.loc[:, 'reference.id'] = pretty_refs
db.to_csv('../../result/database.txt', sep='\t', index=False)
print('Saved to ../../result/database.txt')

import pandas as pd

t = pd.DataFrame(pd.read_table('summary_table.txt', sep = '\t'))
#t.dropna(inplace=True)
#print t
t.to_csv('summary_table.csv', index=False)

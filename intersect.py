import pandas as pd


prots = pd.read_csv('Data/proteins_with_significant_domains.tsv', sep='\t')
data = pd.read_csv('Data/Domains/expanded_data.tsv', sep='\t')

domain_series = prots['SUPFAM']

domains = []
for entry in domain_series:
    #print(entry)
    if not type(entry) == str:
        continue
    for token in entry.split(';'):
        if len(token) > 3:
            domains.append(token[3:])

to_retain = []
for i in range(data.shape[0]):
    entry = data.loc[i, 'exon_associated_superfamilies']
    if str(entry) in domains:
        to_retain.append(i)

data = data.loc[to_retain,:]
data.to_csv('exons_of_interest.tsv', sep='\t')
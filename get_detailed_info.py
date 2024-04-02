import pandas as pd
import numpy as np
from copy import deepcopy

def load_data():
    data = pd.read_csv('Data/Domains/expanded_data.tsv', sep='\t')
    data = data.dropna(subset=['groupID'], how='any')
    return data

def get_families_by_superfamily(superfamily_ids, data):
    orig_data = deepcopy(data)
    result = {}
    superfam_transcripts = {}
    for super_id in superfamily_ids:
        super_id = str(super_id)
        curated_data = data.loc[[super_id in str(x) for x in data['exon_associated_superfamilies']]]
        
        cols_of_interest = ['exon_associated_superfamilies', 'exon_associated_domains', 'exon_associated_domains_meaning']
        curated_data = curated_data[cols_of_interest+['transcript_id']]
        info = {}
        for col in cols_of_interest:
            info[col] = []
        
        superfam_transcripts[super_id] = []
        for i in range(curated_data.shape[0]): 
            superfam_transcripts[super_id].append(curated_data.loc[curated_data.index[i], 'transcript_id'])
            for col in cols_of_interest:
                #print(curated_data.loc[curated_data.index[i], col])
                info[col] = info[col] + str(curated_data.loc[curated_data.index[i], col]).split(';')

        
        result[super_id] = []
        for superfam_id, fam_id, fam_desc in zip(
                                                info['exon_associated_superfamilies'], 
                                                info['exon_associated_domains'], 
                                                info['exon_associated_domains_meaning']):
            value = (fam_id, fam_desc)
            if superfam_id == super_id and value not in result[super_id]:
                result[super_id].append(value)
    return result, superfam_transcripts

if __name__ == '__main__':
    data = load_data()
    result, transcripts = get_families_by_superfamily(['49265', '48065', '109885', '48726', '53067'], data)
    
    for key, value in result.items():
        print(f"{key}: {value}")
    for key in transcripts:
        print(f"{key}: {np.unique(transcripts[key])}")


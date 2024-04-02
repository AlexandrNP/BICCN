import os
import pandas as pd

def get_coord_data(domains_data):
    data_list = []
    for i in range(domains_data.shape[0]):
        residues = domains_data.loc[i, 'residues']
        
        exon = domains_data.loc[i, 'exon']
        exon_id = exon
        tx_id = domains_data.loc[i, 'transcript_id']

        if residues == '-':
            data_list.append( [exon, exon_id, tx_id, -1, -1] )
            continue

        for residue_pair in residues.split(';'):
            start, end = residue_pair.split('-')
            if exon == 'chr10_106903643_106903705_+':
                print([exon, exon_id, tx_id, start, end])
            data_list.append( [exon, exon_id, tx_id, int(start), int(end)] )

    coord_data = pd.DataFrame(data_list, columns=['exon', 'exon_id', 'tx_id', 'start', 'end'])
    coord_data = coord_data.groupby(['exon', 'exon_id', 'tx_id'], as_index=False).agg({'start': min, 'end': max})
    return coord_data


if __name__ == "__main__":
    domain_dir = os.path.join('Data', 'Domains')
    exon_dir = os.path.join('Data', 'Exons')
    domain_data_file = 'domains_human.tsv'
    domains_data = pd.read_csv( os.path.join(domain_dir, domain_data_file), sep='\t' )
    coords = get_coord_data(domains_data)
    coords.to_csv( os.path.join(exon_dir, 'protein_seq_for_exons_coords_human.tsv'), sep='\t', index=None )
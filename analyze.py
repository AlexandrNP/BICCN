from xml import dom
import matplotlib.pyplot as plt
import pingouin as pg
import pandas as pd
import numpy as np
import os
from copy import deepcopy
import gseapy
from statsmodels.stats.multitest import multipletests

DATA_DIR = 'Data'
EXON_DIR = os.path.join(DATA_DIR, 'Exons')
DOMAIN_DIR = os.path.join(DATA_DIR, 'Domains')
RANKED_SET_DIR = 'ranked_sets'
GSEA_DIR = 'GSEA'


def apply_appris_filter(domain_data):
    appris_transcripts = get_appris_transcripts()
    # print(protein_seq_dataset)
    domain_data_size = np.shape(domain_data)[0]
    protein_seq_dataset = pd.read_csv(os.path.join(
        EXON_DIR, 'protein_seq_for_exons_var.tsv'), sep='\t')
    if domain_data_size > 100000:
        protein_seq_dataset = pd.read_csv(os.path.join(
            EXON_DIR, 'protein_seq_for_exons_null.tsv'), sep='\t')
    protein_seq_dataset.columns = [x.lower()
                                   for x in protein_seq_dataset.columns]
    print(f'Number of transcripts in APPRIS:\t{len(appris_transcripts)}')
    print(f'Domain dataset length before filtering:\t{domain_data_size}')
    print(
        f'Protein seq dataset before filtering:\t{protein_seq_dataset.shape[0]}')

    prot_counts = domain_data.groupby(['exon', 'protein_id'])[
        'protein_id'].count()
    # print(prot_counts)
    mean_prot_num = np.mean(prot_counts)
    prot_dev = np.std(prot_counts)
    print(
        f'Mean number of associated proteins before filtering:\t{mean_prot_num}+/-{prot_dev}')
    print(
        f'Percentage of exons with multiple protein associations:\t{float(sum(prot_counts > 1))/len(prot_counts)}')

    protein_seq_dataset_curated = protein_seq_dataset[~protein_seq_dataset['protein_sequence'].isin([
                                                                                                    "-"])]
    prot_seq_counts = protein_seq_dataset_curated.groupby(
        ['exon', 'protein_sequence'])['exon'].count()
    # print(prot_seq_counts)
    mean_prot_seq_num = np.mean(prot_seq_counts)
    prot_seq_dev = np.std(prot_seq_counts)
    print(
        f'Mean number of associated unique proteins sequences before filtering:\t{mean_prot_seq_num}+/-{prot_seq_dev}')
    print(
        f'Percentage of exons with multiple protein sequence associations:\t{float(sum(prot_seq_counts > 1))/len(prot_seq_counts)}')

    N = sum(domain_data['transcript_id'].isin(appris_transcripts))
    print(f'Number of APPRIS transcripts in domain dataset:\t{N}')
    domain_data = domain_data[domain_data['transcript_id'].isin(
        appris_transcripts)]
    print(
        f'Domain dataset length after APPRIS filtering:\t{np.shape(domain_data)[0]}')
    protein_seq_dataset = protein_seq_dataset_curated[protein_seq_dataset_curated['tx_id'].isin(
        appris_transcripts)]
    print(
        f'Protein seq dataset after APPRIS filtering:\t{protein_seq_dataset.shape[0]}')

    prot_counts = domain_data.groupby(['exon', 'protein_id'])[
        'protein_id'].count()
    # print(prot_counts)
    mean_prot_num = np.mean(prot_counts)
    prot_dev = np.std(prot_counts)
    print(
        f'Mean number of associated proteins after APPRIS filtering:\t{mean_prot_num}+/-{prot_dev}')
    print(
        f'Percentage of exons with multiple protein sequence associations:\t{float(sum(prot_counts > 1))/len(prot_counts)}')

    prot_seq_counts = protein_seq_dataset.groupby(
        ['exon', 'protein_sequence'])['exon'].count()
    # print(prot_seq_counts)
    mean_prot_seq_num = np.mean(prot_seq_counts)
    prot_seq_dev = np.std(prot_seq_counts)
    print(
        f'Mean number of associated unique proteins sequences after APPRIS filtering:\t{mean_prot_seq_num}+/-{prot_seq_dev}')
    print(
        f'Percentage of exons with multiple protein sequence associations:\t{float(sum(prot_seq_counts > 1))/len(prot_seq_counts)}')

    return domain_data


def get_var_domains_data(var_domains_filename = 'domains_protein_seq_for_exons_var.tsv',
                        coords_for_exons_filename = 'protein_seq_for_exons_coords_variable.tsv',
                    appris_filter=True):
    import re
    #var_domains_filename = 'domains_protein_seq_for_exons_var.tsv'
    domain_data = pd.read_csv(os.path.join(
        DOMAIN_DIR, var_domains_filename), sep='\t', index_col=None, dtype=str, low_memory=False)
    if appris_filter:
        domain_data = apply_appris_filter(domain_data)
    coords_path = os.path.join(EXON_DIR, coords_for_exons_filename)
    exon_coords = pd.read_csv(coords_path, sep='\t', index_col=None)
    domain_data = pd.merge(domain_data, exon_coords, how='left', left_on=['exon', 'transcript_id'], right_on=['exon', 'tx_id'])
    domain_data['exon_associated_domains'] = np.repeat('-', repeats=domain_data.shape[0])
    domain_data['exon_associated_families'] = np.repeat('-', repeats=domain_data.shape[0])
    domain_data['exon_associated_superfamilies'] = np.repeat('-', repeats=domain_data.shape[0])
    associated_exons_num = 0
    for i in range(domain_data.shape[0]):
        associated_domains = []
        associated_families = []
        associated_superfamilies = []
        start = domain_data.loc[domain_data.index[i], 'start']
        end = domain_data.loc[domain_data.index[i], 'end']
        domains = domain_data.loc[domain_data.index[i], 'domain_name'].split(';')
        families = domain_data.loc[domain_data.index[i], 'domains_family'].split(';')
        superfamilies = domain_data.loc[domain_data.index[i], 'domains_superfamily'].split(';')
        residues = domain_data.loc[domain_data.index[i], 'residues'].split(';')
        for coords, domain, family, superfamily in zip(residues, domains, families, superfamilies):
            for coord in coords.split(','):
                #print(coord)
                if coord == '-':
                    continue
                left, right = [int(x) for x in coord.split('-')]
                if (start > left and end <= right) or (start <= left and end > left) or (start > left and start <= right) or (start < left and end >= right):
                    associated_domains.append(domain)    
                    associated_families.append(family)
                    associated_superfamilies.append(superfamily)
                    associated_exons_num += 1
                if len(associated_domains) > 0:
                    domain_data.loc[domain_data.index[i], 'exon_associated_domains'] = ';'.join(associated_families)
                    domain_data.loc[domain_data.index[i], 'exon_associated_families'] = ';'.join(associated_domains)
                    domain_data.loc[domain_data.index[i], 'exon_associated_superfamilies'] = ';'.join(associated_superfamilies)
        print(f'Associated exons num: {associated_exons_num}')
        print(f'Associated exons percentage: {associated_exons_num/domain_data.shape[0]}')


    print(domain_data)
    return domain_data


def print_col_types(df):
    for col in df:
        col_types = []
        for x in df[col].apply(type):
            col_types.append(str(x))
        print(np.unique(col_types))


def get_appris_transcripts():
    appris_db_file = 'appris_data.principal.txt'
    appris_data = pd.read_csv(os.path.join(
        EXON_DIR, appris_db_file), sep='\t', index_col=None, header=None)
    return appris_data[appris_data.columns[2]]


def get_human_dataset():
    human_data = pd.read_csv(os.path.join(
        EXON_DIR, 'humanHighLowVar_mousePSI_mm10_wGeneNames'), sep='\t',  dtype=str, low_memory=False)
    return human_data

def get_exon_groupings():
    eve_data = pd.read_csv(os.path.join(
        EXON_DIR, 'EVEx_4axes'), sep='\t',  dtype=str, low_memory=False)
    developmental_data = pd.read_csv(os.path.join(
        EXON_DIR, 'exonList_developmentalModalities'), sep='\t',  dtype=str, low_memory=False)
    grouping = pd.read_csv(os.path.join(
        EXON_DIR, 'highlyVariableExons_groupedByVariability'), sep='\t',  dtype=str, low_memory=False)

    low_var_null_data = pd.read_csv(os.path.join(
        EXON_DIR, 'lowVariabilityAltExons_nullSet'), sep='\t',  dtype=str, low_memory=False)

    datasets = [eve_data, developmental_data, grouping]
    for i in range(len(datasets)):
        datasets[i] = pd.concat([datasets[i], low_var_null_data])
        datasets[i] = datasets[i].fillna('LowVarNull')

    merged = None
    # for ds in datasets:
    #    print(len(np.unique(ds['Exon'])))
    #    print(np.shape(ds))
    for i in range(len(datasets)-1):
        if 'Gene' in datasets[i]:
            datasets[i] = datasets[i].drop(['Gene'], axis=1)
        if 'Gene' in datasets[i+1]:
            datasets[i+1] = datasets[i+1].drop(['Gene'], axis=1)
        if merged is None:
            merged = pd.merge(
                datasets[i], datasets[i+1], on='Exon', how='outer')
        else:
            merged = pd.merge(merged, datasets[i+1], on='Exon', how='outer')

    return merged


def add_aggregated_value(data, factor, target_col, function, out_column_name):
    reduced = data.groupby(factor).agg({target_col: [function]})
    new_dataframe = pd.DataFrame()
    new_dataframe[factor] = reduced.index.values
    new_dataframe[out_column_name] = reduced[target_col].values
    merged = pd.merge(data, new_dataframe, on=factor, how='left')
    return merged


def get_all_domains(domain_ids):
    all_domains = []
    for domains in domain_ids:
        all_domains = all_domains + domains.split(';')
    return all_domains


def add_unique_domains_num(data):
    def unique_associated_domain_number(domain_ids):
        all_domains = get_all_domains(domain_ids)
        return len(np.unique(all_domains))

    merged = add_aggregated_value(
        data, 'exon', 'domain_name', unique_associated_domain_number, 'unique_domains_num')
    return merged


def add_all_unique_domains(data):
    def add_unique_domains(domains_ids):
        all_domains = np.unique(get_all_domains(domains_ids))
        return ';'.join(all_domains)

    merged = add_aggregated_value(
        data, 'exon', 'domain_name', add_unique_domains, 'unique_domains')
    return merged


def add_max_number_of_domains_per_exon(data):
    def max_associated_domains(residues):
        max_domains_num = 0
        for exon_residues in residues:
            residues_per_domain = [x for x in exon_residues.split(';')]
            all_domains = []
            for domain_residues in residues_per_domain:
                all_domains = all_domains + domain_residues.split(',')
            filtered_domain_num = 0

            for domain_boundaries in all_domains:
                if domain_boundaries == '' or domain_boundaries == '-':
                    continue
                left, right = domain_boundaries.split('-')
                #print(int(right), int(left))
                if int(right) - int(left) > 10:
                    filtered_domain_num += 1
                    # print(filtered_domain_num)
            if filtered_domain_num > max_domains_num:
                max_domains_num = filtered_domain_num
                # if max_domains_num == 252:
                #    print(all_domains)
                #    print(max_domains_num)
                #    print(len(residues))
        return max_domains_num

    merged = add_aggregated_value(
        data, 'exon', 'residues', max_associated_domains, 'max_domains_num')
    return merged


def add_min_number_of_domains_per_exon(data):
    def min_associated_domains(residues):
        min_domains_num = 1000
        for exon_residues in residues:
            residues_per_domain = [x for x in exon_residues.split(';')]
            all_domains = []
            for domain_residues in residues_per_domain:
                all_domains = all_domains + domain_residues.split(',')
            filtered_domain_num = 0
            for domain_boundaries in all_domains:
                if domain_boundaries == '' or domain_boundaries == '-':
                    continue
                left, right = domain_boundaries.split('-')
                if int(right) - int(left) > 10:
                    filtered_domain_num += 1
            if filtered_domain_num < min_domains_num:
                min_domains_num = filtered_domain_num
        return min_domains_num

    merged = add_aggregated_value(
        data, 'Exon', 'residues', min_associated_domains, 'min_domains_num')
    return merged


def add_mean_number_of_domains_per_exon(data):
    def mean_associated_domains(residues):
        all_domains_num = 0
        for exon_residues in residues:
            residues_per_domain = [x for x in exon_residues.split(';')]
            all_domains = []
            for domain_residues in residues_per_domain:
                all_domains = all_domains + domain_residues.split(',')

            for domain_boundaries in all_domains:
                if domain_boundaries == '' or domain_boundaries == '-':
                    continue
                left, right = domain_boundaries.split('-')
                if int(right) - int(left) > 10:
                    all_domains_num += 1
        n = len(residues)
        n = n if n > 0 else 1
        return float(all_domains_num)/n

    merged = add_aggregated_value(
        data, 'Exon', 'residues', mean_associated_domains, 'mean_domains_num')
    return merged

def add_n_terminus_association(data):
    n = data.shape[0]
    data['n_terminus'] = np.repeat('', n)
    for i in range(n):
        domain_start = data.loc[data.index[i], 'residues'].split(';')[0].split('-')[0]
        start = str(data.loc[data.index[i], 'start'])
        if len(start) > 0:
            if len(domain_start) < 1 or float(start) <= float(domain_start):
                data.loc[data.index[i], 'n_terminus'] = '+'
        
   
    return data

def add_statistic(dataframe):
    expanded_dataframe = add_unique_domains_num(dataframe)
    expanded_dataframe = add_max_number_of_domains_per_exon(expanded_dataframe)
    expanded_dataframe = add_min_number_of_domains_per_exon(expanded_dataframe)
    expanded_dataframe = add_mean_number_of_domains_per_exon(
        expanded_dataframe)
        
    expanded_dataframe = add_all_unique_domains(expanded_dataframe)
    expanded_dataframe['domain_num_discrepancy'] = expanded_dataframe['max_domains_num'] - \
        expanded_dataframe['min_domains_num']
    expanded_dataframe = add_n_terminus_association(expanded_dataframe)
    return expanded_dataframe


def process_numerical_statistics(data, group, variables):
    res = data.groupby([group, 'exon']).first().reset_index()
    for variable in variables:
        # print(res)
        res.boxplot(variable, by=group)
        ax = plt.gca()
        ax.set_ylim([-0.5, 18])
        plt.savefig(f'results/{group}_{variable}_boxplot.pdf')
        anova_res = pg.anova(data=res, dv=variable,
                             between=group, detailed=True)
        
        #pt = pg.pairwise_tukey(data=res, dv=variable, between=group)
        #anova_res.to_csv(f'results/{group}_{variable}_anova.tsv', sep='\t')
        #pt.to_csv(f'results/{group}_{variable}_tukey_hds.tsv', sep='\t')


def process_distributions_statistics(data, group, variables):
    res = data.groupby([group, 'exon']).first().reset_index()

    plot_kinds = ['bar', 'density']
    for plot_kind in plot_kinds:
        for variable in variables:
            common_occurences = {}
            var_res = deepcopy(res)
            var_res[variable] = [x.split(';') for x in var_res[variable]]
            var_res = var_res.explode(variable, ignore_index=True)
            var_res = var_res.fillna('-')
            all_groups = np.unique(var_res[group])

            n = len(all_groups)
            dim_x_size = int(np.ceil(np.sqrt(n)))
            dim_y_size = int(np.ceil(float(n) / dim_x_size))
            print(dim_x_size, dim_y_size, n)

            for current_group in all_groups:
                group_idx = np.array(
                    [x == current_group for x in var_res[variable]])
                group_data = var_res[variable][group_idx]

            domains_numbers = len(np.unique(var_res[variable]))
            print(var_res[variable].astype(str).values)
            series = pd.Series(var_res[variable].astype(
                str).values, index=var_res[group])
            top_10_unified = []
            for idx in all_groups:
                hist_data = series[idx].value_counts(sort=True, ascending=False)
                top_10_unified = top_10_unified + list(hist_data[:5].index)
            top_10_unified = np.unique(top_10_unified)

            # Regular plotting
            f, a = plt.subplots(dim_x_size, dim_y_size, figsize=(48, 20))
            a = np.array(a)
            f.tight_layout(pad=20.0)
            f.suptitle(f'{group} - {variable}', fontsize=24)
            a = a.ravel()
            i = 0
            for idx in all_groups:
                plt.axes(a[i])
                a[i].set_title(idx)
                a[i].tick_params(axis='both', which='major', labelsize=7)
                hist_data = series[idx].value_counts(sort=True, ascending=False)

                for entry in top_10_unified:
                    if entry not in hist_data.index:
                        hist_data[entry] = 0
                hist_data = hist_data.sort_values(ascending=False)
                hist_data = hist_data.iloc[np.lexsort([hist_data.index, hist_data.values])][::-1]
                top_10 = hist_data[top_10_unified]
                top_10 = top_10.iloc[np.lexsort([top_10.index, top_10.values])][::-1]

                pd.DataFrame(top_10).to_csv(
                    f'ranked_sets/top10_united_hist_data_{group}_{variable}_{idx}.tsv', sep='\t')
                top_10_unified = top_10_unified[top_10_unified != '-']
                pd.DataFrame(hist_data).to_csv(
                    f'ranked_sets/ranked_hist_data_{group}_{variable}_{idx}.tsv', sep='\t')
                hist_data = hist_data[hist_data.index != '-']

                if 'Null' not in idx:
                    if plot_kind == 'bar':
                        if 'Group' not in idx and len(idx) > 1:
                            #a[i].set_ylim(0, 10)
                            a[i].set_ylim(0, 150)
                        if 'Group' in idx:
                            a[i].set_ylim(0, 14)
                        if 'Hippocampus' in idx or 'VisCortex' in idx:
                            a[i].set_ylim(0, 16)
                        if len(idx) < 2:
                            a[i].set_ylim(0, 8)
                    else:
                        if 'Group' not in idx and len(idx) > 1:
                            a[i].set_ylim(0, 0.5)
                        if 'Group' in idx:
                            a[i].set_ylim(0, 0.18)
                        if 'Hippocampus' in idx or 'VisCortex' in idx:
                            a[i].set_ylim(0, 0.05)
                        if len(idx) < 2:
                            a[i].set_ylim(0, 0.15)
                

                if plot_kind == 'density':
                    hist_data = hist_data / hist_data.sum()

                hist_data.plot(kind='bar')
                hist_data.to_csv(f'results_data/{group}_{variable}_{plot_kind}_hist.tsv', sep='\t')
                i += 1

            for i in np.arange(i, dim_x_size*dim_y_size):
                f.delaxes(a[i])

            plt.savefig(f'results/{group}_{variable}_{plot_kind}_hist.pdf')

            # Top 10 united plotting
            f, a = plt.subplots(dim_x_size, dim_y_size, figsize=(48, 20))
            a = np.array(a)
            f.tight_layout(pad=20.0)
            f.suptitle(f'{group} - {variable}', fontsize=24)
            a = a.ravel()
            i = 0
            for idx in all_groups:
                plt.axes(a[i])
                a[i].set_title(idx)
                a[i].tick_params(axis='both', which='major', labelsize=7)
                hist_data = series[idx].value_counts(sort=True, ascending=False)
                
                for uni_idx in top_10_unified:
                    if uni_idx not in hist_data.index:
                        hist_data[uni_idx] = 0
                # .sort_values(ascending=False)
                if plot_kind == 'density':
                    hist_data = hist_data[hist_data.index != '-']
                    hist_data = hist_data / hist_data.sum()
                    
                hist_data = hist_data[top_10_unified]
                common_occurences[idx] = {}
                #print('HIST DATA')
                #print(type(hist_data))
                hist_data = hist_data[hist_data.index != '-']
                

                if 'Null' not in idx:
                    if plot_kind == 'bar':
                        if 'Group' not in idx and len(idx) > 1:
                            #a[i].set_ylim(0, 10)
                            a[i].set_ylim(0, 150)
                        if 'Group' in idx:
                            a[i].set_ylim(0, 14)
                        if 'Hippocampus' in idx or 'VisCortex' in idx:
                            a[i].set_ylim(0, 16)
                        if len(idx) < 2:
                            a[i].set_ylim(0, 8)
                    else:
                        if 'Group' not in idx and len(idx) > 1:
                            a[i].set_ylim(0, 0.5)
                        if 'Group' in idx:
                            a[i].set_ylim(0, 0.18)
                        if 'Hippocampus' in idx or 'VisCortex' in idx:
                            a[i].set_ylim(0, 0.05)
                        if len(idx) < 2:
                            a[i].set_ylim(0, 0.15)

                hist_data.plot(kind='bar')
                hist_data.to_csv(f'results_data/{group}_{variable}_{plot_kind}_top10_united_hist.tsv', sep='\t')
                i += 1

            for i in np.arange(i, dim_x_size*dim_y_size):
                f.delaxes(a[i])

            plt.savefig(f'results/{group}_{variable}_{plot_kind}_top10_united_hist.pdf')
            # for idx in common_occurences:


def add_null_set(dataframe, appris_filter=True):
    # , 'domains_protein_seq_for_exons_null_remaining.tsv']
    null_set_files = ['domains_protein_seq_for_exons_null_partial.tsv']
    for null_set_file in null_set_files:
        null_data = pd.read_csv(os.path.join(
            DOMAIN_DIR, null_set_file), sep='\t', index_col=0)
        if appris_filter:
            null_data = apply_appris_filter(null_data)

        coords_path = os.path.join(EXON_DIR, 'protein_seq_for_exons_coords_conserved.tsv')
        exon_coords = pd.read_csv(coords_path, sep='\t', index_col=None)
        domain_data = pd.merge(null_data, exon_coords, how='left', left_on=['exon', 'transcript_id'], right_on=['exon', 'tx_id'])
        domain_data['exon_associated_domains'] = np.repeat('-', repeats=domain_data.shape[0])
        domain_data['exon_associated_families'] = np.repeat('-', repeats=domain_data.shape[0])
        domain_data['exon_associated_superfamilies'] = np.repeat('-', repeats=domain_data.shape[0])
        associated_exons_num = 0
        for i in range(domain_data.shape[0]):
            associated_domains = []
            associated_families = []
            associated_superfamilies = []
            start = domain_data.loc[domain_data.index[i], 'start']
            end = domain_data.loc[domain_data.index[i], 'end']
            domains = domain_data.loc[domain_data.index[i], 'domain_name'].split(';')
            families = domain_data.loc[domain_data.index[i], 'domains_family'].split(';')
            superfamilies = domain_data.loc[domain_data.index[i], 'domains_superfamily'].split(';')
            residues = domain_data.loc[domain_data.index[i], 'residues'].split(';')
            for coords, domain, family, superfamily in zip(residues, domains, families, superfamilies):
                for coord in coords.split(','):
                    #print(coord)
                    if coord == '-':
                        continue
                    left, right = [int(x) for x in coord.split('-')]
                    if (start > left and end <= right) or (start <= left and end > left) or (start > left and start <= right) or (start < left and end >= right):
                        associated_domains.append(domain)    
                        associated_families.append(family)
                        associated_superfamilies.append(superfamily)
                        associated_exons_num += 1
                    if len(associated_domains) > 0:
                        domain_data.loc[domain_data.index[i], 'exon_associated_domains'] = ';'.join(np.unique(associated_families))
                        domain_data.loc[domain_data.index[i], 'exon_associated_families'] = ';'.join(np.unique(associated_domains))
                        domain_data.loc[domain_data.index[i], 'exon_associated_superfamilies'] = ';'.join(np.unique(associated_superfamilies))
            #print(f'Associated exons num: {associated_exons_num}')
            #print(f'Associated exons percentage: {associated_exons_num/domain_data.shape[0]}')


        dataframe = pd.concat([dataframe, domain_data])
        dataframe = dataframe.reset_index()
    dataframe = dataframe.fillna('')
    return dataframe


def load_data(expanded_data_file, skip_var = False):
    data = None
    if not os.path.isfile(expanded_data_file):
        print('Creating new file...')
        exon_data = get_human_dataset() #get_exon_groupings()
        exon_data['exon'] = exon_data['Exon']
        var_data = get_var_domains_data(var_domains_filename='domains_human.tsv',
                                        coords_for_exons_filename='human_coords_exon.tsv') #protein_seq_for_exons_coords_human.tsv
        print(f'Exon dataset shape:\t{np.shape(exon_data)}')
        if not skip_var:
            var_domain_data = get_var_domains_data(appris_filter=True)
            print(f'Domain dataset shape:\t{np.shape(var_domain_data)}')

            var_data = pd.merge(exon_data, var_domain_data,
                            left_on='Exon', right_on='exon', how='left')

            merge_key = ["Exon", "Cluster", "Region", "groupID", "Group"]
            var_data = var_data.groupby(by=merge_key, as_index=False, dropna=False).first()

            var_data = add_null_set(var_data, appris_filter=True)

            var_data = var_data.fillna('')
            var_data = pd.DataFrame(var_data, dtype='str')

        var_data['Exon'] = var_data['exon']
        expanded_var_data = add_statistic(var_data)
            #expanded_var_data = expanded_var_data.groupby(['exon'], sort=False)[
            #    'max_domains_num'].max()

        
        domain_mapping = {}
        

        
        scope_data = pd.read_csv('dir.cla.scope.2.06-stable.txt', sep='\t')
        for i in range(scope_data.shape[0]):
            domain = scope_data.loc[i, scope_data.columns[0]]
            classification = scope_data.loc[i, scope_data.columns[-1]].split(',')
            cl = classification[0].split('=')[-1]
            cf = classification[1].split('=')[-1]
            domain_mapping[domain] = {}
            domain_mapping[domain]['cl'] = cl
            domain_mapping[domain]['cf'] = cf


        scope_des = pd.read_csv('dir.des.scope.2.06-stable.txt', sep='\t', header=None)
        biology_map = {}
        test_file = open('map_test.txt', 'w')
        for i in range(scope_des.shape[0]):
            key1 = str(scope_des.loc[i, scope_des.columns[0]])
            key2 = str(scope_des.loc[i, scope_des.columns[-2]])
            value = scope_des.loc[i, scope_des.columns[-1]]
            biology_map[key1] = value
            biology_map[key2] = value
            test_file.write(key1)
            test_file.write('\n')
        test_file.close()

        default_domains = ['', '-']
        for domain in default_domains:
            domain_mapping[domain] = {}
            domain_mapping[domain]['cl'] = '-'
            domain_mapping[domain]['cf'] = '-'
            biology_map[domain] = domain
       
        cl_columns = ['cl', 'cf', 'cl_associated', 'cf_associated'] 
        biological_columns_map = {'cl': 'cl_meaning', 
                        'cf': 'cf_meaning', 
                        'cl_associated': 'cl_associated_meaning', 
                        'cf_associated': 'cf_associated_meaning',
                        'domain_name': 'domain_meaning', 
                        'domains_family': 'family_meaning', 
                        'domains_superfamily': 'superfamily_meaning',
                        'exon_associated_domains': 'exon_associated_domains_meaning',
                        'exon_associated_families': 'exon_associated_domains_families_meaning',
                        'exon_associated_superfamilies': 'exon_associated_domains_superfamilies_meaning'}
        
        new_columns = cl_columns + list(biological_columns_map.values())
        for col in new_columns:
            expanded_var_data[col] = np.repeat('', expanded_var_data.shape[0])
        

        expanded_var_data = expanded_var_data.reindex()
        index = expanded_var_data.index
        for i in range(expanded_var_data.shape[0]):
            domain_name = expanded_var_data.loc[index[i], 'domain_name']
            exon_associated_domains = expanded_var_data.loc[index[i], 'exon_associated_families']
            expanded_var_data.loc[index[i], 'cl'] = ';'.join([domain_mapping[x]['cl'] for x in domain_name.split(';')])
            expanded_var_data.loc[index[i], 'cf'] = ';'.join([domain_mapping[x]['cf'] for x in domain_name.split(';')])
            expanded_var_data.loc[index[i], 'cl_associated'] = ';'.join([domain_mapping[x]['cl'] for x in exon_associated_domains.split(';')])
            expanded_var_data.loc[index[i], 'cf_associated'] = ';'.join([domain_mapping[x]['cf'] for x in exon_associated_domains.split(';')])

            # Append biological meaning
            for map_from_col, map_to_col in biological_columns_map.items():
                elements = str(expanded_var_data.loc[index[i], map_from_col]).split(';')
                expanded_var_data.loc[index[i], map_to_col] = ';'.join(biology_map[str(x)] for x in elements)

        data = expanded_var_data
        data.to_csv(expanded_data_file, sep='\t', index=None)



    else:
        data = pd.read_csv(expanded_data_file, sep='\t')
        #print(data['Exon'].)
        #exons_mask1 = [len(str(x)) < 2 or str(x) == 'NaN' for x in data['Exon']]
        #exons_mask2 = [len(str(x)) < 2 or str(x) == 'NaN' for x in data['exon']]

        #print(sum(exons_mask1))

    return data


def run_gsea_preranked(ranked_set_dir, group, hist_variables, set_size=10):
    file_names_ranked = {}
    file_names_top10 = {}
    for hist_variable in hist_variables:
        for filename in os.listdir(ranked_set_dir):
            if 'hist_data' not in filename or \
                    group not in filename or \
                    hist_variable not in filename:
                continue
            # print(hist_variable)
            file_path = os.path.join(ranked_set_dir, filename)
            variable = filename.split('.')[0].split('_')[-1]
            gene_set_data = pd.read_csv(file_path, sep='\t',
                                        index_col=0)
            list_id = f'{variable}_{hist_variable}'
            if 'ranked' in filename:
                ranked_lists[list_id] = gene_set_data
                file_names_ranked[list_id] = filename
            if 'top10' in filename:
                file_names_top10[list_id] = filename
                gene_set_data = gene_set_data.sort_values(
                    by=gene_set_data.columns[0], ascending=False)
                sets[list_id] = gene_set_data.index
                scores[list_id] = gene_set_data.values

    for hist_variable in hist_variables:
        # print(hist_variable)
        for var1 in ranked_lists:
            if hist_variable not in var1:
                continue
            gene_set = {}
            for var2 in ranked_lists:
                if var1 == var2 or hist_variable not in var2:
                    continue


                gene_set[var2] = list(
                    ranked_lists[var2].index[:set_size].values)

                preranked = gseapy.prerank(rnk=ranked_lists[var1],
                                           gene_sets=gene_set,
                                           min_size=2,
                                           max_size=10000,
                                           outdir=os.path.join(
                    GSEA_DIR, f'{var1}-{var2}'),
                    verbose=False)
                preranked.res2d.head()
    set_size = 10
    file_names_ranked = {}
    file_names_top10 = {}
    for hist_variable in hist_variables:
        for filename in os.listdir(ranked_set_dir):
            if 'hist_data' not in filename or \
                    group not in filename or \
                    hist_variable not in filename:
                continue
            # print(hist_variable)
            file_path = os.path.join(ranked_set_dir, filename)
            variable = filename.split('.')[0].split('_')[-1]
            gene_set_data = pd.read_csv(file_path, sep='\t',
                                        index_col=0)
            list_id = f'{variable}_{hist_variable}'
            if 'ranked' in filename:
                ranked_lists[list_id] = gene_set_data
                file_names_ranked[list_id] = filename
            if 'top10' in filename:
                file_names_top10[list_id] = filename
                gene_set_data = gene_set_data.sort_values(
                    by=gene_set_data.columns[0], ascending=False)
                sets[list_id] = gene_set_data.index
                scores[list_id] = gene_set_data.values

    for hist_variable in hist_variables:
        # print(hist_variable)
        for var1 in ranked_lists:
            if hist_variable not in var1:
                continue
            for var2 in ranked_lists:
                if var1 == var2 or hist_variable not in var2:
                    continue

                gene_set = {}
                gene_set[var2] = list(
                    ranked_lists[var2].index[:set_size].values)

                preranked = gseapy.prerank(rnk=ranked_lists[var1],
                                           gene_sets=gene_set,
                                           min_size=3,
                                           max_size=10000,
                                           outdir=os.path.join(
                    GSEA_DIR, f'{var1}-{var2}'),
                    verbose=False)
                preranked.res2d.head()


def correct_gsea(gsea_dir):
    gsea_dataset_list = []
    for dirname in os.listdir(gsea_dir):
        if dirname in ['.', '..', '.DS_Store']:
            continue
        dirpath = os.path.join(gsea_dir, dirname)
        filename = os.path.join(dirpath, 'gseapy.gene_set.prerank.report.csv')
        data = pd.read_csv(filename, index_col=None)
        ranked_list, gene_set = dirname.split('-')
        data['RankedList'] = ranked_list
        data['GeneSet'] = gene_set
        gsea_dataset_list.append(data)
    gsea = pd.concat(gsea_dataset_list)
    gsea = gsea.reset_index()
    for i in range(len(gsea['NOM p-val'])):
        p_val = gsea.loc[gsea.index[i], 'NOM p-val']
        #print(p_val)
        #print(type(p_val))
        if p_val < 1E-20:
            gsea.loc[gsea.index[i], 'NOM p-val'] = 1E-3
    gsea = gsea.sort_values(by='NOM p-val', axis=0, ascending=True)
    corrected_pval = multipletests(
        gsea['NOM p-val'], method='fdr_bh', is_sorted=True)
    # print(multipletests)
    gsea['FDR_Benjamini-Hochberg'] = corrected_pval[1]
    gsea.to_csv('GSEA_all_sets.tsv', sep='\t', index=None)


if __name__ == "__main__":
    #expanded_data_file = os.path.join(DOMAIN_DIR, 'expanded_data.tsv')
    data = load_data(os.path.join(DOMAIN_DIR, 'human_data.tsv'), skip_var=True) #expanded_data_file
    data['Group'] = np.repeat('Human', data.shape[0])
    data['group'] = np.repeat('Human', data.shape[0])

    groups = ['group'] #['Cluster', 'Region', 'groupID', 'Group', ]
    num_variables = ['max_domains_num', 'min_domains_num',
                     'mean_domains_num', 'domain_num_discrepancy', 'unique_domains_num']
    #hist_variables = ['domains_family',
    #                  'domains_superfamily', 'domain_name', 'unique_domains']
    #hist_variables = ['exon_associated_superfamilies', 'exon_associated_families', 'exon_associated_domains']
    hist_variables = ['exon_associated_domains_superfamilies_meaning', 
                        'cf_associated_meaning', 
                        'exon_associated_domains_meaning', 
                        'exon_associated_domains_families_meaning']

    for group in groups:
        ranked_lists = {}
        sets = {}
        scores = {}

        # Hiding for Human dataset
        #data[group] = [x if x != '' and x != None and x != np.nan and x !=
        #               'nan' and type(x) == str else 'ConsNull' for x in data[group]]


        # print(np.unique(data[group]))
        # print(len(data[group]))

        process_numerical_statistics(data, group, num_variables)
        process_distributions_statistics(data, group, hist_variables)
        #run_gsea_preranked(RANKED_SET_DIR, group, hist_variables, set_size=10)
        #correct_gsea(GSEA_DIR)

    # [reduced.index.transpose(), reduced['domain_name'].values], columns=['Exon', 'unique_domains_num']
    #unique_exons = get_largest_number_of_domains(var_data)
    # print(unique_exons)
    # print(len(np.unique(var_data['Exon'])))

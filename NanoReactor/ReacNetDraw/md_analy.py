#!/usr/bin/env python3

import re
import ast

from collections import Counter, defaultdict
from operator import itemgetter

import pandas as pd

from tqdm import tqdm
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from _ioToolkit import getLocalFolderList, extendFolderList


atom_list = ['H', 'C', 'N', 'O', 'S']

folder_list = getLocalFolderList()
csv_list = extendFolderList(folder_list, '_ref_relation.csv')

replac_dict = {
    '[][H[]][][H[]][]': '[HH]', 
    '[][H[]-[]][]': '[H]'
    }


def alpha_count(string, alpha=None):
    if alpha != None:
        return string.count(alpha)
    else:
        return len(re.findall('[a-zA-Z]', string))


def getReacRelat():
    df_list = [pd.read_csv(csv) for csv in tqdm(csv_list, desc='Analysing reaction data')]
    df_concat = pd.concat(df_list)

    reac_list = [i for i in df_concat['reaction'].to_list()]
    reac_traj_list = [idx for idx, df in enumerate(df_list) 
        if df.empty == False]

    ### build a list containing all realation ###
    relats = [ast.literal_eval(i) 
        for i in df_concat['relation'].to_list()]
    relat_list = [(replac_dict.get(i[0], i[0]), replac_dict.get(i[1], i[1]))
        for relat in relats 
            for i in relat 
                if replac_dict.get(i[0], i[0]) \
                    != replac_dict.get(i[1], i[1])]
    
    return relat_list, reac_list, reac_traj_list


def filterRare(l, freq_criter=2):
    ### get main reaction list ###
    freq = Counter(l)
    main_list = list(filter(
        lambda x: freq[x] > freq_criter, l))
    return main_list

top_left = '┌'
top_right = '┐'
bottom_left = '└'
bottom_right = '┘'
horizontal = '─'
vertical = '│'


def printReacInfo(reac_list, main_reac_list, reac_traj_list):
    reac_n = len(reac_list)
    reac_type_n = len(set(reac_list))
    main_reac_n = len(main_reac_list)
    main_reac_type_n = len(set(main_reac_list))
    reac_traj_n = len(reac_traj_list)
    
    content_width = 40
    label_width = 30
    number_width = content_width - label_width
    print('\n * Reactions Analysis Report')
    print(top_left + horizontal * content_width + top_right)
    print(f'{vertical}{" Reactive trajectory number:":<{label_width}}{reac_traj_n:>{number_width}}{vertical}')
    print(f'{vertical}{" Reaction number:":<{label_width}}{reac_n:>{number_width}}{vertical}')
    print(f'{vertical}{" Reaction class number:":<{label_width}}{reac_type_n:>{number_width}}{vertical}')
    print(f'{vertical}{" Main Reaction number:":<{label_width}}{main_reac_n:>{number_width}}{vertical}')
    print(f'{vertical}{" Main Reaction class number:":<{label_width}}{main_reac_type_n:>{number_width}}{vertical}')
    print(bottom_left + horizontal * content_width + bottom_right)
    print('\n')
    with open('reac_main_eve.log', 'w') as f:
        for i, count in Counter(main_reac_list).items():
            f.write(f'{i}\n')

def sortReacRelat(relat_list):
    ### sort the relation list by its freq and atom nums ###
    freq = Counter(relat_list)
    sort_relat_list = sorted(relat_list, key=lambda x: 
        (freq[x], -len(re.findall("[a-zA-Z]", x[0] + x[1]))), 
            reverse=True)
    return sort_relat_list


def clearRelat(sort_relat_list):
    ### make sure a way to node
    judge = True
    while judge:
        old = len(sort_relat_list)
        prd_set = set(list(map(itemgetter(1), sort_relat_list)))
        for relat in sort_relat_list:
            if relat[0] not in prd_set:
                sort_relat_list.remove(relat)
        new = len(sort_relat_list)
        if new == old:
            judge = False
    return sort_relat_list


def filterRelatList_C(relat_list, NC_range=None):
    ### remove relation of large C molecules
    smi_list = []
    filt_rela_list = []
    for relat in relat_list:
        if NC_range is None or \
            (NC_range[0] <= alpha_count(relat[0], 'C') <= NC_range[1] \
                and NC_range[0] <= alpha_count(relat[1], 'C') <= NC_range[1]):
            smi_list.append(relat[0])
            smi_list.append(relat[1])
            filt_rela_list.append(relat)
    return smi_list, filt_rela_list


def filterRelatList_atom(relat_list, Natom_range=None):
    ### remove relation of large molecules
    smi_list = []
    filt_rela_list = []
    for relat in relat_list:
        if Natom_range is None or \
            (Natom_range[0] <= len(re.findall("[a-zA-Z]", relat[0])) <= Natom_range[1]\
                and Natom_range[0] <= len(re.findall("[a-zA-Z]", relat[1])) <= Natom_range[1]):
            smi_list.append(relat[0])
            smi_list.append(relat[1])
            filt_rela_list.append(relat)
    return smi_list, filt_rela_list

def filterSMIList_atom(relat_list, Natom_range=None):
    ### remove relation of large molecules
    smi_list = []
    for relat in relat_list:
        if Natom_range == None:
            break
        if Natom_range[0] <= len(re.findall("[a-zA-Z]", relat[0])) <= Natom_range[1]:
            smi_list.append(relat[0])
        if Natom_range[0] <= len(re.findall("[a-zA-Z]", relat[1])) <= Natom_range[1]:
            smi_list.append(relat[1])
    
    ### sort smi list by its freq and atom nums ###
    freq = Counter(smi_list)
    sort_smi_list = sorted(
        smi_list, key=lambda x: 
        (freq[x], -len(re.findall("[a-zA-Z]", x))), 
        reverse=True)
        
    ### build a smi and idx dict based on the smi list ###
    ### key: smi, value: idx                           ###
    smi_dict = defaultdict(lambda: len(smi_dict) + 1)
    for smi in sort_smi_list:
        smi_dict[smi]
    return smi_dict

def getIdxList(smi_list, rela_list):
    ### sort smi list by its freq and atom nums ###
    freq = Counter(smi_list)
    sort_smi_list = sorted(
        smi_list, key=lambda x: 
        (freq[x], -len(re.findall("[a-zA-Z]", x))), 
        reverse=True)
        
    ### build a smi and idx dict based on the smi list ###
    ### key: smi, value: idx                           ###
    smi_dict = defaultdict(lambda: len(smi_dict) + 1)
    for smi in sort_smi_list:
        smi_dict[smi]
    
    idx_rela_list = [(smi_dict[rela[0]], smi_dict[rela[1]]) 
        for rela in rela_list]

    return smi_dict, smi_list, idx_rela_list


def anylseMD(freq_criter=2, NC_range=[0,5], Natom_range=None):
    relat_list, reac_list, reac_traj_list = getReacRelat()
    
    main_reac_list = filterRare(reac_list, freq_criter)
    printReacInfo(reac_list, main_reac_list, reac_traj_list)
    
    sort_relat_list = sortReacRelat(relat_list)
    sort_relat_list = filterRare(sort_relat_list, freq_criter)
    # sort_relat_list = clearRelat(sort_relat_list)

    filter_smi_dict = filterSMIList_atom(sort_relat_list, Natom_range)
    smi_list, filter_relat_list = filterRelatList_C(sort_relat_list, NC_range)
    smi_list, filter_relat_list = filterRelatList_atom(filter_relat_list, Natom_range)
    smi_dict, smi_list, idx_rela_list = getIdxList(smi_list, filter_relat_list)
    return smi_dict, smi_list, idx_rela_list, filter_smi_dict


if __name__ == '__main__':
    smi_dict, smi_list, idx_rela_list, filter_smi_dict = anylseMD()
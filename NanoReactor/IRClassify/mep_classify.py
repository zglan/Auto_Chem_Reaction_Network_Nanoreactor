#!/usr/bin/env python3

import numpy as np
import networkx as nx

from rdkit import Chem
from rdkit.Chem import Draw
from matplotlib.gridspec import GridSpec

from gau_api import *
from ioToolkit import *
from link_judge import linkJudge
from rewrite import SMILES, getBondfromcrd, getSMILESfromcrd

folder = Path('.')
log_list = list(folder.glob('*%s' % ('.log')))


def listCut(l, idx_list):
    newl = []
    for idx in idx_list:
        newl.append((l[idx]))
    return newl

def listSameCheck(lst1, lst2):
    if len(lst1) != len(lst2):
        return False
    for item in lst1:
        if item not in lst2:
            return False
    return True

def pathUnite(path_dict):
    same_list = []
    for i in path_dict.keys():
        irsmi_list, ipsmi_list, ileft_Ets, iright_Ets = path_dict[i]
        for j in path_dict.keys():
            if i != j:
                jrsmi_list, jpsmi_list, jleft_Ets, jright_Ets = path_dict[j]
                if listSameCheck(irsmi_list, jrsmi_list) or listSameCheck(irsmi_list, jpsmi_list):
                    if listSameCheck(ipsmi_list, jrsmi_list) or listSameCheck(ipsmi_list, jpsmi_list):
                        same_list.append((i, j))
    G = nx.Graph()
    G.add_edges_from(same_list)

    ulog_list = []
    Ealog_dict = {}
    connected = list(nx.connected_components(G))
    for log_set in connected:
        Ea_list = []
        for log in log_set:
            Ea_list.append(abs(path_dict[log][2] - path_dict[log][3]))
            ulog_list.append(log)

        Ea = np.average(Ea_list)
        Ealog_dict[Ea] = log_set
    
    for log in path_dict.keys():
        if log not in ulog_list:
            Ealog_dict[abs(path_dict[log][2] - path_dict[log][3])] = log

    sort_Ealog_list = [(Ea, Ealog_dict[Ea]) for Ea in sorted(Ealog_dict.keys())] 
    
    for idx, Ealog in enumerate(sort_Ealog_list):
        new_dir = Path('%.2f' % (Ealog[0]))
        mkdir(new_dir)
        if type(Ealog[1]) == set:
            for log in Ealog[1]:
                shutil.move(log, new_dir.joinpath(log))
        else:
            shutil.move(Ealog[1], new_dir.joinpath(Ealog[1]))

path_dict = {}
discard_dir = Path('discard')
mkdir(discard_dir)

for log in log_list:
    info = ReadIRCJob(log)
    try:
        atom_list = info.atomlist
        print(log)
        energyGap, energyPath, coordPath = info.getPath()

        right = coordPath[-1][-1]
        left = coordPath[1][-1]
        right_link, right_bonds = linkJudge(right)
        left_link, left_bonds = linkJudge(left)

        left_Ets = energyPath[0][-1] - energyPath[1][-1]
        right_Ets = energyPath[0][-1] - energyPath[-1][-1]

        rsmi_list = []
        for i in right_link:
            l = listCut(right, i)
            rsmi = getSMILESfromcrd(l).split('\t\n')[0]
            rsmi = SMILES.rewriteSMILES(rsmi, atom_list)
            rsmi_list.append(rsmi)  

        psmi_list = []
        for i in left_link:
            l = listCut(left, i)
            psmi = getSMILESfromcrd(l).split('\t\n')[0]
            psmi = SMILES.rewriteSMILES(psmi, atom_list)
            psmi_list.append(psmi)

        if listSameCheck(rsmi_list, psmi_list):
            print('No links changed: ' + str(log))
            shutil.move(log, discard_dir.joinpath(log))
        else:
            path_dict[log] = (rsmi_list, psmi_list, left_Ets, right_Ets)
    except:
        shutil.move(log, discard_dir.joinpath(log))

pathUnite(path_dict)


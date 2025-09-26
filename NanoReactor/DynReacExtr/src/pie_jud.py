from _smiles import *
from _toolkit import *
from _attr_mole import *
from _trajOper import *
from collections import Counter, defaultdict
from pathlib import Path

def getSMILES(xyz):
    with open(xyz) as f:
        file = f.readlines()
        atom_coord = OperTraj.oneFrame(file)
        bonds, bonds_order = getBondfromcrd(atom_coord)
        for index, line in enumerate(file):
            s = line.split()
            if index == 0:
                N = int(line.strip())
                atom_list = []
            elif index > N + 1:
                break
            elif index > 1:
                atom_list.append(s[0])
            atom_list = atom_list
        smi_list = [SMILES.node2SMILES_fbond(atom_list, fnodes, fbonds) for (fnodes, fbonds) in getConnection(bonds)]
        if smi_list == []:
            raise ValueError
        return smi_list

def classify(smi_collec):
    size_dict = {}
    for item in Counter(smi_collec).items():
        if item[0] == '[H]CN':
            size_dict['[H]CN'] = item[1]
        elif item[0] == '[H]O[H]':
            size_dict['[H]O[H]'] = item[1]
        else:
            atom_num = SMILES.molAtomNum(item[0])
            if atom_num <= 4:
                size_dict.setdefault('MolSize < 4 atoms', []).append(item[1])
            if 5 <= atom_num <= 9:
                size_dict.setdefault('MolSize = 5~9 atoms', []).append(item[1])
            if 9 < atom_num <= 15:
                size_dict.setdefault('MolSize = 10~15 atoms', []).append(item[1])
            if 15 <= atom_num <= 24:
                size_dict.setdefault('MolSize = 16~24 atoms', []).append(item[1])
            if atom_num >= 25:
                size_dict.setdefault('MolSize > 25 atoms', []).append(item[1])
    for i in size_dict.items():
        if type(i[1]) == list:
            size_dict[i[0]] = sum(i[1])
    print(size_dict)
    return size_dict

folder_list = getLocalFolderList()
subfolder_list = [getSubFolderList(folder) for folder in folder_list]
xyz_array = [getSuffixFileList(subfolder, 'xyz') for subfolder in subfolder_list]
group = {}
for t_group in xyz_array:
    for bias_group in t_group:
        typ = str(Path(bias_group[0]).parent).split('/')
        if typ[1] == '0':
            mark = 'MD (no bias)'# % (typ[0])# - %sps
        else:
            mark = 'MTD (ki = %s Eh)' % (typ[1])#, typ[0])# - %sps
        smi_collec = []
        for xyz in bias_group:
            smi_list = getSMILES(xyz)
            smi_collec += smi_list
        
        # for i in Counter(smi_collec).items():
        #     print(SMILES.molAtomNum(i[0]), i[0], i[1])
        group[mark] = classify(smi_collec)
        print()

print(group)
import matplotlib.pyplot as plt
import numpy as np

import numpy as np
import matplotlib.pyplot as plt

results = {}
category_names = ['MolSize < 4 atoms', 'MolSize = 5~9 atoms', 
    'MolSize = 10~15 atoms', 'MolSize = 16~24 atoms', 'MolSize > 25 atoms']
for item in group.items():
    results[item[0]] = [item[1][name]
        if name in item[1].keys() else 0
            for name in category_names]
    print(results)


def survey(results, category_names):
    """
    Parameters
    ----------
    results : dict
        A mapping from question labels to a list of answers per category.
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_names*.
    category_names : list of str
        The category labels.
    """
    labels = list(results.keys())
    data = np.array(list(results.values()))
    print()
    sum = [np.sum(data[i]) for i in range(len(data))]
    rate = np.array([data[i] / sum[i] for i in range(len(data))])
    print(rate)
    data_cum = data.cumsum(axis=1)
    rate_cum = rate.cumsum(axis=1)
    print(rate_cum)
    category_colors = plt.colormaps['RdYlGn'](
        np.linspace(0.15, 0.85, data.shape[1]))

    fig, ax = plt.subplots(figsize=(10, 3))
    ax.invert_yaxis()
    # ax.xaxis.set_visible(False)
    ax.set_xticks(np.arange(0, np.sum(data, axis=1).max(), 5))
    plt.xticks(fontsize=12)

    ax.set_xlim(0, np.sum(data, axis=1).max())
    ax.set_xlabel('Number of Species', fontsize=16)
    plt.yticks(fontsize=16, fontweight='bold', fontstyle='italic')

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        rects = ax.barh(labels, widths, left=starts, height=0.7,
                        label=colname, color=color)
        print(rects)
        # tmp_widths = rate[:, i]
        # tmp_starts = rate_cum[:, i] - tmp_widths
        # tmp_rects = ax.barh(labels, tmp_widths, left=tmp_starts, height=0.7,
        #                 label=colname, color=color)
        r, g, b, _ = color
        text_color = 'black' if r * g * b < 0.5 else 'black'
        ax.bar_label(rects, label_type='center', color=text_color,fontsize=14)
    ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1),
              loc='lower left', fontsize=18)

    return fig, ax


survey(results, category_names)
plt.show()
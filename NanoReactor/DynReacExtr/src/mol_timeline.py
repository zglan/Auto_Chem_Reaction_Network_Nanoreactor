from collections import Counter,defaultdict

import cmocean

import numpy as np

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from tqdm import tqdm

from _smiles import SMILES
from _toolkit import *

folder_list = getLocalFolderList()
data_list = extendFolderList(folder_list, '.ini_analy.npy')

interval = 2
fs2ps = 1000

bounds = [1, 2, 5, 10, 20, 35]
cmap = cmocean.cm.matter

norm = mpl.colors.Normalize(vmin=1, vmax=35)

# cmap = plt.colormaps["plasma"]
# norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='neither')

n = 1

graph_data_list = []

for data in data_list:
    data = np.load(data, allow_pickle=True)
    atom_list = data[0]
    line_num = data[1]
    bond_datas = data[2]
    step_num = data[3]

    if step_num == 15000:

        if n < 2:
            n += 1 
        else:
            break

        stepNum_list = [i / fs2ps for i in range(0, step_num * interval, interval)]

        maxWt_list = []
        maxAtomNum_list = []
        molNum_list = []

        tot_atomNum_list = []

        for bond_data in tqdm(bond_datas):
            bonds = [i[0] for i in bond_data]
            bond_orders = [i[1] for i in bond_data]

            smi_list = SMILES.convertSMILES(atom_list, bonds, bond_orders).split('.')

            atomNum_list = [SMILES.molAtomNum(smi) for smi in smi_list]
            atomNum = [[k,v] for k,v in Counter(atomNum_list).items()]
            tot_atomNum_list.append(atomNum)

            # wt_list = [SMILES.molWt(smi) for smi in smi_list]
            # maxWt = max([SMILES.molWt(smi) for smi in smi_list])
            # maxAtomNum = max([SMILES.molAtomNum(smi) for smi in smi_list])

            # maxAtomNum_list.append(maxAtomNum)
            # maxWt_list.append(maxWt)
            # molNum_list.append(len(smi_list))

        # graph_data = zip(tot_atomNum_list, stepNum_list)
        graph_data = tot_atomNum_list
        graph_data_list.append(graph_data)

def merge_lists(lst):
    d = defaultdict(int)
    for l in lst:
        for k, v in l:
            d[k] += v
    return [[k, v] for k, v in d.items()]

uni_graph_data = []

for step in list(range(0, step_num)):
    tmp_data = merge_lists(
        [graph_data_list[graph_idx][step] 
            for graph_idx in list(range(0, len(graph_data_list)))]
            )
    uni_graph_data.append(tmp_data)

fig, ax = plt.subplots()

def jugde(num):
    if num < 3:
        return 1.5
    elif num < 6:
        return 4.5
    elif num < 9:
        return 7.5
    elif num < 12:
        return 10.5
    elif num < 15:
        return 13.5
    elif num < 22:
        return 18.5
    elif num < 30:
        return 26
    elif num < 40:
        return 35
    
# size_dict = {}
# for atomNum_list, stepNum in tqdm(zip(uni_graph_data, stepNum_list)):
#     for atom, num in atomNum_list:
#         jnum = jugde(num)
#         size_dict.setdefault(jnum, []).append(stepNum)
#         ax.plot(stepNum, atom, 'o', color=cmap(num/30/n))


for atomNum_list, stepNum in tqdm(zip(uni_graph_data, stepNum_list)):
    for atom, num in atomNum_list:
        ax.plot(stepNum, atom, '.', color=cmap(num/15/n))

color_bar = plt.colorbar(
    cm.ScalarMappable(cmap=cmap, norm=norm),
    ax=ax,
    ticks=bounds
    # pad=0
)
color_bar.ax.tick_params(
    labelsize=16
)
color_bar.set_label(
    'Number of molecules ', 
    size=20,
    loc='center'
)

ax.tick_params(
    direction='in', 
    labelsize=16
)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlim(0, 30)
ax.set_ylim(0, 40)
ax.set_yticks([1,5,10,15,20,25,30,35,40])
ax.set_xlabel("Time (ps)", size=22)
ax.set_ylabel("Size of molecules (number of atoms)", size=22)

fig.tight_layout()

plt.show()
r"""
Reference: 

N. M. O'Boyle, A. L. Tenderholt, K. M. Langner, cclib: a library for package-independent computational chemistry algorithms, 
J. Comp. Chem. 29 (5), pp. 839-845, 2008

cclib package
"""

import os, re
import glob
import cclib
import numpy as np

from rdkit import Chem
from rdkit.Chem import Draw

import openbabel
from openbabel import pybel

def save_xyz(coords, elements, filename):
    with open(filename, 'w') as f:
        f.write(f"{len(elements)}\n\n")
        for elem, xyz in zip(elements, coords):
            f.write(f"{elem} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}\n")

# irc_files = glob.glob("*_*-irc.log")
irc_files = glob.glob('irc/**/Ener-*/Freq-*/*irc.log', recursive=True)

all_coords = []
elements = None

# print(irc_files)
# exit()

for file in irc_files:
    data = cclib.io.ccopen(file).parse()
    Ener_id = re.search(r"Ener-(\d+)", file).group(1)
    Freq_id = re.search(r"Freq-(\d+)", file).group(1)
    print(Ener_id, Freq_id)
    filename = f"irc-{Ener_id}_{Freq_id}"
    
    if not elements:
        elements = [data.atomnos[i] for i in range(data.natom)]
    all_coords.append(data.atomcoords)

    irc_coordinates = np.concatenate(all_coords, axis=0)

    reverse_end = len(all_coords[0]) - 1
    forward_end = len(irc_coordinates) - 1

    # print(irc_coordinates)
    # exit()

    ts_structure = irc_coordinates[0]          # TS结构
    reactant_structure = all_coords[0][reverse_end]  # 反向路径终点为反应物
    product_structure = irc_coordinates[forward_end] # 正向路径终点为产物

    # 保存结构
    save_xyz(ts_structure, elements, f"{filename}-ts.xyz")
    save_xyz(reactant_structure, elements, f"{filename}-reactant.xyz")
    save_xyz(product_structure, elements, f"{filename}-product.xyz")


def xyz2smiles(filename: str) -> str:
    mol = next(pybel.readfile("xyz", filename))

    smi = mol.write(format="smi")

    return smi.split()[0].strip()

xyz_files = glob.glob('./*.xyz', recursive=True)

# smi_list = []
mol_list = []
for xyz in xyz_files:
    smi = xyz2smiles(xyz)
    # smi_list.append(smi)
    mol = Chem.MolFromSmiles(smi)
    mol.SetProp("_Name", os.path.basename(xyz))
    mol_list.append(mol)


img = Draw.MolsToGridImage(mol_list, molsPerRow=8,subImgSize=(200,200),legends=[x.GetProp("_Name") for x in mol_list])  
img.save('molecular.png')

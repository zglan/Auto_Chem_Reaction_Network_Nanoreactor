#!/usr/bin/env python3

import re

from tqdm import tqdm
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
RDLogger.DisableLog('rdApp.*')

from md_analy import anylseMD


def stepGroup(l, step):
    group = [l[i: i + step] for i in range(0, len(l), step)]
    return group

def format_chemical_formula(formula):
    """Convert numbers in chemical formulas to subscripts"""
    subs = '₀₁₂₃₄₅₆₇₈₉'
    return re.sub(r'\d+', lambda m: ''.join(subs[int(d)] for d in m.group()), formula)

def drawGrid(mols, legends, img_name='', img_size=(1500, 900), font_size=150):
    # nRows = len(mols) // molsPerRow
    # if len(mols) % molsPerRow:
    #     nRows += 1
    nRows = 10
    fullSize = (molsPerRow * img_size[0], nRows * img_size[1])
    d2d = rdMolDraw2D.MolDraw2DSVG(fullSize[0], fullSize[1], img_size[0], img_size[1])
    # d2d.drawOptions().maxFontSize = 50
    # d2d.drawOptions().addAtomIndices = True
    # d2d.drawOptions().scalingFactor=100
    # d2d.drawOptions().centreMoleculesBeforeDrawing = True
    # d2d.drawOptions().useBWAtomPalette()
    # d2d.drawOptions().legendFraction = 0.3
    
    d2d.drawOptions().minFontSize = 50
    d2d.drawOptions().legendFontSize=font_size
    d2d.drawOptions().annotationFontScale = 16
    d2d.drawOptions().scaleBondWidth = True
    d2d.drawOptions().padding = 0.0001
    d2d.drawOptions().bondLineWidth = 6

    # change color no work
    # opts = d2d.drawOptions()
    # opts.legendFontSize = font_size
    # opts.legendColour = (1.0, 0.0, 0.0)

    # change legend no work
    # formatted_legends = [format_chemical_formula(legend) for legend in legends]

    d2d.DrawMolecules(mols, legends=legends)
    
    d2d.FinishDrawing()
    p = d2d.GetDrawingText()
    open(f'mols{img_name}-RDkit.svg','w+').write(p)


def getMolsLegends(smis, group=0):
    mols = []
    legends = []
    atom_nums = []
    for idx, smi in enumerate(smis):
        # atom_nums.append(len(re.findall("[a-zA-Z]", smi)))        # [a-zA-Z] can not deal with the situation like Ca.
        mol = Chem.MolFromSmiles(smi)
        atom_count = mol.GetNumAtoms()          # atom numbers obtained by rdkit
        atom_nums.append(atom_count)
        if isinstance(mol, Chem.rdchem.Mol):
            legend = CalcMolFormula(mol)
            mols.append(mol)
            legends.append(f'({idx + 1 + group * maxImgNum})   {legend}')
        else:
            legend = ''
    print(f' * Max atom number:  {max(atom_nums)}')
    return mols, legends, max(atom_nums)

def judgeImgSize(max_atom_nums):
    if max_atom_nums <= 15:
        subImgSize = (400, 400)
        legendFontSize = 300
    elif max_atom_nums <= 25:
        subImgSize = (600, 600)
        legendFontSize = 300
    elif max_atom_nums <= 35:
        subImgSize = (1600, 1000)
        legendFontSize = 300
    elif max_atom_nums <= 45:
        subImgSize = (1700, 1000)
        legendFontSize = 300
    elif max_atom_nums <= 60:
        subImgSize = (2000, 1300)
        legendFontSize = 300
    elif max_atom_nums <= 200:
        subImgSize = (3000, 2500)
        legendFontSize = 300
    print(f' * Sub-Grid Image Size:  {subImgSize}')
    return subImgSize, legendFontSize


smi_dict, smi_list, idx_rela_list, filter_smi_dict \
    = anylseMD(freq_criter=0, NC_range=None, Natom_range=[0,20])

# smi_dict.update(d2)
# del smi_dict['key']

"""Remove SMILES with incorrect atomic valence bonds"""
def standardize_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return smiles
    except:
        return None 

smi_dict_new = {}
idx = 1
for smi in smi_dict.keys():
    standardized = standardize_smiles(smi)
    if standardized:
        smi_dict_new[standardized] = idx
        idx = idx + 1

idx_mapping = {}
remove_idx = []
remove_smiles = []
for smiles, old_idx in smi_dict.items():
    if smiles in smi_dict_new.keys():
          new_idx = smi_dict_new[smiles]
          idx_mapping[old_idx] = new_idx
    else:
         remove_idx.append(old_idx)
         remove_smiles.append(smiles)

smi_list_new = []
for smi in smi_list:
    if smi not in remove_smiles:
        smi_list_new.append(smi)

idx_rela_list_new = []
for le, re in idx_rela_list:
     if le in remove_idx or re in remove_idx:
          pass
     else:
          le_new = idx_mapping[le]
          re_new = idx_mapping[re]
          idx_rela_list_new.append((le_new, re_new))

smi_dict = smi_dict_new
smi_list = smi_list_new
idx_rela_list = idx_rela_list_new



from _network import ReacEveNet
R = ReacEveNet(smi_dict, smi_list, idx_rela_list)
R.setLayoutParam({'k':0.6})
R.drawNet()


molsPerRow = 8
# maxImgNum = molsPerRow * 16
maxImgNum = molsPerRow * 10

tot_smis = [smi for smi in smi_dict.keys()]
if len(tot_smis) > maxImgNum:
    smis_list = stepGroup(tot_smis, maxImgNum)
    for group, smis in enumerate(tqdm(smis_list)):
        print(f'Drawing _Grid_Image_{group}......')
        mols, legends, max_atom_nums = getMolsLegends(smis, group)
        subImgSize, legendFontSize = judgeImgSize(max_atom_nums)
        drawGrid(mols, legends, f'{group + 1}', subImgSize, legendFontSize)
else:
    print(f'Drawing _Grid_Image......')
    mols, legends, max_atom_nums = getMolsLegends(tot_smis)
    subImgSize, legendFontSize = judgeImgSize(max_atom_nums)
    drawGrid(mols, legends, img_size=subImgSize, font_size=legendFontSize)


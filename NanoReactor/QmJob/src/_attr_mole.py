"""
Obtain molecluar attributes \
  only based on atom types and atom coordinates
      Attributes include:
        - default electronic numbers, 
        - molecular connectivity (bond, bond order),
        - molecular corresponding SMILES codes.
Caculate molecular attributes using GFN2-xTB method
      Attributes include:
        - Charge, 
        - Bond Order (wbo).
Refine molecular attributes by cumstom rules.
      Attributes include: 
        - Charge, 
        - Spin Multiplicity,
        - Bond Order (wbo).
"""

import os

from rdkit import Chem
from openbabel import openbabel

from ._toolkit import mkdir, coordStr2fxyz, atom_coord2str
from ._api_xtb import Xtb_Set
from ._setting_share import SharedSetting

#========================================================#
#                Obtain Info from atom_coord             #
#========================================================#

metal_sym = ['Li','Na','K']

def getAtomEleInfo(atom_coord, idx_list=None):
    ele_nums = 0
    atom_nums = 0
    atom_list = []
    for idx, atom in enumerate(atom_coord):
        if idx_list != None and idx in idx_list:
            ele_nums += Chem.Atom(atom[0].split()[0]).GetAtomicNum()
            atom_list.append(atom[0])
            atom_nums += 1
        elif idx_list == None:
            ele_nums += Chem.Atom(atom[0].split()[0]).GetAtomicNum()
            atom_list.append(atom[0])
            atom_nums += 1
    return atom_list, ele_nums, atom_nums

def getBondfromcrd(atom_coord):
    """
    Modified from ReacNetGenerator
    -----
    https://github.com/tongzhugroup/reacnetgenerator
    Doi: 10.1039/C9CP05091D

    :return bond -> array:
        [[atom1, atom2], [atom2, atom3]...]
    :return bond_order -> list:
    
    """
    mol = openbabel.OBMol()
    mol.BeginModify()
    for idx, (sym, position) in enumerate(atom_coord):
        a = mol.NewAtom(idx)
        if sym not in metal_sym:
            a.SetAtomicNum(Chem.Atom(sym).GetAtomicNum())
            a.SetVector(*position)
    # Adds single bonds based on atom proximity. 
    mol.ConnectTheDots()
    # Attempts to perceive multiple bonds based on geometries. 
    # This method uses bond angles and geometries from current connectivity 
    #   to guess atom types and then filling empty valences with multiple bonds. 
    #   It currently has a pass to detect some frequent functional groups. 
    #   It still needs a pass to detect aromatic rings to "clean up." 
    #   AssignSpinMultiplicity(true) is called at the end of the function. 
    #   The true states that there are no implict hydrogens in the molecule. 
    mol.PerceiveBondOrders() 
    mol.EndModify()

    bond = []
    bond_order = []
    for b in openbabel.OBMolBondIter(mol):
        s1 = b.GetBeginAtom().GetId()
        s2 = b.GetEndAtom().GetId()
        bond.append((s1, s2))
        order = b.GetBondOrder()
        if order == 5:
            # Aromatic: in openbabel order is 5, but in rdkit is 12
            order = 12
        bond_order.append(order)
    return bond, bond_order


def getSMILESfromcrd(atom_coord):
    con = openbabel.OBConversion()
    con.SetOutFormat("smi")
    con.SetOptions('-c', con.OUTOPTIONS)
    con.SetOptions('-h', con.OUTOPTIONS)

    mol = openbabel.OBMol()
    mol.BeginModify()

    for idx, (sym, position) in enumerate(atom_coord):
        a = mol.NewAtom(idx)
        if sym not in metal_sym:
            a.SetAtomicNum(Chem.Atom(sym).GetAtomicNum())
            a.SetVector(*position)
            
    mol.ConnectTheDots()
    mol.PerceiveBondOrders() 
    mol.EndModify()
    smi = con.WriteString(mol)
    return smi


#========================================================#
#                  Attributes Calculation                #
#========================================================#

class AttrCalc_xTB(SharedSetting):
    def __init__(self, tatom_coord, node_list, t_reac, t):
        super().__init__()
        self.tatom_coord = tatom_coord
        self.node_list = node_list
        self.t_reac, self.t = t_reac, t

    def getFragmentInfo(self):
        """
        Calculate fragment molecular attributes 
          using xTB program.
        
        :return wbo_list -> array:
            [ [[atom_1, atom_2], wbo_value -> float],
              [[atom_3, atom_4], wbo_value -> float]...]
        :return chrg_list -> list:
            fragment conrresponding atomic charge list 
                [chrg1 -> float, chrg2 -> float, ...]
        :return chrg -> float:
            framgent's total charge.
        :return ele_nums -> int:
            framgent's total default electron numbers.

        """
        tatom_list, tele_nums, tatom_nums = getAtomEleInfo(
            self.tatom_coord)
        tcoord_str = atom_coord2str(self.tatom_coord)
        atom_list, ele_nums, atom_nums = getAtomEleInfo(
            self.tatom_coord, self.node_list)

        src = os.getcwd()
        dst = 'sp/%s/' % (self.t_reac)
        xyz = '%s.xyz' % (self.t)

        mkdir(dst)
        coordStr2fxyz(dst + xyz, tcoord_str, tatom_nums)
        
        os.chdir(dst)
        xtb_spin = self.tspin - 1 # xtb spin transform
        wbo_list, chrg_list, chrg = Xtb_Set(
            xyz, 
            chrg=self.tchrg, 
            spin=xtb_spin
        ).getFragmentInfo(self.node_list)
        os.chdir(src)
        
        return wbo_list, chrg_list, chrg, ele_nums


class AttrCalc_ase:
    ### Exists Problem, Banned !!! ###
    def __init__(self, atom_list, coord):
        from ase import Atoms
        self.mole = Atoms(atom_list)
        self.mole.set_positions(coord)
        self.atomsCharge = self.getAtomsCharge()

    def getAtomsCharge(self):
        from xtb.ase.calculator import XTB
        self.mole.calc = XTB(
            method="GFN2-xTB", max_iterations=300)
        return self.mole.get_charges()
    
    def getFragmentCharge(self, idx_list):
        import numpy as np
        fragmentCharge = np.sum(
            [self.atomsCharge[idx] for idx in idx_list])
        return fragmentCharge


#========================================================#
#               Jugde molecular attributes               #
#========================================================#

def setChrg(chrg):
    """
    Judge charge based on custom rule.

    :param chrg -> float:
    :return chrg -> int:
 
    """
    rules = [
        (-0.45,  0.45,  0), 
        ( 0.45,  1.5,  1), 
        ( 1.5,  2.5,  2), 
        (-1.5, -0.45, -1), 
        (-2.5, -1.5, -2)
    ]
    for rule in rules:
        if rule[0] < chrg <= rule[1]:
            return rule[2]
    return int(chrg)


def setSpin(spin):
    """
    Judge spin based on custom rule.

    :param spin -> float:
    :return chrg -> int:
 
    """
    rules = [
        (-0.5,  0.5,  1), 
        ( 0.5,  1.5,  2), 
        (-1.5, -0.5,  2),
        ( 1.5,  2.5,  3), 
        (-2.5, -1.5,  3)
    ]
    for rule in rules:
        if rule[0] < spin <= rule[1]:
            return rule[2]  
    return int(spin)


def setChrgSpin_ele(ele_nums, 
        chrg=None, spin=None):
    """
    Judge charge and spin based on custom rule \
        with electron numbers.

    :param elenums -> int: 
        electron numbers
    :param chrg -> float/none/int:
    :param spin -> float/none/int:

    :return chrg -> int:
    :return spin -> int:
    
    """
    if chrg != None and spin != None:
        chrg = int(chrg)
        spin = int(spin)
    elif chrg != None:
        chrg = setChrg(chrg)
        if (ele_nums - chrg) % 2 == 0:
            spin = 1
        elif (ele_nums - chrg) % 2 != 0:
            spin = 2
    return chrg, spin


def setBO_Wbo(wbo_list):
    """
    Judge wbo bond order based on custom rule.

    :param wbo_list: 
        [ [[atom_1, atom_2], wbo_value -> float],
          [[atom_3, atom_4], wbo_value -> float]...]
    :return nwbo_list:
        [ [[atom_1, atom_2], wbo_value -> int],
          [[atom_3, atom_4], wbo_value -> int]...]
    """
    nwbo_list = []
    for wbo in wbo_list:
        if 0.5 <= wbo[-1] < 1.5:
            nwbo_list.append([wbo[0], 1])
        elif 1.5 <= wbo[-1] < 2.40:
            nwbo_list.append([wbo[0], 2])
        elif 2.40 <= wbo[-1] < 3.40:
            nwbo_list.append([wbo[0], 3])
    return nwbo_list
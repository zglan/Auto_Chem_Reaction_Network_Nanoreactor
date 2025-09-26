import re


from rdkit import Chem
from rdkit.Chem import Draw

from _attr_mole import setChrg


class SMILES:
   
    def node2SMILES_fbond(tatom_list, fnodes, fbonds):
        """
        Convert fragment node to SMILES code using
          entire atom list, fragment nodes, fragment bond connection
        """
        node2idx = {}
        fatom_list, fbond_list  = [], []
        for idx, node in enumerate(fnodes):
            node2idx[node] = idx
            fatom_list.append(tatom_list[node])
        for bond in fbonds:
            fbond_list.append(
                (node2idx[bond[0]], node2idx[bond[1]]))
        return SMILES.convertSMILES(fatom_list, fbond_list)

    def node2SMILES_wbo(tatom_list, fnodes, 
                        twbo_list, tchrg_list):
        """
        Convert fragment node to SMILES code using
          entire atom list, entire wbo bond order list, 
          entire atomic chrg list
        """
        node2idx = {}
        fatom_list, fwbo_list, fchrg_list  = [], [], []
        for idx, node in enumerate(fnodes):
            node2idx[node] = idx
            fatom_list.append(tatom_list[node])
            for chrg in tchrg_list:
                if node == chrg[-1]:
                    fchrg_list.append([
                        chrg[0], node2idx[chrg[-1]]
                            ])
        for idx, node in enumerate(fnodes):
            for wbo in twbo_list:
                if node in wbo[0] and wbo[0][0] in fnodes and wbo[0][-1] in fnodes:
                    fwbo_list.append((
                        (node2idx[wbo[0][0]], node2idx[wbo[0][-1]]), 
                            wbo[-1]))
        fwbo_list = list(set(fwbo_list))
        return SMILES.rewriteSMILES(
            SMILES.convertSMILES_re(fatom_list, fwbo_list, fchrg_list), 
            tatom_list)

    def fixSMILES(smi):
        repl_smi = {
            '[HH]': '[H]', 
            '[NaH]': '[Na]', 
            '[SH]': '[S]'
        }
        for ini_smi in repl_smi:
            smi = smi.replace(ini_smi, repl_smi[ini_smi])
        return smi

    def rewriteSMILES(smi, atom_list):
        atom_type = set(atom_list)
        Satom = sorted(atom_type, key = lambda i:len(i), reverse=True)
        elements = "|".join([(
            (atom.upper() + "|" + atom.lower()) 
                if len(atom)==1 else atom) 
            for atom in Satom if atom != 'H'])
        smi = re.sub(r'(?<!\[)(' + elements + r')(?!H)', r'[\1]', smi)
        return smi
    
    def removeSMILES_BO(smi):
        return smi.replace('=', '').replace('#', '')

    def convertSMILES(atom_list, bonds, bond_orders=None):
        m = Chem.RWMol(Chem.MolFromSmiles(''))
        for atom in atom_list:    
            m.AddAtom(Chem.Atom(atom))
        if bond_orders == None:
            for bond in bonds:
                m.AddBond(bond[0], bond[1], Chem.BondType(1))
        else:
            for bond, bond_orders in zip(bonds, bond_orders):
                m.AddBond(bond[0], bond[1], Chem.BondType(bond_orders))
        mole = Chem.MolToSmiles(m)
        return mole
    
    def convertSMILES_re(atom_list, wbo_list, chrg_list):
        pt = Chem.GetPeriodicTable()
        m = Chem.RWMol(Chem.MolFromSmiles(''))

        fchrg = 0
        for chrg, idx in chrg_list:
            fchrg += chrg  
        fchrg = setChrg(fchrg)
        for atom in atom_list:    
            m.AddAtom(Chem.Atom(atom))
        else:
            for bond, bond_order in wbo_list:
                m.AddBond(bond[0], bond[1], Chem.BondType(bond_order))
        m.UpdatePropertyCache(strict=False)

        atom_rwmol = [atom for atom in m.GetAtoms()]
        for idx, atom in enumerate(atom_rwmol):
            if atom.GetExplicitValence() > pt.GetDefaultValence(atom.GetAtomicNum()):
                gap = int(atom.GetExplicitValence() - pt.GetDefaultValence(atom.GetAtomicNum()))
                m.GetAtomWithIdx(idx).SetFormalCharge(gap)
            if atom.GetExplicitValence() < pt.GetDefaultValence(atom.GetAtomicNum()) and fchrg <= -1:
                # gap = int(atom.GetExplicitValence() - pt.GetDefaultValence(atom.GetAtomicNum()))
                # m.GetAtomWithIdx(idx).SetFormalCharge(gap)
                m.GetAtomWithIdx(idx).SetFormalCharge(fchrg)
                print(m)
                break
        smi = Chem.MolToSmiles(m)
        
        if smi == '[HH]':
            smi = '[H]'

        return smi
    
    def drawSMILES(smi, loc, sanitize=False):
        rdmol = Chem.MolFromSmiles(smi, sanitize) 
        Draw.MolToFile(rdmol, loc, size=(500, 500), kekulize=True)

    def molWt(smi):
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        wt = 0
        for atom in mol.GetAtoms():
            wt += Chem.GetPeriodicTable().GetAtomicWeight(atom.GetAtomicNum())
        return wt
    
    def molAtomNum(smi):
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        num = len(mol.GetAtoms())
        return num
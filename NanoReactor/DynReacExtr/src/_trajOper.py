import os
import itertools

import numpy as np
import pandas as pd
import networkx as nx

from tqdm import tqdm
from rdkit import Chem

from _api_gau import G16_Set
from _api_orca import Orca_Set
from _api_xtb import Xtb_Set
from _setting_share import SharedSetting
from _smiles import SMILES
from _attr_mole import getAtomEleInfo
from _toolkit import (
        getConnection, 
        atom_coord2array, atom_coord2str, array2atom_coord, 
        mkdir, coordStr2fxyz, 
        diffNetworks
)
from _attr_mole import (
        AttrCalc_xTB, 
        getBondfromcrd,
        setChrgSpin_ele, setBO_Wbo
)

class OperTraj(SharedSetting):

    def __init__(self):
        SharedSetting.__init__(self)
        if os.path.isfile(self.loc_config['ini_read_data']):
            self.loadIniData()
        else:
            self.readLineAtom()
            self.readFullTraj()

    def loadIniData(self):
        self.atom_list, self.step_linenum, \
            self.bond_data, self.step, self.smi_nodes = np.load(
                self.loc_config['ini_read_data'], allow_pickle=True)

    def readLineAtom(self):
        """
        Modified from ReacNetGenerator
        -----
        https://github.com/tongzhugroup/reacnetgenerator
        Doi: 10.1039/C9CP05091D
        """

        with open(self.inputfile) as f:
            for index, line in enumerate(f):
                s = line.split()
                if index == 0:
                    N = int(line.strip())
                    atom_list = []
                elif index > N + 1:
                    break
                elif index > 1:
                    atom_list.append(s[0])
            self.step_linenum = N + 2
            self.atom_list = atom_list
    
    def iterObject(self, f):
        """
        Modified from ReacNetGenerator
        -----
        https://github.com/tongzhugroup/reacnetgenerator
        Doi: 10.1039/C9CP05091D
        """
        
        obj = itertools.zip_longest(*[f]*self.step_linenum)
        obj = itertools.islice(obj, 0, None, self.interval)
        obj = enumerate(obj, 0)
        
        return obj

    @staticmethod
    def oneFrame(lines):
        """
        Modified from ReacNetGenerator
        -----
        https://github.com/tongzhugroup/reacnetgenerator
        Doi: 10.1039/C9CP05091D
        """
        
        atom_coord = []      
        for index, line in enumerate(lines):
            s = line.split()
            if index > 1:
                atom_coord.append((s[0], [float(x) for x in s[1:4]]))
        return atom_coord

    def readFullTraj(self):
        with open(self.inputfile) as f:
            bond_data = []
            smi_nodes = {}
            obj = self.iterObject(f)
            while True:
                try:
                    step, lines = next(obj)
                    atom_coord = self.oneFrame(lines)
                    bonds, bonds_order = getBondfromcrd(atom_coord)
                    bond_data.append(list(zip(bonds, bonds_order)))
                    for (fnodes, fbonds) in getConnection(bonds):
                        smi = SMILES.node2SMILES_fbond(
                            self.atom_list, fnodes, fbonds)
                        smi_node = f'''{smi}_{str(fnodes)}'''
                        smi_nodes.setdefault(smi_node, []).append(step)
                except StopIteration:
                    break
        
        self.bond_data = bond_data
        self.step = step + 1
        self.smi_nodes = smi_nodes
        
        save_data = np.array(
            [self.atom_list, self.step_linenum,
            self.bond_data, self.step, self.smi_nodes], 
            dtype=object)

        np.save(
            self.loc_config['ini_read_data'], 
            save_data)
    
    @staticmethod
    def extractOneFrame(lines, job):
        ele_nums = 0
        atom_nums = 0
        coord_str = ''
        for idx, line in enumerate(lines):
            if idx - 2 in job[0][1]:
            # if idx > 1:
                ele_nums += Chem.Atom(line.split()[0]).GetAtomicNum()
                atom_nums += 1
                coord_str += line
        return coord_str, ele_nums, atom_nums
    
    def extractG16Job(self, job, t_list, type):
        info_list = []
        chrg_list = []
        with open(self.inputfile) as f:
            obj = self.iterObject(f)
            while True:
                ### list [step, lines] ###
                inf = next(obj, None)
                if inf == None:
                    break
                elif inf[0] in t_list:
                    atom_coord = self.oneFrame(inf[-1])
                    coord_str, \
                    ele_nums, atom_nums = self.extractOneFrame(inf[-1], job)

                    wbo_list, fchrg_list, \
                    fchrg, fele_nums = AttrCalc_xTB(
                        atom_coord, job[0][1], job[-1], inf[0]
                    ).getFragmentInfo()
                        
                    chrg_list.append(fchrg)
                    coord_str = atom_coord2str(atom_coord, job[0][1])
                    
                    info_list.append((
                        '%s/%s/%s.gjf' % (type, job[-1], inf[0]), type, 
                        fchrg, fele_nums, coord_str))

        fchrg = np.average(chrg_list)
        fchrg, fspin = setChrgSpin_ele(fele_nums, fchrg)

        for info in info_list:
            G16_Set(
                gjf=info[0], 
                job_type=info[1],
                coord=info[4],
                chrg=fchrg,
                spin=fspin
            ).setGjf()

    def runXtbOptJob(self, job, t_list, type):
        info_list = []
        chrg_list = []
        with open(self.inputfile) as f:
            obj = self.iterObject(f)
            while True:
                ### list [step, lines] ###
                data = next(obj, None)
                if data == None:
                    break
                elif data[0] in t_list:
                    t = data[0]
                    xyz_str = data[-1]

                    atom_coord = self.oneFrame(xyz_str)
                    coord_str, ele_nums, atom_nums = self.extractOneFrame(xyz_str, job)
                    bond_ori, bond_order = getBondfromcrd(atom_coord)

                    wbo_list, fchrg_list, fchrg, fele_nums = AttrCalc_xTB(
                        atom_coord, job[0][1], job[-1], t).getFragmentInfo()
                    chrg_list.append(fchrg)
                    coord_str = atom_coord2str(atom_coord, job[0][1])
                    fatom_list, ele_nums, fatom_nums = getAtomEleInfo(atom_coord)

                    if t in t_list[:3] or t in t_list[-3:]:
                        info_list.append((
                            '%s.xyz' % (t), type, bond_ori, bond_order,
                        fchrg, fele_nums, coord_str, atom_nums, fatom_list))

        fchrg = np.average(chrg_list)
        fchrg, fspin = setChrgSpin_ele(fele_nums, fchrg)
        xtb_spin = fspin - 1

        dst = '%s/%s/' % (type, job[-1])
        os.chdir(dst)

        data = []
        for idx, info in enumerate(info_list):
            xyz_traj, bond_ori, bond_order = info[0], info[2], info[3]
            coord_str, atom_nums, fatom_list = info[-3], info[-2], info[-1]
            coordStr2fxyz(xyz_traj, coord_str, atom_nums)

            with open(xyz_traj) as opt:
                f = opt.readlines()
                atom_coord = self.oneFrame(f)
                atom_list, ele_nums, atom_nums = getAtomEleInfo(atom_coord)
            bond_ori, bond_order_ori = getBondfromcrd(atom_coord)

            smi = SMILES.convertSMILES(atom_list, bond_ori)

            opt_xyz = Xtb_Set(
                xyz=xyz_traj,
                jobType='opt',
                chrg=fchrg,
                spin=xtb_spin
            ).getOptStruct()
            with open(opt_xyz) as opt:
                f = opt.readlines()
                atom_coord_opt = self.oneFrame(f)
            bond_opt, bond_order_opt = getBondfromcrd(atom_coord_opt)
            smi_opt = SMILES.convertSMILES(atom_list, bond_opt)

            G1 = nx.Graph(bond_ori)
            G1.add_nodes_from([i for i in range(atom_nums)])
            G2 = nx.Graph(bond_opt)
            G2.add_nodes_from([i for i in range(atom_nums)])

            G = diffNetworks(G1, G2)

            status = False if len(G.edges) > 0 else True
            row = {'Job_index': job[-1], 'Index': idx+1, 'Smi_equal': status, 
                        'Xyz': xyz_traj, 'Opt_Xyz': opt_xyz, 
                        'Chrg': fchrg, 'Spin': fspin, 
                        'SMI': smi, 'SMI_Opt': smi_opt}
            data.append(row)
        return data


    def refReaction(self, job, t_list):

        def crd2Info(inf, job):
            atom_coord = self.oneFrame(inf[-1])
            fgt_wbo_list, fgt_chrg_list, \
            fgt_chrg, fele_nums = AttrCalc_xTB(
                atom_coord, job[0][1], job[-1], inf[0]
                ).getFragmentInfo()
            fgt_wbo_list = setBO_Wbo(fgt_wbo_list)
            G = nx.Graph()
            for fgt_wbo in fgt_wbo_list:
                if fgt_wbo[0][0] in job[0][1] and fgt_wbo[0][1] in job[0][1]:
                    G.add_edges_from([fgt_wbo[0]])
            for node in job[0][1]:
                G.add_node(node)
            return fgt_wbo_list, fgt_chrg_list, fgt_chrg, G
                    
        with open(self.inputfile) as f:
            obj = self.iterObject(f)
            while True:
                ### list [step, lines] ###
                inf = next(obj, None)
                if inf == None:
                    break
                elif inf[0] == t_list[0]:
                    rfgt_wbo_list, rfgt_chrg_list,\
                    rfgt_chrg, G_r = crd2Info(inf, job)
                elif inf[0] == t_list[1]:
                    pfgt_wbo_list, pfgt_chrg_list, \
                    pfgt_chrg, G_p = crd2Info(inf, job)
           
            G_diff = diffNetworks(G_r, G_p)
            node_relation = []
            mole_relation = []

            for node in G_diff.nodes:
                r_nodes = tuple(sorted(
                    nx.node_connected_component(G_r, node)))
                p_nodes = tuple(sorted(
                    nx.node_connected_component(G_p, node)))

                r_mole = SMILES.node2SMILES_wbo(self.atom_list, r_nodes, 
                        rfgt_wbo_list, rfgt_chrg_list)
                p_mole = SMILES.node2SMILES_wbo(self.atom_list, p_nodes, 
                        pfgt_wbo_list, pfgt_chrg_list)

                if (r_nodes, p_nodes) not in node_relation:
                    node_relation.append((r_nodes, p_nodes))
                    mole_relation.append((r_mole, p_mole))

            return mole_relation

    def extractPath(self, job, t_list):
        coord_list = []
        chrg_list = []
        ele_list = []

        mid = job[-1]
        st_list = set(
            list(np.arange(t_list[mid], t_list[0], -self.jobset[-1])) + 
            list(np.arange(t_list[mid], t_list[-1], self.jobset[-1])))
        
        with open(self.inputfile) as f:
            obj = self.iterObject(f)
            while True:
                inf = next(obj, None)
                if inf == None:
                    break
                elif inf[0] in t_list:
                    atom_coord = self.oneFrame(inf[-1])
                    atom_array, coord_array = atom_coord2array(
                        atom_coord, job[0][1])
                    coord_list.append(coord_array)
                    
                    if inf[0] in st_list:
                        wbo_list, fchrg_list, \
                        chrg, ele_nums = AttrCalc_xTB(
                            atom_coord, job[0][1], job[-1], inf[0]
                        ).getFragmentInfo()
                        chrg_list.append(chrg)
                        ele_list.append(ele_nums)
        
        ele_nums = np.sum(ele_list)
        coord_array = np.array(coord_list)
        chrg = np.average(chrg_list)
                    
        return atom_array, coord_array, chrg, ele_nums

    def makePathPoints(self, job_type, job, t_list):

        r_t, p_t = job[0][2][0] - t_list[0], job[0][2][-1]- t_list[0]
        
        atom_array, coord_array, \
        fchrg, ele_nums = self.extractPath(job, t_list)

        if abs(fchrg) > 0.4 and abs(fchrg) < 0.6:
            os.system('echo \'%s\' >> ./chrgException' % (job[0][0]))
        
        fchrg, fspin = setChrgSpin_ele(ele_nums, fchrg)
        print(fchrg)
        atom_nums = len(atom_array)

        if self.smooth:
            try:
                from nebterpolator.path_operations import smooth_internal
                scoord_array = [coord_array, None] \
                    if atom_nums < 5 \
                        else smooth_internal(coord_array, atom_array, 5)
            except:
                print('* Need nebterpolator !')
        else:
            scoord_array = [coord_array, None]
        satom_coord_str = atom_coord2str(
            array2atom_coord(atom_array, scoord_array[0][r_t]))
        eatom_coord_str = atom_coord2str(
            array2atom_coord(atom_array, scoord_array[0][p_t]))

        sxyz, exyz = '%s/%s/%s.xyz' % (job_type, job[-1], 'start'), \
                '%s/%s/%s.xyz' % (job_type, job[-1], 'end')
        xyz = [sxyz.split('/')[-1], exyz.split('/')[-1]]

        mkdir('%s/%s/' % (job_type, job[-1]))
        coordStr2fxyz(sxyz, satom_coord_str, atom_nums)
        coordStr2fxyz(exyz, eatom_coord_str, atom_nums)

        if job_type == 'neb':
            Orca_Set(
                '%s/%s/%s.inp' % (job_type, job[-1], job_type), 
                job_type, 
                xyz, 
                chrg=fchrg, 
                spin=fspin).setInp()
        if job_type == 'path':
            xtb_spin = fspin - 1
            Xtb_Set(
                sxyz, 
                chrg=fchrg, 
                spin=xtb_spin, 
                jobType='path').makeInp()
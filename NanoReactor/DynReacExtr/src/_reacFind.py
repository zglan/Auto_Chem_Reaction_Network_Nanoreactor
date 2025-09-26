from collections import Counter

import numpy as np
import networkx as nx
import pandas as pd

from rdkit import Chem

from _setting_share import SharedSetting
from _smiles import SMILES
from _toolkit import diffNetworks, limitList

class FindReaction(SharedSetting):
    
    @staticmethod
    def formReaction(reactant, product):
        return SMILES.removeSMILES_BO(
            '+'.join(reactant) + '->' + '+'.join(product))
    
    @staticmethod
    def findMole(atom_list, bond_data, node_list, dict_nodeT):
        """
        Convert nodes_list to SMILES
        """

        smi_list = []
        for time in dict_nodeT[str(node_list)]:
            m = Chem.RWMol(Chem.MolFromSmiles(''))
            bond_involved = []
            for node in node_list:
                m.AddAtom(Chem.Atom(atom_list[node]))
                for b, bo in bond_data[time]:              
                    if node in b and b[0] not in bond_involved:
                        bond_involved.append((b, bo))
            bond_involved = set(bond_involved)
            for b, bo in bond_involved:
                m.AddBond(
                    node_list.index(b[0]), 
                    node_list.index(b[1]), 
                    Chem.BondType(bo))
            smi_list.append(tuple(
                    SMILES.fixSMILES(
                        Chem.MolToSmiles(m))
                    .split('.')))
        return tuple(sorted(
                    Counter(smi_list).most_common()[0][0]))

    def findConectNode(self, time, node_list, 
                            atom_list, bond_data, mode):
        all_node_list = []
        dict_nodeT = {}

        if mode == 'from reac':
            t_list = np.arange(-self.scale, -2, 2) + int(time)
        elif mode == 'from prod':
            t_list = np.arange(2, self.scale, 2) + int(time)
        t_list = limitList(t_list, len(bond_data) - 1, 0)

        for t in t_list:
            G = nx.Graph()
            G.add_nodes_from(range(len(atom_list)))
            G.add_edges_from([b for b, bo in bond_data[t]])
            
            node_set = set()
            for node in node_list:
                node_set |= nx.node_connected_component(G, node)
            node_set = sorted(list(node_set))
            all_node_list.append(node_set)
            dict_nodeT.setdefault(str(node_set), []).append(t)

        max_node_set = max(all_node_list, key=all_node_list.count)
        if set(node_list) != max_node_set:
            node_list = sorted(list(max_node_set))
        return node_list, dict_nodeT    

    def findAllNode(self, r_d_t_node, p_d_t_node, 
                    time, atom_list, bond_data):
        if self.mode == 'reac':
            node_list = r_d_t_node[time]
        elif self.mode == 'prod':
            node_list = p_d_t_node[time]

        contition = True
        while contition:
            if self.mode == 'reac':
                prodNode, prodDict_nodeT = self.findConectNode(time, node_list,
                        atom_list, bond_data, mode='from prod')
                reacNode, reacDict_nodeT = self.findConectNode(time, prodNode, 
                        atom_list, bond_data, mode='from reac')
                if reacNode == prodNode:
                    contition = False
                node_list = reacNode

            elif self.mode == 'prod':
                reacNode, reacDict_nodeT = self.findConectNode(time, node_list, 
                        atom_list, bond_data, mode='from reac')
                prodNode, prodDict_nodeT = self.findConectNode(time, reacNode, 
                        atom_list, bond_data, mode='from prod')
                if prodNode == reacNode:
                    contition = False
                node_list = prodNode
        return node_list, prodDict_nodeT, reacDict_nodeT

    def findRelation(self, prodDict_nodeT, reacDict_nodeT, node_list, 
                     atom_list, bond_data):
        
        def node2SMILES(node, connect):
            atom = []
            atom_connect = []
            node_dict = {}
            for idx, node in enumerate(node):
                node_dict[node] = idx
                atom.append(atom_list[node])
                for c in connect:
                    if node in c:
                        atom_connect.append(c)
            atom_connect = [(node_dict[c[0]], node_dict[c[1]]) 
                            for c in list(set(atom_connect))]
            smi = SMILES.removeSMILES_BO(
                    SMILES.rewriteSMILES(
                        SMILES.convertSMILES(
                            atom, atom_connect), atom_list))
            return smi

        def getSubgraph(bond_data, t):
            tbond_data = bond_data[t]
            G = nx.Graph()
            for bond in tbond_data:
                if bond[0][0] in node_list or bond[0][1] in node_list:
                    G.add_edge(bond[0][0], bond[0][1])
            for node in node_list:
                if node not in G.nodes():
                    G.add_node(node)
            return G

        rt = reacDict_nodeT[str(node_list)][0]
        pt = prodDict_nodeT[str(node_list)][-1]

        Gr_sub = getSubgraph(bond_data, rt)
        Gp_sub = getSubgraph(bond_data, pt)
        G_new = diffNetworks(Gr_sub, Gp_sub)
        
        mole_relation = []
        node_relation = []
        for node in G_new.nodes:
            r_node = tuple(sorted(
                nx.node_connected_component(Gr_sub, node)))
            p_node = tuple(sorted(
                nx.node_connected_component(Gp_sub, node)))         

            r_mole = node2SMILES(r_node, Gr_sub.edges())
            p_mole = node2SMILES(p_node, Gp_sub.edges())
            
            if (r_node, p_node) not in node_relation:
                node_relation.append((r_node, p_node))
                mole_relation.append((r_mole, p_mole))
            
        return mole_relation

    def findReaction(self, r_d_t_node, p_d_t_node, 
                        atom_list, bond_data):
        reaction_mole, reactions, job_list, relation_list = [], [], [], []

        if self.mode == 'reac':
            d_t_node = r_d_t_node
        elif self.mode == 'prod':
            d_t_node = p_d_t_node

        for time in d_t_node.keys():
            node_list, prodDict_nodeT, reacDict_nodeT = self.findAllNode(
                    r_d_t_node, p_d_t_node, time, 
                    atom_list, bond_data)
            r_time, p_time = reacDict_nodeT[
                str(node_list)][0], prodDict_nodeT[str(node_list)][-1]
            reac_mole = self.findMole(
                atom_list, bond_data, 
                node_list, reacDict_nodeT)
            prod_mole = self.findMole(
                atom_list, bond_data, 
                node_list, prodDict_nodeT)
            
            reaction = self.formReaction(
                reac_mole, prod_mole)
            relation = self.findRelation(
                prodDict_nodeT, reacDict_nodeT, 
                node_list, atom_list, bond_data)
            reaction_mole.append(
                (reac_mole, prod_mole))
            
            if reac_mole != prod_mole:
                judge = True
                if len(job_list) == 0:
                    job_list.append(
                        (reaction, node_list, (r_time, p_time), int(time)))
                    reactions.append(reaction)
                    relation_list.append((reaction, relation))
                else:
                    for job in job_list:
                        if node_list == job[1] \
                              and abs(int(time)- job[-1]) < 20 \
                              and reaction == job[0]:
                            judge = False
                    if judge:
                        job_list.append(
                            (reaction, node_list, 
                             (r_time, p_time), int(time)))
                        reactions.append(reaction)
                        relation_list.append((reaction, relation))
            
        reactions = Counter(reactions).items()
        for reaction in reactions:
            print('   *  %s times  |  %s'%(reaction[1], reaction[0]))

        columns = ['reaction', 'relation']
        df = pd.DataFrame(columns=columns, data=relation_list)
        df.to_csv('_relation.csv')
        
        columns = ['reaction', 'atom-idx', 'start-to-end', 'time']
        df = pd.DataFrame(columns=columns, data=job_list)
        df.to_csv('_job.csv')

        return job_list
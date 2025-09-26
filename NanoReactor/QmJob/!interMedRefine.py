#!/usr/bin/env python3

import ast
import subprocess

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import PIL
import dash
import dash_cytoscape as cyto

from dash import html, dcc
from dash.dependencies import Input, Output
from openbabel import openbabel
from openbabel import pybel

try:
    from rdkit import Chem
    from rdkit.Chem import Draw, rdMolDescriptors
except:
    print('need rdkit')

from src._toolkit import *
from src._logger import MyLog
from src._api_gau import G16_Set, G16_Read, G16_Ropt

from src._setting_share import SharedSetting

from src._set_pbs import buildJobPbs
from src._submit import submitJobs
from src._csv import recordStateCSV, recordJobCSV, csv2List

myLogger = MyLog().logger

class OPT_Job:

    def __init__(self, datas, idx):
        self.idx = idx
        self.reactants = datas[0]
        self.products = datas[1]
        self.reactants_no_eq = self.str2list(datas[2])
        self.products_no_eq = self.str2list(datas[3])
        self.reactants_eq = self.str2list(datas[4])
        self.products_eq = self.str2list(datas[5])
        self.counts = int(datas[-2])
        self.spin_chrg = ast.literal_eval(datas[-1])

    def str2list(self, data):
        data = ast.literal_eval(data)
        data = [(i[0], int(i[1]), int(i[2]))
            for i in data] if len(data) > 0 else []
        return data

    def pre_process(self, data_eq, data_no_eq):
        if data_no_eq == []:
            max_item = max(data_eq, key=lambda x: x[2])
            eq = True
        elif data_eq == []:
            max_item = max(data_no_eq, key=lambda x: x[2])
            eq = False
        else:
            if sum([x[2] for x in data_eq]) > \
                sum([x[2] for x in data_no_eq]):
                max_item = max(data_eq, key=lambda x: x[2])
                eq = True
            else:
                max_item = max(data_no_eq, key=lambda x: x[2])
                eq = False
        return max_item, eq

    def main(self):
        reactants, r_eq = self.pre_process(
            self.reactants_eq, self.reactants_no_eq)
        products, p_eq = self.pre_process(
            self.products_eq, self.products_no_eq)
        
        info = [self.idx, self.spin_chrg, self.counts,
            reactants[1], r_eq, products[1], p_eq]
        return  info


class G16_SP(SharedSetting):
    
    def __init__(self, info_list, step):
        super().__init__()
        self.info_list = info_list
        self.step = step
        self.type = 'sp'
        self.workspace = Path('sp_workspace')
        self.job_csv = 'sp_job.csv'
        self.state_csv = 'sp_state.csv'

    def initialize(self):
        rmf(self.job_csv)
        rmf(self.state_csv)
        mkdir(self.workspace)
        
        selected_info = []
        folder_list = getLocalFolderList()
        for i in range(len(folder_list)):
            info = self.info_list[i]
            if not info[-3]:
                r_folder = folder_list[i].joinpath('_Re_Ori').joinpath(str(1))
            else:
                r_folder = folder_list[i].joinpath('_Re_Opt_Eq').joinpath(str(info[-4]))

            if not info[-1]:
                p_folder = folder_list[i].joinpath('_Pr_Ori').joinpath(str(1))
            else:
                p_folder = folder_list[i].joinpath('_Pr_Opt_Eq').joinpath(str(info[-2]))
            spin_chrg = info[1]
            idx = info[0]
            selected_info.append([idx, spin_chrg, r_folder, p_folder])
        
        df = pd.DataFrame(selected_info)
        df.to_csv('sp.csv',index=None)
    
    def buildSPGjf(self):
        infos = csv2List('sp.csv')
        for info in infos:
            spin, chrg = ast.literal_eval(info[1])
            r_folder = Path(info[2])
            p_folder = Path(info[3])
            idx = info[0]

            r_xyz = [file for file in r_folder.glob('*.xyz')][0]
            p_xyz = [file for file in p_folder.glob('*.xyz')][0]
            r_coord = xyz2coordStr(r_xyz)
            p_coord = xyz2coordStr(p_xyz)

            r_gjf = str(idx) + '/' + '-'.join(str(r_xyz).strip('.xyz').split('/')[1:]) + '.gjf'
            p_gjf = str(idx) + '/' + '-'.join(str(p_xyz).strip('.xyz').split('/')[1:]) + '.gjf'
            G16_Set(r_gjf, job_type='sp', coord=r_coord, spin=spin, chrg=chrg).setGjf()
            G16_Set(p_gjf, job_type='sp', coord=p_coord, spin=spin, chrg=chrg).setGjf()

    def prepareQM(self, gjf_list=None, workspace=None):
        if gjf_list == None:
            folder_list = getLocalFolderList()
            gjf_list = compressList(getSuffixFileList(folder_list, 'gjf'))
        gjfIdx_list = makeIdxList(gjf_list)
        group_list = splitSublists(gjfIdx_list, self.step)

        recordJobCSV(self.type, gjf_list)
        buildJobPbs(self.comp_config["nproc"], 
            self.type, group_list, self.workspace)

        recordStateCSV(self.type, self.workspace)

        moveGjfList(gjf_list, self.workspace, self.type)

        myLogger.info('Finished initialize job')

    def precheck(self):
        localFolder_list = getLocalFolderList()
        sp_log_list = getSuffixFileList(localFolder_list, 'log', iterdir=False)
        err_list = []

        for logs in sp_log_list:
            for log in logs:
                spins = G16_Read(log).extractSpinDensity()
                if identify_abnormal_spins(spins):
                    err_list.append(log.parent)
        return set(err_list)

    def refineSP(self, err_list):
        infos = csv2List('sp.csv')
        refine_info = []
        print(err_list)
        for [idx, spin_chrg, r_folder, p_folder] in infos:
            spin, chrg = ast.literal_eval(spin_chrg)
            if Path(r_folder).parents[1] in err_list or Path(p_folder).parents[1] in err_list:
                spin += 2
                spin_chrg = (spin, chrg)
                refine_info.append([idx, spin_chrg])
        df = pd.DataFrame(refine_info)
        df.to_csv('_refine_sp.csv',index=None)

    def rebuildGjf(self):
        infos = csv2List('_refine_sp.csv')
        gjf_list = []
        for [folder_idx, spin_chrg] in infos:
            spin, chrg = ast.literal_eval(spin_chrg)
            log_list = [file for file in Path(str(folder_idx)).glob('*.log') 
                if file.exists()]
            for log in log_list:            
               gjf = str(log).replace('.log', '.gjf')
               gjf_list.append(Path(gjf))
               coord = atom_coord2str(G16_Read(log).extractCoord()[0])
               G16_Set(gjf, job_type='sp', coord=coord, spin=spin, chrg=chrg).setGjf()
        return gjf_list

    def main(self):
        # self.initialize()
        # self.buildSPGjf()
        # self.prepareQM()
        # submitJobs(self.type, self.state_csv, self.workspace)
        # job_list = csv2List(self.job_csv)
        # log_list = getIdxLogList(self.workspace)
        # moveLogList(job_list, log_list)
        # err_list = self.precheck()
        # self.refineSP(err_list)
        # gjf_list = self.rebuildGjf()
        self.type = 're_sp'
        self.workspace = Path('re_sp_workspace')
        # mkdir(self.workspace)
        # self.prepareQM(gjf_list)
        self.job_csv = 're_sp_job.csv'
        self.state_csv = 're_sp_state.csv'
        job_list = csv2List(self.job_csv)
        log_list = getIdxLogList(self.workspace)
        # submitJobs(self.type, self.state_csv, self.workspace)
        moveLogList(job_list, log_list)

def getBondfromcrd(atom_coord):
    mol = openbabel.OBMol()
    mol.BeginModify()
    for idx, (sym, position) in enumerate(atom_coord):
        a = mol.NewAtom(idx)
        a.SetAtomicNum(Chem.Atom(sym).GetAtomicNum())
        a.SetVector(*position)
    mol.ConnectTheDots()
    mol.PerceiveBondOrders() 
    mol.EndModify()

    bond = []
    for b in openbabel.OBMolBondIter(mol):
        s1 = b.GetBeginAtom().GetId()
        s2 = b.GetEndAtom().GetId()
        bond.append((s1, s2))
    return bond


def build_reaction_relations(data, file_path):
    relations = []

    reactants = {key: set(value) for key, value in data.items() if 'Re' in key}
    products = {key: set(value) for key, value in data.items() if 'Pr' in key}

    for reactant_file, reactant_atoms in reactants.items():
        for product_file, product_atoms in products.items():
            common_atoms = set(reactant_atoms).intersection(product_atoms)
            if common_atoms:
                relations.append({
                    'Reactant': reactant_file,
                    'Product': product_file
                })

    df = pd.DataFrame(relations)
    df.to_csv(file_path, index=False)

def identify_abnormal_spins(spins, threshold=0.01):
    """
    params:
    spins (list of float): spin value list。
    threshold (float): 。

    return:
    bool: if spin abnormal, return True; else False。
    """
    non_zero_spins = [abs(spin) for spin in spins if abs(spin) > threshold]
    return len(non_zero_spins) < 2

# from src._attr_mole import setChrg, setSpin
class G16_SP_OPT(SharedSetting):
    
    def __init__(self, step):
        super().__init__()
        self.pre_space = Path('pre_space')
        self.step = step
        self.type = 'opt'
        self.work_space = Path('opt_workspace')
        self.job_csv = 'opt_job.csv'
        self.state_csv = 'opt_state.csv'

    def initialize(self):

        mkdir(self.work_space)
        
        localFolder_list = getLocalFolderList()
        workspaceFolder_list = [folder.joinpath(self.pre_space) for folder in localFolder_list]
        sp_log_list = getSuffixFileList(localFolder_list, 'log', iterdir=False)
        for folder in workspaceFolder_list:
            mkdir(folder)
        
        for logs in sp_log_list:
            data = {}
            for log in logs:
                coords = G16_Read(log).extractCoord()[0]
                spins = G16_Read(log).extractSpinDensity()
                chrgs = G16_Read(log).extractAtomChrg()
                bonds = getBondfromcrd(coords)
                for idx, fragments in enumerate(getConnection(bonds)):
                    spin = setSpin(sum([spins[idx] for idx in fragments[0]]))
                    chrg = setChrg(sum([chrgs[idx] for idx in fragments[0]]))
                    coord_str = atom_coord2str([coords[idx] for idx in fragments[0]])
                    name = str(log.parent.joinpath(self.pre_space).joinpath(f'%s-%s.log' % (str(log).split('_')[1], idx)))
                    data[name] = fragments[0]
                    gjf = log.parent.joinpath(self.pre_space).joinpath(f'%s-%s.gjf' % (str(log).split('_')[1], idx))
                    G16_Set(gjf=str(gjf),job_type='opt',coord=coord_str,spin=spin, chrg=chrg).setGjf()

                    print(gjf, spin, chrg, fragments[0])
                    print(sum([spins[idx] for idx in fragments[0]]))
                    print()
                    # print(coord_str)
                
                build_reaction_relations(data, log.parent.joinpath('_relation.csv'))
                

    def prepareQM(self):       
        folder_list = extendFolderList(getLocalFolderList(),self.pre_space)
        gjf_list = compressList(getSuffixFileList(folder_list, 'gjf'))

        gjfIdx_list = makeIdxList(gjf_list)
        group_list = splitSublists(gjfIdx_list, self.step)

        recordJobCSV(self.type, gjf_list)
        buildJobPbs(self.comp_config["nproc"], 
            self.type, group_list, self.work_space)

        recordStateCSV(self.type, self.work_space)

        moveGjfList(gjf_list, self.work_space, self.type)

        myLogger.info('Finished initialize job')
    
    def visualize(self):
        folder_list = extendFolderList(getLocalFolderList(),self.pre_space)
        log_list = compressList(getSuffixFileList(folder_list, 'log'))
        log2Smiles, log2formula, log2img = {}, {}, {}
        for log in log_list:
            print(log, G16_Ropt(log).extractTerminType())
            smiles, formula = convert_gau_output_to_img_by_obabel(log)
            log_str = str(log)
            log2Smiles[log_str] = smiles
            log2formula[log_str] = formula
            log2img[log_str] = log_str.replace('.log', '.svg')
        return log2Smiles, log2formula, log2img
    
    def collect_relations(self):
        localFolder_list = getLocalFolderList()
        csv_list = extendFolderList(localFolder_list, '_relation.csv')
        dataframes = [pd.read_csv(file) for file in csv_list]
        combined_df = pd.concat(dataframes, ignore_index=True)
        combined_df.to_csv('combined_csv.csv', index=False)

    def draw_net(self, log2Smiles, log2formula, log2img):
        images = {k: PIL.Image.open(fname) for k, fname in log2img.items()}
        relations = csv2List('combined_csv.csv')
        nodes = list(set(compressList(relations)))
        print(nodes)
        # Generate the computer network graph
        G = nx.Graph()
        for node in nodes:
            G.add_node(log2Smiles[node], image=images[node])
        for relation in relations:
            G.add_edge(log2Smiles[relation[0]], log2Smiles[relation[1]])
        # Get a reproducible layout and create figure
        pos = nx.spring_layout(G)
        fig, ax = plt.subplots()
        # Note: the min_source/target_margin kwargs only work with FancyArrowPatch objects.
        # Force the use of FancyArrowPatch for edge drawing by setting `arrows=True`,
        # but suppress arrowheads with `arrowstyle="-"`
        nx.draw_networkx_edges(
            G,
            pos=pos,
            ax=ax,
            arrows=True,
            arrowstyle="-",
            min_source_margin=15,
            min_target_margin=15,
        )
        # Transform from data coordinates (scaled between xlim and ylim) to display coordinates
        tr_figure = ax.transData.transform
        # Transform from display to figure coordinates
        tr_axes = fig.transFigure.inverted().transform

        # Select the size of the image (relative to the X axis)
        icon_size = (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.025
        icon_center = icon_size / 2.0

        # Add the respective image to each node
        for n in G.nodes:
            xf, yf = tr_figure(pos[n])
            xa, ya = tr_axes((xf, yf))
            # get overlapped axes and plot icon
            a = plt.axes([xa - icon_center, ya - icon_center, icon_size, icon_size])
            a.imshow(G.nodes[n]["image"])
            a.axis("off")
        plt.show()
    
    def build_network(self, log2Smiles, log2formula, log2img):
        relations = csv2List('combined_csv.csv')
        nodes = list(set(compressList(relations)))
        print(nodes)
        G = nx.DiGraph()
        for node in nodes:
            G.add_node(log2Smiles[node], image=log2img[node])
        for relation in relations:
            G.add_edge(log2Smiles[relation[0]], log2Smiles[relation[1]])
        pos = nx.spring_layout(G)
        return G


    def main(self):
        # self.initialize()
        # self.prepareQM()
        # submitJobs(self.type, self.state_csv, self.work_space)
        # job_list = csv2List(self.job_csv)
        # log_list = getIdxLogList(self.work_space)
        # moveLogList(job_list, log_list)
        log2Smiles, log2formula, log2img = self.visualize()
        print(log2Smiles, log2formula, log2img)
        self.collect_relations()
        G = self.build_network(log2Smiles, log2formula, log2img)
        return G
        # self.draw_net(log2Smiles, log2formula, log2img)

def generate_cyto_elements(G):
    """将 NetworkX 图转换为 Cytoscape 元素"""
    return [
        {'data': {'id': n, 'label': n}} for n in G.nodes()
    ] + [
        {'data': {'source': u, 'target': v}} for u, v in G.edges()
    ]

def convert_gau_output_to_img_by_obabel(file_path):
    file_path = str(file_path)

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("g16", "sdf")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, file_path)

    pybel_mol = pybel.Molecule(mol)
    smiles = pybel_mol.write("smi").strip().split("\t")[0]
    formula = pybel_mol.formula
    # sdf_path = file_path.replace('.log', '.sdf')
    # pybel_mol.write("sdf", sdf_path, overwrite=True)

    # img_path = file_path.replace('.log', '.svg')
    # subprocess.run(["obabel", sdf_path, "-d", "-O", img_path, "-xP", "200", "-xp", "200", "-xd"])

    return smiles, formula

def convert_gau_output_to_img(file_path):
    file_path = str(file_path)
    
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("g16", "sdf")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, file_path)

    sdf_output = file_path.replace('.log', '.sdf')
    obConversion.WriteFile(mol, sdf_output)
    
    rdkit_mol = Chem.SDMolSupplier(sdf_output)[0]
    smiles = Chem.MolToSmiles(rdkit_mol)
    formula = rdMolDescriptors.CalcMolFormula(rdkit_mol)

    img_path = file_path.replace('.log', '.png')
    img = Draw.MolToImage(rdkit_mol)
    img.save(img_path)
    return smiles, formula


G = G16_SP_OPT(2).main()
app = dash.Dash(__name__)
app.layout = html.Div([
        cyto.Cytoscape(
            id='cytoscape-graph',
            elements=generate_cyto_elements(G),
            style={'width': '100%', 'height': '400px'},
            layout={'name': 'breadthfirst'}
        ),
dcc.Dropdown(
        id='node-dropdown',
        options=[{'label': n, 'value': n} for n in G.nodes()],
        value='Node 1'
    ),
    html.Button('Delete Node', id='delete-button')
])
@app.callback(
    Output('cytoscape-graph', 'elements'),
    Input('delete-button', 'n_clicks'),
    [dash.dependencies.State('node-dropdown', 'value')]
)
def delete_node(n_clicks, node_to_delete):
    """删除选定的节点"""
    if n_clicks:
        G.remove_node(node_to_delete)
        return generate_cyto_elements(G)
    return dash.no_update

if __name__ == '__main__':
    # l = csv2List('react_concise_data.csv')
    # info_list = [OPT_Job(i, idx+1).main() for idx, i in enumerate(l)]

    # G16_SP(info_list, 2).main()

    app.run_server(debug=True, port=8080)
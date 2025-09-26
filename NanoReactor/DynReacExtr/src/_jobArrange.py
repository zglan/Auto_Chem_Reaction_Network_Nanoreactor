import os
import re
import ast
import PIL

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.gridspec import GridSpec
from tqdm import tqdm

from _toolkit import mkdir, limitList
from _setting_share import SharedSetting
from _trajOper import OperTraj
from _smiles import SMILES

class Arrangejob(SharedSetting):

    @staticmethod
    def cleanJob(job_list):

        new_job_list = [] ; count = 1
        for job in job_list:
            reaction = job[0]
            reactant, product = reaction.split('->')
            try:
                # ban only H element reaction
                eval(re.sub(r'\[H\]', '1', reaction.split('->')[0]))
            except:
                # make sure reactant no same as product
                if reactant != product:
                    # ban ternary reaction
                    if len(reactant.split('+')) < 4:
                        new_job_list.append((job, count))
                        count += 1

        return new_job_list

    def setG16Job(self, job_list, job_type, step_linenum):
        
        job_list = self.cleanJob(job_list)
        
        if len(job_list) == 0:
            pass
        else:
            mkdir(job_type)
            for job in job_list:
                mkdir('%s/%s' % (job_type, job[-1]))
                t_list = np.arange(
                    -self.jobset[0], self.jobset[0] + 1, self.jobset[1]) + job[0][-1]
                t_list = limitList(t_list, step_linenum, 0)
                
                OperTraj().extractG16Job(job, t_list, job_type)

                with open('%sjob.info' % job_type, 'a') as info:
                    info.write('%s___%s___%s___%s\n' 
                        % (job[-1], job[0][0], job[0][-1], t_list.tolist()))
    
    def setXtbOptJob(self, job_list, job_type, step_linenum):
        
        job_list = self.cleanJob(job_list)
        
        if len(job_list) == 0:
            pass
        else:
            mkdir(job_type)
            src = os.getcwd()

            columns = ['Job_index', 'Index', 'Smi_equal',
                       'Xyz', 'Opt_Xyz', 'Chrg', 'Spin', 
                       'SMI', 'SMI_Opt']
            df = pd.DataFrame(columns=columns) 
            df.to_csv('_opt_xtb_data.csv', index=False)
            
            for job in job_list:
                mkdir('%s/%s' % (job_type, job[-1]))
                t_list = np.arange(
                    -self.jobset[0], self.jobset[0] + 1, self.jobset[1]) + job[0][-1]
                t_list = limitList(t_list, step_linenum, 0)
                data = OperTraj().runXtbOptJob(job, t_list, job_type)
                os.chdir(src)

                df = df.append(data)
                df.to_csv('_opt_xtb_data.csv', index=False)

    
    def setOrcaNEBJob(self, job_list, job_type, step_linenum):
        job_list = self.cleanJob(job_list)

        if len(job_list) == 0:
            pass
        else:
            mkdir(job_type)
            for job in job_list:
                mkdir('%s/%s' % (job_type, job[-1]))

                t_list = np.array(np.arange(-self.jobset[0], self.jobset[0] + 1) + job[0][-1])
                t_list = limitList(t_list, step_linenum, 0).tolist()

                OperTraj().makePathPoints(job_type, job, t_list)

                with open('%sjob.info' % job_type, 'a') as info:
                    info.write('%s___%s___%s___%s\n' 
                    % (job[-1], job[0][0], job[0][-1], t_list))
    
    def setXtbPathJob(self, job_list, job_type, step_linenum):
        job_list = self.cleanJob(job_list)
        
        if len(job_list) == 0:
            pass
        else:
            mkdir(job_type)
            for job in job_list:
                mkdir('%s/%s' % (job_type, job[-1]))

                t_list = np.array(
                    np.arange(-self.jobset[0], self.jobset[0] + 1) + job[0][-1])
                t_list = limitList(t_list, step_linenum, 0).tolist()

                OperTraj().makePathPoints(job_type, job, t_list)

                with open('%sjob.info' % job_type, 'a') as info:
                    info.write('%s___%s___%s___%s\n' 
                    % (job[-1], job[0][0], job[0][-1], t_list))
    
    def drawReacs(self, job_list, atom_list):

        ### draw SMILES pics ###
        smi_folder = self.loc_config['reaction_smi_folder']
        mkdir(smi_folder)

        png_dict = {}
        loc_list = []
        loc_max = 0
        smi_idx = 1
        for (idx, job) in enumerate(tqdm(job_list)):
            reaction = job[0]
            smi_group = [s.split('+') for s in reaction.split('->')]
            loc_idx = 0
            loc_list.append([idx, 0, -1]) ### ordered num
            loc_list.append([idx, len(smi_group[0]) + 1, 0]) ### reaction arrow
            for smi_list in smi_group:
                loc_idx += 1
                for smi in smi_list:
                    loc = smi_folder.joinpath(
                        '%s.png' % (smi_idx))
                    SMILES.drawSMILES(
                        SMILES.rewriteSMILES(smi, atom_list), loc)
                    loc_list.append([idx, loc_idx, smi_idx])
                    png_dict[smi_idx] = loc
                    smi_idx += 1
                    loc_idx += 1
            loc_max = loc_idx if loc_idx > loc_max else loc_max
        imgs = {k: PIL.Image.open(fname) for k, fname in png_dict.items()}
        
        ### draw reaction pics ###
        ax = plt.gca()
        ax.set_axis_off()

        n = len(job_list)
        fig = plt.figure(figsize=(2 * loc_max, 2 * n))
        gs = GridSpec(n, loc_max + 1, figure=fig)

        for [idx, loc_idx, png_idx] in tqdm(loc_list):
            str_idx = 'ax%s' % (idx)
            locals()[str_idx] = fig.add_subplot(
                gs[idx, loc_idx])
            if loc_idx == 0:
                locals()[str_idx].text(
                  0.5, 0.5, str(idx),
                  horizontalalignment='center', 
                  verticalalignment='center')
            elif png_idx == 0:
                locals()[str_idx].text(
                  0.5, 0.5, '==>',
                  horizontalalignment='center', 
                  verticalalignment='center')
            else:
                locals()[str_idx].imshow(imgs[png_idx])
            locals()[str_idx].axis('off')

        fig.tight_layout()
        plt.savefig(self.loc_config['reaction_pic'])

    def drawReacNet(self, atom_list):

        def loadRealat():
            try:
                if os.path.isfile(self.loc_config['relation_file']):
                    if os.path.isfile(self.loc_config['ref_relation_file']):
                        df = pd.read_csv(self.loc_config['ref_relation_file'])
                        ref_relat = True
                    else:
                        df = pd.read_csv(self.loc_config['relation_file'])
                        ref_relat = False
                else:
                    raise FileNotFoundError('No species relation file exists')
            except FileNotFoundError as err:
                print(err)
            relats = [ast.literal_eval(i) for i in df['relation'].to_list()]
            return relats, ref_relat

        def getInfo(relats):
            n = 1
            smi_dict = {}
            idx_relat_list, smi_list = [], []

            for relat_list in relats:
                for relat in relat_list:
                    for smi in relat:
                        if smi not in smi_dict.keys():
                            smi_dict[smi] = n
                            smi_list.append(smi)
                            n += 1
                    # self loop remove
                    if smi_dict[relat[0]] != smi_dict[relat[1]]:
                        idx_relat_list.append((smi_dict[relat[0]], smi_dict[relat[1]]))
            return smi_dict, smi_list, idx_relat_list
        
        def drawNet(smi_dict, smi_list, idx_relat_list):
            from _network import ReacEveNet, upadteParam, default_net_param
            n_config = upadteParam(default_net_param, 
                {'fig_name': '_network','spec_label':{}})
            net = ReacEveNet(smi_dict, smi_list, idx_relat_list, n_config)
            net.drawNet()

        def drawMols(smi_dict, ref_relat):
            mkdir(self.loc_config['network_smi_folder'])
            png_dict = {}
            for smi in smi_dict.keys():
                idx = smi_dict[smi]
                loc = self.loc_config['network_smi_folder'].joinpath(
                    '%s.png' % (idx))
                png_dict[idx] = loc
                SMILES.drawSMILES(
                    SMILES.rewriteSMILES(smi, atom_list), loc, ref_relat)
            png_list = list(
                zip(list(png_dict.keys()), list(png_dict.values())))
            return png_dict, png_list

        def mergePic(png_dict, png_list):
            imgs = {k: PIL.Image.open(fname) 
                for k, fname in png_dict.items()}
            net_img = PIL.Image.open(self.loc_config['network_pic'])

            step = 7
            group = [
                png_list[i: i + step] 
                for i in range(0, len(png_list), step)]
            n = len(group) + step + 1

            fig = plt.figure(constrained_layout=True)
            gs = GridSpec(n, 7, figure=fig)

            ax = fig.add_subplot(gs[0:step, :])
            ax.set_aspect('auto')
            ax.axis('off')
            ax.imshow(net_img)

            for idx, g in enumerate(group):
                for i in g:
                    str_idx = 'ax%s' % (idx)               
                    locals()[str_idx] = fig.add_subplot(
                        gs[idx + step, i[0] - idx * step - 1])
                    locals()[str_idx].axis('off')
                    locals()[str_idx].set_title(i[0], fontsize=4)
                    locals()[str_idx].imshow(imgs[i[0]])
            plt.savefig(
                self.loc_config['reaction_network_pic'], dpi=800)

        relats, ref_relat = loadRealat()
        smi_dict, smi_list, idx_relat_list = getInfo(relats)
        drawNet(smi_dict, smi_list, idx_relat_list)
        png_dict, png_list = drawMols(smi_dict, ref_relat)
        mergePic(png_dict, png_list)

    def setRefineReaction(self, job_list):
        job_list = self.cleanJob(job_list)
        relation_list = []
        
        if len(job_list) == 0:
            pass
        else:
            for job in job_list:
                t_list = job[0][-2]
                relation = OperTraj().refReaction(job, t_list)
                relation_list.append((job[0][0], relation))
        
        columns = ['reaction', 'relation']
        df = pd.DataFrame(columns=columns, data=relation_list)
        df.to_csv(self.loc_config['ref_relation_file'])
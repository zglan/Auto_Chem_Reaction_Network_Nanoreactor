#!/usr/bin/env python3

"""
DynReacExtr: Automatic reaction event extraction tool and
  reaction assistant analysis tool
  for reactive molecular dynamics simulation.

============================
Features
============================
  - Analyze reaction events
  - Draw reaction network
  - Extract TS guess
============================
Simple example
============================
DynReacExtr can process XYZ format trajectory files
$ DynReacExtr.py -i traj.xyz

You can running the following script for help:
$ DynReacExtr.py -h
"""

from _setting_share import SharedSetting
from _trajOper import OperTraj
from _hmm import FilterNoise
from _reacFind import FindReaction
from _jobArrange import Arrangejob


class DynReacExtr(SharedSetting):
    
    def operateTraj(self):
        if self.load:
            import numpy as np
            TrajInfo = np.load(
                self.loc_config['ini_read_data'], allow_pickle=True)
            self.atom_list, self.step_linenum, \
             self.bond_data,  self.step, self.smi_nodes = TrajInfo
        else:
            TrajInfo = OperTraj()
            self.atom_list = TrajInfo.atom_list
            self.step_linenum = TrajInfo.step_linenum
            self.step = TrajInfo.step
            self.bond_data = TrajInfo.bond_data
            self.smi_nodes = TrajInfo.smi_nodes
        
    def getTimeNodeDict(self):
        if self.load:
            pass
        else:
            Filter = FilterNoise()
            r_list_t_node, p_list_t_node = Filter.getFilteredSig(
                self.smi_nodes, self.step)
            self.r_dict_t_node = Filter.uniteSig(r_list_t_node)
            self.p_dict_t_node = Filter.uniteSig(p_list_t_node)
    
    def findReactions(self):
        if self.load:
            import json
            import pandas as pd
            self.job_list = pd.read_csv('_job.csv').values.tolist()
            for idx, job in enumerate(self.job_list):
                self.job_list[idx] = [
                    job[1], json.loads(job[2]), 
                    tuple(json.loads(
                        job[3].replace('(', '[').replace(')', ']'))), job[4]]
        else:
            self.job_list =FindReaction().findReaction(
                self.r_dict_t_node, self.p_dict_t_node, 
                self.atom_list, self.bond_data)

    def arrangeJobs(self):
        if self.ts:
            Arrangejob().setG16Job(self.job_list, 'ts', self.step)
        if self.opt:
            Arrangejob().setG16Job(self.job_list, 'opt', self.step)
        if self.neb:
            Arrangejob().setOrcaNEBJob(self.job_list, 'neb', self.step)
        if self.xtbpath:
            Arrangejob().setXtbPathJob(self.job_list, 'path', self.step)
        if self.xtbopt:
            Arrangejob().setXtbOptJob(self.job_list, 'xtbopt', self.step)
        if self.refine:
            Arrangejob().setRefineReaction(self.job_list)
        if self.draw:
            Arrangejob().drawReacs(self.job_list, self.atom_list)
        if self.network:
            Arrangejob().drawReacNet(self.atom_list)

    def main(self):
        self.operateTraj()
        self.getTimeNodeDict()
        self.findReactions()
        self.arrangeJobs()



if __name__ == '__main__':
    DynReacExtr().main()
    
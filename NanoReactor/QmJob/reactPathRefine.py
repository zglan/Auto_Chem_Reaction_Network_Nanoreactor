#!/usr/bin/env python3

import os

from pathlib import Path

import numpy as np

from tqdm import tqdm
from sklearn.cluster import AgglomerativeClustering

from src._toolkit import *
from src._logger import MyLog
from src._api_orca import Orca_Set, Orca_Read
from src._api_gau import G16_Rts, G16_Ropt, G16_Rirc, G16_Set

from src._setting_share import SharedSetting
from src._set_pbs import buildJobPbs
from src._submit import submitJobs
from src._csv import recordStateCSV, recordJobCSV, csv2List


myLogger = MyLog().logger


#========================================================#
#                  Clustering   functions                #
#========================================================#

def cluster(l, threshold):
    if len(l) == 1:
        labels = [0]
    else:
        l = np.array(l).reshape(-1, 1)
        clustered = AgglomerativeClustering(
            n_clusters=None, distance_threshold=threshold).fit(l)
        labels = clustered.labels_
    return labels


def clusterEnerFreq_G16():
    folder_list = getLocalFolderList()
    log_list = getSuffixFileList(folder_list, 'log')
    folderlog_list = packList(folder_list, log_list)

    for folder, logList in folderlog_list:
        eleE_list = []
        freq_list = []
        logPath_list = []

        workSpace = Path(folder.joinpath('workspace'))
        mkdir(workSpace)
        
        for log in tqdm(logList, desc=f'  Loading reaction-{folder}'):
            freq, freqError = G16_Rts(log).extractFreq()
            energy = G16_Rts(log).extractEnergy()
            if freqError == False:
                eleE_list.append(energy)
                freq_list.append(freq[0])
                logPath_list.append(log)
        myLogger.info('  Bulit folder-%s workspace. Got all frequency and energy' % (folder))

        if len(eleE_list) > 0:
            labelList = cluster(eleE_list, 3)
            ungrouped = packList(labelList, eleE_list, freq_list, logPath_list)
            groups = groupByIdx(ungrouped, 0)

            for idx, group in enumerate(groups):
                mkdir(workSpace.joinpath('Ener-' + str(idx)))
                
                newFreqList = getIdxElemList(group, 2)
                newLogPathList = getIdxElemList(group, 3)
                newLabelList = cluster(newFreqList, 100)
                newungrouped = packList(newLabelList, newLogPathList)
                newgroups = groupByIdx(newungrouped, 0)

                for i, g in enumerate(newgroups):
                    mkdir(workSpace.joinpath('Ener-' + str(idx)).joinpath('Freq-' + str(i)))
                    
                    newlogPaths = getIdxElemList(g, 1)
                    for log in newlogPaths:

                        copyFile(log, workSpace
                        .joinpath('Ener-' + str(idx))
                        .joinpath('Freq-' + str(i))
                        .joinpath((str(log).split('/')[1] + '.log')))

            myLogger.info('Finished clustering. Finished folder moving job')
        else:
            myLogger.info('No need to cluster')


def clusterEnerFreq_G16_ref():
    folder_list = getSelectedFolderList(Path('ts_ref'))
    log_list = getSuffixFileList(folder_list, 'log')
    folderlog_list = packList(folder_list, log_list)

    for folder, logList in folderlog_list:
        eleE_list = []
        freq_list = []
        logPath_list = []

        workSpace = Path(folder.joinpath('workspace'))
        mkdir(workSpace)
        
        for log in tqdm(logList, desc=f'  Loading reaction-{folder}'):
            freq, freqError = G16_Rts(log).extractFreq()
            energy = G16_Rts(log).extractEnergy()
            if freqError == False:
                eleE_list.append(energy)
                freq_list.append(freq[0])
                logPath_list.append(log)
        myLogger.info('  Bulit folder-%s workspace. Got all frequency and energy' % (folder))

        if len(eleE_list) > 0:
            labelList = cluster(eleE_list, 3)
            ungrouped = packList(labelList, eleE_list, freq_list, logPath_list)
            groups = groupByIdx(ungrouped, 0)

            for idx, group in enumerate(groups):
                mkdir(workSpace.joinpath('Ener-' + str(idx)))
                
                newFreqList = getIdxElemList(group, 2)
                newLogPathList = getIdxElemList(group, 3)
                newLabelList = cluster(newFreqList, 100)
                newungrouped = packList(newLabelList, newLogPathList)
                newgroups = groupByIdx(newungrouped, 0)

                for i, g in enumerate(newgroups):
                    mkdir(workSpace.joinpath('Ener-' + str(idx)).joinpath('Freq-' + str(i)))
                    
                    newlogPaths = getIdxElemList(g, 1)
                    for log in newlogPaths:

                        copyFile(log, workSpace
                        .joinpath('Ener-' + str(idx))
                        .joinpath('Freq-' + str(i))
                        .joinpath((str(log).split('/')[1] + '.log')))

            myLogger.info('Finished clustering. Finished folder moving job')
        else:
            myLogger.info('No need to cluster')


def clusterEnerFreq_Orca():
    folder_list = getLocalFolderList()
    out_list = getSuffixFileList(folder_list, 'out')
    folderout_list = packList(folder_list, out_list)
    
    for folder, out_list in folderout_list:
        workSpace = Path(folder.joinpath('workspace'))
        mkdir(workSpace)

        eleE_list = []
        freq_list = []
        path_list = []

        for out in out_list:
            folder = out.parent
            eleE, freq, error = Orca_Read(folder, 'ts', 'ts').readTsOut()
            
            if error:
                eleE_list.append(eleE)
                freq_list.append(freq)
                path_list.append(folder)

        myLogger.info('Bulit folder-%s workspace. Got all frequency and energy' % (folder))

        if len(eleE_list) > 0:
            label_list = cluster(eleE_list, 5)
            ungrouped = packList(label_list, eleE_list, freq_list, path_list)
            grouped = groupByIdx(ungrouped, 0)

            for idx, group in enumerate(grouped):
                mkdir(workSpace.joinpath('Ener-' + str(idx)))
                
                newFreq_list = getIdxElemList(group, 2)
                newPath_list = getIdxElemList(group, 3)
                newLabel_list = cluster(newFreq_list, 100)
                newungrouped = packList(newLabel_list, newPath_list)
                newgrouped = groupByIdx(newungrouped, 0)

                for i, g in enumerate(newgrouped):
                    mkdir(workSpace.joinpath('Ener-' + str(idx)).joinpath('Freq-' + str(i)))
                    
                    newPaths = getIdxElemList(1, g)
                    
                    for ip, path in enumerate(newPaths):
                        coverFolder(path, 
                            workSpace.joinpath('Ener-' + str(idx))
                                .joinpath('Freq-' + str(i))
                                .joinpath(str(ip))
                                )
            
            myLogger.info('Finished clustering. Finished folder moving job')
        
        else:
            myLogger.info('No need to cluster')

##########################################

vmdTcl = '''
mol load xyz %s.xyz
display depthcue off
display projection Orthographic
display rendermode GLSL
display height 12
%s
render snapshot %s.png
quit'''

def xyz2array(xyz):
    arr = []
    for l in xyz:
        s = l.split()
        arr.append([float(x) for x in s[1:4]])
    return np.array(arr)

def getDirection(arr):
    l = []
    for i in range(0, len(arr)):
        for j in range(i + 1, len(arr)):
            c = np.abs(arr[i] - arr[j])
            l.append(c / np.linalg.norm(c))
    return l

def checkDirection(d):
    if len([i for i in d 
        if any(i == np.array([0,0,1]))
            ]) == len(d):
                return True
    else:
        return False

#############################################

class Orca_Neb:
    
    def __init__(self, num, step):
        self.num = num
        self.step = step
        self.type = 'orcaNeb'
        self.workspace = Path('orcaNeb_workspace')
        self.job_csv = 'orcaNeb_job.csv'
        self.state_csv = 'orcaNeb_state.csv'

    def initialize(self):
        mkdir(self.workspace)
        rmf(self.job_csv)
        rmf(self.state_csv)
        
        folder_list = getLocalFolderList()
        neb_list = compressList(
            [getSubFolderList(folder)[:self.num] for folder in folder_list])
        nebIdx_list = makeIdxList(neb_list)
        group_list = splitSublists(nebIdx_list, self.step)

        buildJobPbs(8, self.type, group_list, self.workspace)
        recordStateCSV(self.type, self.workspace)
        moveFolderList(neb_list, self.workspace)

        myLogger.info('Finished initialize job')
        
    def main(self):
        self.initialize()


class Orca_TS:

    def __init__(self, step):
        self.step = step
        self.type = 'orcaTs'
        self.workspace = Path('orcaTs_workspace')
        self.job_csv = 'orcaTs_job.csv'
        self.state_csv = 'orcaTs_state.csv'

    def initialize(self):
        mkdir(self.workspace)
        rmf(self.job_csv)
        rmf(self.state_csv)

        folder_list = getLocalFolderList()
        ts_list = compressList(getNameFileList(folder_list, 'xtbpath_ts.xyz'))
        folder_list = [ts.parent for ts in ts_list]
        chrgspin_list = self.getChrgSpin(folder_list)
        folderIdx_list = makeIdxList(folder_list)
        group_list = splitSublists(folderIdx_list, self.step)

        self.makeTsInp(chrgspin_list, folder_list)

        buildJobPbs(8, self.type, group_list, self.workspace)
        recordJobCSV(self.type, folder_list)
        recordStateCSV(self.type, self.workspace)
        moveFolderList(folder_list, self.workspace)

        myLogger.info('Finished initialize job')

    @staticmethod
    def makeTsInp(chrgspin_list, folder_list):
        for (ch, sp), folder in packList(chrgspin_list, folder_list):
            Orca_Set(folder.joinpath('ts.inp'), 'ts', 'xtbpath_ts.xyz', chrg=ch, spin=sp)

    @staticmethod
    def getChrgSpin(folder_list):
        chrgspin_list = []
        for folder in folder_list:
            with open(Path(folder).joinpath('path.inp')) as inp:
                lines = inp.readlines()
                for l in lines:
                    if 'spin' in l:
                        spin = int(l.split()[-1]) + 1
                    if 'chrg' in l:
                        chrg = int(l.split()[-1])
                chrgspin_list.append((chrg, spin))
        return chrgspin_list
    
    def moveTsFolder(self, folder_list, job_list):
        for folder in folder_list:
            for idx, path, state in job_list:
                if folder == Path('orcaTs_workspace/%s' % (idx)):
                    coverFolder(folder, path)
        
        myLogger.info('Finished Ts folder moving job')


    def main(self):
        self.initialize()
        submitJobs(self.type, self.state_csv, self.workspace)
        job_list = csv2List(self.job_csv).values.tolist()
        folder_list = getSubFolderList(self.workspace)
        self.moveTsFolder(folder_list, job_list)
        clusterEnerFreq_Orca()


class Orca_Irc:

    def __init__(self, step):
        self.step = step
        self.type = 'orcaIrc'
        self.tmp = Path('tmp')
        self.irc = Path('orcaIrc')
        self.workspace = Path('orcaIrc_workspace')
        self.job_csv = 'orcaIrc_job.csv'
        self.state_csv = 'orcaIrc_state.csv'

    def initialize(self):
        mkdir(self.tmp)
        mkdir(self.irc)
        mkdir(self.workspace)

        localFolder_list = getLocalFolderList()
        workspaceFolder_list = [folder.joinpath('workspace') for folder in localFolder_list]

        for folder in workspaceFolder_list:
            copyFolder(folder, self.tmp.joinpath(folder.parts[0]))

        tmpFolder_list = getSelectedFolderList(self.tmp)
        tmp_xyz_list = getNameFileList(tmpFolder_list, 'ts.xyz')

        sel_list = []
        for xyz_list in tmp_xyz_list:
            if len(xyz_list) > 0:
                d = {}
                for xyz in xyz_list:
                    if str(xyz.parents[1]) not in d.keys():
                        d[str(xyz.parents[1])] = None
                        with open(xyz.parent.joinpath('ts.inp')) as inp:
                            lines = inp.readlines()
                            for line in lines:
                                if '* XYZfile' in line:
                                    chrg, spin = line.split()[2], line.split()[3]
                        sel_list.append((xyz, chrg, spin))
        
        move_list = []

        for xyz, chrg, spin in sel_list:
            inp = xyz.parent.joinpath('irc.inp')
            Orca_Set(inp, 'irc', 'ts.xyz', chrg, spin)
            move_list.append((xyz, inp))

        for idx, (xyz, inp) in enumerate(move_list):
            drc = self.workspace.joinpath(str(idx + 1))
            mkdir(drc)
            copyFile(xyz, drc)
            copyFile(inp, drc)

        Idx_list = makeIdxList(move_list)
        group_list = splitSublists(Idx_list, self.step)

        recordJobCSV(self.type, move_list)
        buildJobPbs(8, self.type, group_list, self.workspace)
        recordStateCSV(self.type, self.workspace)
    
    def moveIrc(self):
        ircFolder_list = getSelectedFolderList(self.workspace)

        for ircFolder in ircFolder_list:
            ePath, error = Orca_Read(ircFolder, 'irc', 'irc').readIrcOut()
            if error:
                eGap = abs(ePath[0] - ePath[-1])
                if eGap >= 2:
                    copyFolder(ircFolder, self.irc.joinpath(str(ircFolder).split('/')[-1]))
        
        myLogger.info('Finished moveing irc log job')

    def main(self):
        # self.initialize()
        # submitAll(self.type, self.state_csv, self.workspace)
        self.moveIrc()

#############################################

class Xtb_Path:
    
    def __init__(self, num, step):
        self.num = num
        self.step = step
        self.type = 'xtbPath'
        self.workspace = Path('xtbPath_workspace')
        self.job_csv = 'xtbPath_job.csv'
        self.state_csv = 'xtbPath_state.csv'

    def initialize(self):
        mkdir(self.workspace)
        rmf(self.job_csv)
        rmf(self.state_csv)
        
        folder_list = getLocalFolderList()
        path_list = compressList(
            [getSubFolderList(folder)[:self.num] for folder in folder_list])
        pathIdx_list = makeIdxList(path_list)
        group_list = splitSublists(pathIdx_list, self.step)

        buildJobPbs(4, self.type, group_list, self.workspace)
        recordJobCSV(self.type, path_list)
        recordStateCSV(self.type, self.workspace)
        moveFolderList(path_list, self.workspace)

        myLogger.info('Finished initialize job')

    def moveTsXyz(self, ts_list, job_list):
        for ts in ts_list:
            for idx, path, state in job_list:
                if ts.parent == Path('xtbPath_workspace/%s' % (idx)):
                    copyFile(ts, path)

    def main(self):
        self.initialize()
        submitJobs(self.type, self.state_csv, self.workspace)
        job_list = csv2List(self.job_csv)
        ts_list = extendFolderList(getSubFolderList(self.workspace), 'xtbpath_ts.xyz')
        self.moveTsXyz(ts_list, job_list)

#############################################

class G16_TS(SharedSetting):
    
    def __init__(self, num, step):
        super().__init__()
        self.num = num
        self.step = step
        self.type = 'ts'
        self.work_space = Path('ts_workspace')
        self.job_csv = 'ts_job.csv'
        self.state_csv = 'ts_state.csv'

    def initialize(self):
        mkdir(self.work_space)
        rmf(self.job_csv)
        rmf(self.state_csv)
        
        folder_list = getLocalFolderList()
        gjf_list = getSuffixNumsFileList(folder_list, 'gjf', self.num)
        gjfIdx_list = makeIdxList(gjf_list)
        group_list = splitSublists(gjfIdx_list, self.step)

        recordJobCSV(self.type, gjf_list)
        buildJobPbs(self.server_config["nproc"], self.type, group_list, self.work_space)

        recordStateCSV(self.type, self.work_space)

        moveGjfList(gjf_list, self.work_space, self.type)

        myLogger.info('Finished initialize job')
        
    def main(self):
        self.initialize()
        # submitJobs(self.type, self.state_csv, self.work_space)
        # job_list = csv2List(self.job_csv)
        # log_list = getIdxLogList(self.work_space)
        # moveLogList(job_list, log_list)
        # clusterEnerFreq_G16()


class G16_TS_Ref(SharedSetting):
    
    def __init__(self, step):
        super().__init__()
        self.step = step
        self.type = 'ts_ref'
        self.ts_ref = Path('ts_ref')
        self.ts_tmp = Path('ts_tmp')
        self.work_space = Path('ts_ref_workspace')
        self.job_csv = 'ts_ref_job.csv'
        self.state_csv = 'ts_ref_state.csv'
    
    def initialize(self):
        mkdir(self.ts_tmp)
        mkdir(self.work_space)

        localFolder_list = getLocalFolderList()
        workspaceFolder_list = [folder.joinpath('workspace') for folder in localFolder_list]

        for folder in workspaceFolder_list:
            copyFolder(folder, self.ts_tmp.joinpath(folder.parts[0]))

        tmpFolder_list = getSelectedFolderList(self.ts_tmp)
        tmp_TsLog_list = getSuffixFileList(tmpFolder_list, 'log')

        selectedLog_list = []
        for tsLog_list in tmp_TsLog_list:
            if len(tsLog_list) > 0:
                d = {}
                for tsLog in tsLog_list:
                    if str(tsLog.parent) not in d.keys():
                        d[str(tsLog.parent)] = None
                        selectedLog_list.append(tsLog)

        for log in selectedLog_list:
            G16_Rts(str(log)).tsRef()

        tmp_TsGif_list = getSuffixFileList(tmpFolder_list, 'gjf')
        for TsGif_list in tmp_TsGif_list:
            for tsGif in TsGif_list:
                drc = self.ts_ref.joinpath('/'.join(tsGif.parts[1:-1]))
                mkdir(drc)
                copyFile(tsGif, drc)
        
        tsfolder_list = getSelectedFolderList(self.ts_ref)
        tsgjf_list = compressList(getSuffixFileList(tsfolder_list, 'gjf'))
        gjfIdx_list = makeIdxList(tsgjf_list)
        group_list = splitSublists(gjfIdx_list, self.step)

        recordJobCSV(self.type, tsgjf_list)
        buildJobPbs(self.server_config["nproc"], self.type, group_list, self.work_space)
        recordStateCSV(self.type, self.work_space)
        moveGjfList(tsgjf_list, self.work_space, self.type)

        myLogger.info('Finished initialize job')


    def main(self):
        self.initialize()
        submitJobs(self.type, self.state_csv, self.work_space)
        job_list = csv2List(self.job_csv)
        log_list = getIdxLogList(self.work_space)
        moveLogList(job_list, log_list)
        exit()
        clusterEnerFreq_G16()


class G16_xtbTS:
    
    def __init__(self, step):
        self.step = step
        self.type = 'ts'
        self.workspace = Path('ts_workspace')
        self.job_csv = 'ts_job.csv'
        self.state_csv = 'ts_state.csv'
    
    def initialize(self):
        mkdir(self.workspace)
        rmf(self.job_csv)
        rmf(self.state_csv)

        folder_list = getLocalFolderList()
        ts_list = compressList(getNameFileList(folder_list, 'xtbpath_ts.xyz'))
        folder_list = [ts.parent for ts in ts_list]
        chrgspin_list = self.getChrgSpin(folder_list)

        self.makeTSGjf(ts_list, chrgspin_list, folder_list)

        gjf_list = compressList(getSuffixFileList(folder_list, 'gjf'))
        gjfIdx_list = makeIdxList(gjf_list)
        group_list = splitSublists(gjfIdx_list, self.step)

        recordJobCSV(self.type, gjf_list)
        buildJobPbs(8, self.type, group_list, self.workspace)
        recordStateCSV(self.type, self.workspace)
        moveGjfList(gjf_list, self.workspace, self.type)

        myLogger.info('Finished initialize job')

    @staticmethod
    def makeTSGjf(ts_list, chrgspin_list, folder_list):
        coord_list = [xyz2coord(ts) for ts in ts_list]
        for idx, ((ch, sp), folder) in enumerate(packList(chrgspin_list, folder_list)):
            G16_Set(
                jobname=str(folder.joinpath('ts')), 
                nproc=4,
                jobtype='ts',
                charge=ch,
                multiplicity=sp,
                coord=coord_list[idx]).setJob()

    @staticmethod
    def getChrgSpin(folder_list):
        chrgspin_list = []
        for folder in folder_list:
            with open(Path(folder).joinpath('path.inp')) as inp:
                lines = inp.readlines()
                for l in lines:
                    if 'spin' in l:
                        spin = int(l.split()[-1]) + 1
                    if 'chrg' in l:
                        chrg = int(l.split()[-1])
                chrgspin_list.append((chrg, spin))
        return chrgspin_list
    
    def moveTSFolder(self, folder_list, job_list):
        for folder in folder_list:
            for idx, path, state in job_list:
                if folder == Path('ts_workspace/%s' % (idx)):
                    coverFolder(folder, path)
        
        myLogger.info('Finished TS folder moving job')
    
    def main(self):
        self.initialize()
        submitJobs(self.type, self.state_csv, self.workspace)
        job_list = csv2List(self.job_csv)
        log_list = getIdxLogList(self.workspace)
        moveLogList(job_list, log_list)
        clusterEnerFreq_G16()


class G16_Irc(SharedSetting):
    
    def __init__(self, step):
        super().__init__()
        self.step = step
        self.type = 'irc'
        self.irc_tmp = Path('irc_tmp')
        self.irc = Path('irc')
        self.work_space = Path('irc_workspace')
        self.job_csv = 'irc_job.csv'
        self.state_csv = 'irc_state.csv'

    def initialize(self):
        mkdir(self.irc_tmp)
        mkdir(self.irc)
        mkdir(self.work_space)
        
        localFolder_list = getLocalFolderList()
        workspaceFolder_list = [folder.joinpath('workspace') for folder in localFolder_list]

        for folder in workspaceFolder_list:
            copyFolder(folder, self.irc_tmp.joinpath(folder.parts[0]))

        tmpFolder_list = getSelectedFolderList(self.irc_tmp)
        tmp_TsLog_list = getSuffixFileList(tmpFolder_list, 'log')

        selectedLog_list = []
        for tsLog_list in tmp_TsLog_list:
            if len(tsLog_list) > 0:
                d = {}
                for tsLog in tsLog_list:
                    if str(tsLog.parent) not in d.keys():
                        d[str(tsLog.parent)] = None
                        selectedLog_list.append(tsLog)

        for log in selectedLog_list:
            G16_Rts(str(log)).ts2irc()
    
        tmp_IrcGif_list = getSuffixFileList(tmpFolder_list, 'gjf')
        for ircGif_list in tmp_IrcGif_list:
            for ircGif in ircGif_list:
                drc = self.irc.joinpath('/'.join(ircGif.parts[1:-1]))
                mkdir(drc)
                copyFile(ircGif, drc)
        
        ircfolder_list = getSelectedFolderList(self.irc)
        ircgjf_list = compressList(getSuffixFileList(ircfolder_list, 'gjf'))
        gjfIdx_list = makeIdxList(ircgjf_list)
        group_list = splitSublists(gjfIdx_list, self.step)

        recordJobCSV(self.type, ircgjf_list)
        buildJobPbs(self.server_config["nproc"], self.type, group_list, self.work_space)
        recordStateCSV(self.type, self.work_space)
        moveGjfList(ircgjf_list, self.work_space, self.type)

        myLogger.info('Finished initialize job')
    
    def initialize_ref(self):
        mkdir(self.irc_tmp)
        mkdir(self.irc)
        mkdir(self.work_space)

        folder_list = getSelectedFolderList(Path('ts_ref'))
        workspaceFolder_list = [folder.joinpath('workspace') for folder in folder_list]
        
        for folder in workspaceFolder_list:
            # print(self.irc_tmp.joinpath(folder.parts[1]))
            copyFolder(folder, self.irc_tmp.joinpath(folder.parts[1]))

        tmpFolder_list = getSelectedFolderList(self.irc_tmp)
        tmp_TsLog_list = getSuffixFileList(tmpFolder_list, 'log')

        selectedLog_list = []
        for tsLog_list in tmp_TsLog_list:
            if len(tsLog_list) > 0:
                d = {}
                for tsLog in tsLog_list:
                    if str(tsLog.parent) not in d.keys():
                        d[str(tsLog.parent)] = None
                        selectedLog_list.append(tsLog)

        for log in selectedLog_list:
            G16_Rts(str(log)).ts2irc()
    
        tmp_IrcGif_list = getSuffixFileList(tmpFolder_list, 'gjf')
        for ircGif_list in tmp_IrcGif_list:
            for ircGif in ircGif_list:
                drc = self.irc.joinpath('/'.join(ircGif.parts[1:-1]))
                mkdir(drc)
                copyFile(ircGif, drc)
        
        ircfolder_list = getSelectedFolderList(self.irc)
        ircgjf_list = compressList(getSuffixFileList(ircfolder_list, 'gjf'))
        gjfIdx_list = makeIdxList(ircgjf_list)
        group_list = splitSublists(gjfIdx_list, self.step)

        recordJobCSV(self.type, ircgjf_list)
        buildJobPbs(self.server_config["nproc"], self.type, group_list, self.work_space)
        recordStateCSV(self.type, self.work_space)
        moveGjfList(ircgjf_list, self.work_space, self.type)

        myLogger.info('Finished initialize job')

    def moveIrc(self):
        ircFolder_list = getSelectedFolderList(self.irc)
        irc_IrcLog_list = getSuffixFileList(ircFolder_list, 'log')
        for ircLog_list in irc_IrcLog_list:
            energyGap_list = []
            for ircLog in ircLog_list:
                energyGap, energyPath, coordPath = G16_Rirc(ircLog).getPath()
                energyGap_list.append(energyGap)
            for idx, energyGap in enumerate(energyGap_list):
                if energyGap != None and energyGap >= 2:
                    part = ircLog_list[idx].parts
                    copyFile(ircLog_list[idx], Path(part[1]).joinpath(part[-1]))
        
        myLogger.info('Finished moveing irc log job')

    
    def main(self):
        # self.initialize_ref()
        submitJobs(self.type, self.state_csv, self.work_space)
        job_list = csv2List(self.job_csv)
        log_list = getIdxLogList(self.work_space)
        moveLogList(job_list, log_list)
        # self.moveIrc()

class G16_Opt:

    def __init__(self):
        self.type = 'opt'
        self.opt = Path('opt')
        self.workspace = Path('opt_workspace')
        self.job_csv = 'opt_job.csv'
        self.state_csv = 'opt_state.csv'
    
    def initialize(self):
        mkdir(self.workspace)
        mkdir(self.opt)
        rmf(self.job_csv)
        rmf(self.state_csv)
        
        folder_list = getLocalFolderList()
        log_list = [list(folder.glob('*-irc.log')) 
            for folder in folder_list 
                if len(list(folder.glob('*.log'))) > 0]
        
        for l in log_list:
            for idx, log in enumerate(l):
                dst = self.opt.joinpath(log.parts[0]).joinpath(str(idx + 1)).joinpath(log.parts[-1])
                mkdir(dst)
                copyFile(log, dst)
                G16_Rirc(dst).irc2opt()
        
        gjf_list = compressList(
            getSuffixFileList(getSelectedFolderList(self.opt), 'gjf'))
        gjfIdx_list = makeIdxList(gjf_list)
        group_list = splitSublists(gjfIdx_list, 8)

        recordJobCSV(self.type, gjf_list)
        buildJobPbs(4, self.type, group_list, self.workspace)
        recordStateCSV(self.type, self.workspace)
    
    def sperateMol(self):
        folder_list = getSelectedFolderList(self.opt)
        log_list = compressList(
            [list(folder.glob('**/*-opt.log')) 
                for folder in folder_list ])
        
        for log in log_list:
            G16_Ropt(log).main()

    def drawPng(self):
        folder_list = getSelectedFolderList(self.opt)
        xyz_list = compressList(
            [list(folder.glob('**/*.xyz')) 
                for folder in folder_list ])
        
        for xyz in xyz_list:
            name = str(xyz).split('.')[0]
            xyz_l = loadFile2List(xyz)[2:]

            if len(xyz_l) == 2:
                writeFile('vmd.tcl', 'w' ,vmdTcl % (name, 'rotate x by 90', name))
            elif len(xyz_l) <= 4:
                if checkDirection(getDirection(xyz2array(xyz_l))):
                    writeFile('vmd.tcl', 'w' ,vmdTcl % (name, 'rotate x by 90', name))
                else:
                    writeFile('vmd.tcl', 'w' ,vmdTcl % (name, '', name))            
            else:
                writeFile('vmd.tcl', 'w' ,vmdTcl % (name, '', name))

            os.system('vmd -e vmd.tcl')

    def main(self):
        self.initialize()
        submitJobs(self.type, self.state_csv, self.workspace)
        job_list = csv2List(self.job_csv)
        log_list = getIdxLogList(self.workspace)
        moveLogList(job_list, log_list)
        self.sperateMol()
        self.drawPng()


if __name__ == '__main__':
    # Xtb_Path(10,10).main()
    # G16_xtbTS(5).main()
    G16_TS(10,5).main()
    # G16_TS_Ref(5).main()
    # G16_Irc(2).main()

    # Orca_Ts(3).main()
    # Orca_Irc(3).main()
    pass
#========================================================#
#                  Job related functions                 #
#========================================================#

import os
import time

import pandas as pd

from ._setting_share import SharedSetting

def getJobList(mode, name):
    ### Need to modify according to the way to inquire the cluster job
    job_list = os.popen(
        "jobe | grep -w %s | grep %s" % (mode, name)).readlines()
    return job_list
    
def handleInfo(work_space, job_list):
    return [['%s/' % str(work_space) 
        + i.split()[-1]
            + '.pbs', i.split()[0]]
                for i in job_list]

class JobSubmit(SharedSetting):

    def __init__(self, csv, workspace):
        super().__init__()
        self.csv = csv
        self.workspace = workspace

    def submit(self):
        start_t = time.ctime()

        while True:
            df = pd.read_csv(self.csv)
            df.set_index(['Path'], inplace=True)

            if len(list(df[df['State'] == 'R'].index)) <= self.server_config['qsub_num']:
                unsub_list = list(df[df['State'] == 'U'].index)
                for job in unsub_list:
                    ### Need to modify according to the way to submit the cluster job
                    os.system('qsub %s >> job_list' % job)

                    df.loc[job, 'State'] = 'R'
                    if len(list(df[df['State'] == 'R'].index)) >= self.server_config['qsub_num']:
                        break
                df.to_csv(self.csv)
            
            run_list = handleInfo(self.work_space, getJobList('R', type))
            fin_list = handleInfo(self.work_space, getJobList('C', type))
            
            for job in run_list:
                df.loc[job[0], ('State', 'Id')] = ['R', job[1]]
            for job in fin_list:
                df.loc[job[0], 'State'] = 'C'
            
            df.to_csv(self.csv)

            fin_nums = len(df[df['State'] == 'C']['State'].tolist())
            tot_nums = len(df['State'].tolist())
            
            current_t = time.ctime()
            printout = '[%s/%s]  Start Time: %s; Current Time: %s' % \
                (fin_nums, tot_nums, start_t, current_t)
            print(printout)
            
            state_list = list(df.loc[:, 'State'])
            if state_list == ['C'] * len(state_list):
                break
            
            time.sleep(self.server_config['sleep_t'])

def submitJobs(type, csv, work_space, qsub_num=25, sleep_t=120):
    start_t = time.ctime()

    while True:
        df = pd.read_csv(csv)
        df.set_index(['Path'], inplace=True)

        if len(list(df[df['State'] == 'R'].index)) <= qsub_num:
            unsub_list = list(df[df['State'] == 'U'].index)
            for job in unsub_list:
                ### Need to modify according to the way to submit the cluster job
                os.system('qsub %s >> job_list' % job)

                df.loc[job, 'State'] = 'R'
                if len(list(df[df['State'] == 'R'].index)) >= qsub_num:
                    break
            df.to_csv(csv)
        
        run_list = handleInfo(work_space, getJobList('R', type))
        fin_list = handleInfo(work_space, getJobList('C', type))
        
        for job in run_list:
            df.loc[job[0], ('State', 'Id')] = ['R', job[1]]
        for job in fin_list:
            df.loc[job[0], 'State'] = 'C'
        
        df.to_csv(csv)

        fin_nums = len(df[df['State'] == 'C']['State'].tolist())
        tot_nums = len(df['State'].tolist())
        
        current_t = time.ctime()
        printout = '[%s/%s]  Start Time: %s; Current Time: %s' % \
            (fin_nums, tot_nums, start_t, current_t)
        print(printout)
        
        state_list = list(df.loc[:, 'State'])
        if state_list == ['C'] * len(state_list):
            break
        
        time.sleep(sleep_t)
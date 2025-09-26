#========================================================#
#                  CSV related functions                 #
#========================================================#

import re
import pandas as pd


def csv2List(csv):
    return pd.read_csv(csv).values.tolist()


def recordJobCSV(type, path_list):
    data = []
    for path in path_list:
        data.append([path, False])

    df = pd.DataFrame(data, columns=['Path', 'Submitted'], 
        index=[i + 1 for i in range(len(data))])
    df.to_csv('%s_job.csv' % (type))


def recordStateCSV(type, workspace):
    pbs_list = sorted(
        [str(pbs) for pbs in workspace.glob('*.pbs')],
            key = lambda l : int(re.findall('\d+', str(l))[0]))
    jobState = []
    for pbs in pbs_list:
        jobState.append([pbs, 'U', None])
    df = pd.DataFrame(jobState, columns=['Path', 'State', 'Id'])
    df.set_index(['Path'], inplace=True)
    df.to_csv('%s_state.csv' % (type))
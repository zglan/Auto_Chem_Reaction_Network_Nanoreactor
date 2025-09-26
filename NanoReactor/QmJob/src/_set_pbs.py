#========================================================#
#                  PBS related functions                 #
#========================================================#

from ._setting_default import default_setting
pbsContent = default_setting['pbs_content']
orcaLoc = default_setting['loc_config']['orca_loc']
nproc = default_setting['comp_config']["nproc"]
    

def writePbs(pbsPath, job_type, idx, pal, snum, enum, workspace, command):
    with open(pbsPath, 'w') as pbs:
        pbs.write(pbsContent % 
            (job_type, idx, pal, 
                workspace, snum, enum, command))

def buildJobPbs(pal, job_type, group_list, workspace):
    if job_type == 'ts':
        command = 'g16 *.gjf >> error 2>&1'
    elif job_type == 'ts_ref':
        command = 'g16 *.gjf >> error 2>&1'
    elif job_type ==  'irc':
        command = 'g16 *.gjf >> error 2>&1'
    elif job_type ==  'opt':
        command = 'g16 *.gjf >> error 2>&1'
    elif job_type == 'sp' or job_type == 're_sp':
        command = 'g16 *.gjf >> error 2>&1'

    elif job_type == 'orcaNeb':
        command = f'{orcaLoc} neb.inp > neb.out'
    elif job_type == 'orcaTs':
        command = f'{orcaLoc} ts.inp > ts.out'
    elif job_type == 'orcaIrc':
        command = f'{orcaLoc} irc.inp > irc.out'

    elif job_type == 'xtbPath':
        command = f'export OMP_NUM_THREADS={nproc}\nxtb start.xyz --path end.xyz --input path.inp > path.log'
    
    for idx, job in enumerate(group_list):
        pbsPath = workspace.joinpath('%s-%s.pbs' %(job_type, idx + 1))
        writePbs(pbsPath, job_type, idx + 1, 
                pal, job[0][-1], job[-1][-1], 
                workspace, command)
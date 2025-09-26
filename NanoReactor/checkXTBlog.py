#!/usr/bin/env python3
# check XTB log


from pathlib import Path
from ioToolkit import *


def checkError(lines):
    errorList = []
    time = '0d 0h 0min 0sec'
    for idx, line in enumerate(lines):
        if 'xtb_calculator_singlepoint: Electronic structure method terminated' in line:
            errorList.append('single point / electronic structure')
        if 'scf: Setup of Coulomb evaluator failed' in line:
            errorList.append('scf / coulomb evaluator')
        if 'scf: Self consistent charge iterator did not converge' in line:
            errorList.append('scf / consistent charge')
        if 'type_latticepoint_update: Could not generate lattice points' in line:
            errorList.append('latticepoint update / not generate')
        if 'MD is unstable, emergency exit' in line:
            errorList.append('md / unstable')
        if 'no convergence in svdcmp' in line:
            errorList.append('svdcmp / no convergence')
    
        if 'total:' in line:
            time = ' '.join(''.join(lines[idx + 1].split()[2: ]).split(','))
    return errorList, time


def checkAll():    
    fail_folder = Path('./fail')
    folder_list = getLocalFolderList()
    log_list = extendFolderList(folder_list, 'log')
    
    for log in log_list:
        folder = log.parents[0]
        
        lines = loadFile2List(log)
        error, time = checkError(lines)
        
        print(str(log) + ':  ' + time)
        if time == '0d 0h 0min 0sec' or len(error) > 0:
            mkdir(fail_folder)
            moveFolderList(folder, Path('./fail').joinpath(folder))


checkAll()
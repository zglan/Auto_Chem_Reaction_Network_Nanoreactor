#!/usr/bin/env python3

import argparse
from pathlib import Path
from collections import Counter

import pandas as pd
from tqdm import tqdm

from src._toolkit import *


def unite_forwardBack(d):
    for key_1 in list(d.keys()):
        for key_2 in list(d.keys()):
            if key_1 == '->'.join(key_2.split('->')[::-1]):
                key_new = key_1.replace('->', '<->')
                d[key_new] = d[key_1] + d[key_2]
                d.pop(key_1)
                d.pop(key_2)
    return d


def sort_by_times(d):
    return dict(
        sorted(d.items(), key=lambda x: len(x[1]), reverse=True))


def sort_by_atomNums(d):
    time_dict = {}
    for type in d.keys():
        atom = [x for x in 
             re.split("[^a-zA-Z]*", type) if x]
        n = int(len("".join(atom)) / 2)
        time = len(d[type])
        if time not in time_dict.keys():
            time_dict[time] = [(n, type, d[type])]
        elif time in time_dict.keys():
            time_dict[time].append((n, type, d[type]))
    
    new_d = {}
    for time in time_dict.keys():
        time_dict[time] = sorted(time_dict[time], key=lambda x: x[0])
        for i in time_dict[time]:
            reaction, info = i[1], i[-1]
            new_d[reaction] = info
    return new_d


def ts_preprocess(job_type, step, unite, sort):
    folder_list = getLocalFolderList()  
    if step == None:
        step = len(folder_list)
    else:
        if step >= len(folder_list):
            step = len(folder_list)
            print(f'Folder number does not fit, correct to {step}')
    
    mkdir('./%s' % (job_type))
    
    folderList_group = [
        folder_list[i: i + step] 
            for i in range(0, len(folder_list), step)]
    infoList_group = [
        extendFolderList(folder_list, '%sjob.info' % (job_type)) 
            for folder_list in folderList_group]
    
    for index, group in enumerate(infoList_group):
        group_path = Path(
            './%s/%s' % (job_type, infoList_group.index(group) + 1))
        
        print('*** Operating Grouping Job-%s ***' % (index + 1))
        mkdir(group_path)

        reac_dict = {}
        for info in group:
            lines = loadFile2List(info)
            for line in lines:
                idx, type, time, l = line.split('___')
                reac_dict.setdefault(type,[(info.parent, time, idx)])\
                    .append((info.parent, time, idx))
        
        if unite:
            reac_dict = unite_forwardBack(reac_dict)
        
        if sort:
            reac_dict = sort_by_times(reac_dict)
            reac_dict = sort_by_atomNums(reac_dict)

        for idx, type in enumerate(list(reac_dict.keys())):
            num = len(reac_dict[type])
            
            writeFile(
                Path('./%s/%s/info' % (job_type, infoList_group.index(group) + 1)),
                  'a', type + '   ' + str(num) + '\n')
            mkdir(group_path.joinpath(str(idx + 1)))

            for reacinfo in reac_dict[type]:
                srcFolder = reacinfo[0]\
                    .joinpath('%s' % (job_type))\
                    .joinpath(reacinfo[-1])
                dstFolder = group_path\
                    .joinpath(str(idx + 1))\
                    .joinpath(str(reacinfo[0]) + '_' + reacinfo[-1])
                copyFolder(srcFolder, dstFolder)
                writeFile(dstFolder.parent.joinpath('info'), 'w', type)
                print(f'''  {dstFolder}''')

#########################################

def update_dict(exist_dict, new_dict):
    for k, v in new_dict.items():
        if k not in exist_dict:
            exist_dict[k] = v
        else:
            exist_dict[k].extend(v)


def opt_path(smi_eq, smi_opt, xyz_opt, opt_n, opt_eq):
    if not smi_eq:
        opt_n.setdefault(smi_opt, []).append(xyz_opt)
    else:
        opt_eq.setdefault(smi_opt, []).append(xyz_opt)
    return opt_n, opt_eq


def get_ori_info(info, parent, job_idx):
    rows = info[info['Index'].isin([1, 6])]

    for idx, row in rows.iterrows():
        t_idx = row['Index']
        spin = row['Spin']
        chrg = row['Chrg']
        
        info = {'job_idx': 
                ['-'.join(map(str, (parent, job_idx)))], 
            'spin_chrg': [(spin, chrg)]}
        smi, smi_opt, smi_eq  = row['SMI'], row['SMI_Opt'], row['Smi_equal']
        xyz, xyz_opt = joinPaths(
            parent, str('xtbopt'), str(job_idx), row['Xyz']), \
                joinPaths(
                    parent, str('xtbopt'), str(job_idx), row['Opt_Xyz'])
        
        if t_idx == 1:
            r_smi = smi
            r_opt_n = {}
            r_opt_eq = {}
            ori_r = {}
            ori_r.setdefault(smi, []).append(xyz)
            r_opt_n, r_opt_eq = opt_path(
                smi_eq, smi_opt, xyz_opt, r_opt_n, r_opt_eq)

        elif t_idx == 6:
            p_smi = smi
            p_opt_n = {}
            p_opt_eq = {}
            ori_p = {}
            ori_p.setdefault(smi, []).append(xyz)
            p_opt_n, p_opt_eq = opt_path(
                smi_eq, smi_opt, xyz_opt, p_opt_n, p_opt_eq)
            
    return r_opt_n, p_opt_n,  \
            ori_r, ori_p, r_smi, p_smi, \
            r_opt_eq, p_opt_eq, info


def update_info(r_smi, p_smi, 
            info, parent, job_idx, 
            r_opt_n, p_opt_n,
            ori_r, ori_p,
            r_opt_eq, p_opt_eq):
    rows = info[info['Index'].isin([2, 3, 4, 5])]

    for idx, row in rows.iterrows():
        t_idx = row['Index']
        smi, smi_opt, smi_eq  = row['SMI'], row['SMI_Opt'], row['Smi_equal']
        xyz, xyz_opt = joinPaths(
            parent, str('xtbopt'), str(job_idx), row['Xyz']), \
                joinPaths(
                    parent, str('xtbopt'), str(job_idx), row['Opt_Xyz'])
        
        if 1 < t_idx <= 3:                  
            if r_smi != smi:
                print('  %s/%s/%s::reactant error'%(str(parent), job_idx, t_idx))
                # print(info)
            else:
                ori_r.setdefault(smi, []).append(xyz)
                r_opt_n, r_opt_eq = opt_path(
                    smi_eq, smi_opt, xyz_opt, r_opt_n, r_opt_eq)
                    
        elif 4 <= t_idx < 6:
            if p_smi != smi:
                print('  %s/%s/%s::product error'%(str(parent), job_idx, t_idx))
                # print(info)
            else:
                ori_p.setdefault(smi, []).append(xyz)
                p_opt_n, p_opt_eq = opt_path(
                    smi_eq, smi_opt, xyz_opt, p_opt_n, p_opt_eq)

    return r_opt_n, p_opt_n, ori_r, ori_p, r_opt_eq, p_opt_eq


def integrate_status(r_opt_n, p_opt_n, ori_r, ori_p,
        r_opt_eq, p_opt_eq, info):
    status = {}
    status['reactant_no_equal'] = r_opt_n
    status['product_no_equal'] = p_opt_n
    status['reactant_equal'] = r_opt_eq
    status['product_equal'] = p_opt_eq
    status['reactant'] = ori_r
    status['product'] = ori_p
    status['info'] = info
    return status


def update_react(react_dict, r_smi, p_smi, status):
    key = (r_smi, p_smi)
    if key in react_dict:
        exist_status = react_dict[key]
        update_dict(
            exist_status["reactant_no_equal"], status['reactant_no_equal'])
        update_dict(
            exist_status["product_no_equal"], status['product_no_equal'])

        update_dict(
            exist_status["reactant_equal"], status['reactant_equal'])
        update_dict(
            exist_status["product_equal"], status['product_equal'])

        update_dict(exist_status["reactant"], status['reactant'])
        update_dict(exist_status["product"], status['product'])

        update_dict(exist_status["info"], status['info'])
    else:
        react_dict[key] = status


def integrate_xtb_opt(job_type, step):
    folder_list = getLocalFolderList()  
    if step == None:
        step = len(folder_list)
    else:
        if step >= len(folder_list):
            step = len(folder_list)
            print(f'Folder number does not fit, correct to {step}')
    
    mkdir('./%s' % (job_type))
    folderList_group = [
        folder_list[i: i + step] 
            for i in range(0, len(folder_list), step)]
    infoList_group = [
        extendFolderList(folder_list, '_opt_xtb_data.csv') 
            for folder_list in folderList_group]

    for group in infoList_group:
        group_path = Path('./%s/%s' % (
            job_type, infoList_group.index(group) + 1))
        mkdir(group_path)

        react_dict = {}
        for info in group:
            parent = info.parent
            df = pd.read_csv(info)
            job_indexs = df['Job_index'].unique()
            for job_idx in job_indexs:
                info = df[df['Job_index']==job_idx]

                r_opt_n, p_opt_n, \
                    ori_r, ori_p, \
                    r_smi, p_smi, \
                    r_opt_eq, p_opt_eq, count  = get_ori_info(
                        info, parent, job_idx)
                
                r_opt_n, p_opt_n, \
                    ori_r, ori_p,\
                    r_opt_eq, p_opt_eq = update_info(
                        r_smi, p_smi, 
                        info, parent, job_idx, 
                        r_opt_n, p_opt_n,
                        ori_r, ori_p,
                        r_opt_eq, p_opt_eq)

                status = integrate_status(
                    r_opt_n, p_opt_n, 
                    ori_r, ori_p,
                    r_opt_eq, p_opt_eq, count)
                
                update_react(react_dict, r_smi, p_smi, status)

    return react_dict, group_path


def operate_xtb_opt(react_dict, group_path, freq):
    
    def copy_xyzs(xyz_dict, path, limit=3):
        for idx, key in enumerate(xyz_dict.keys()):
            dst = joinPaths(path, str(idx + 1))
            mkdir(dst)

            for xyz in xyz_dict[key][:limit]:
                copyFile(xyz, dst)

    def built_dir(path, name, xyz_dict):
        if len(xyz_dict.keys()) != 0:
            sec_path = joinPaths(path, name)
            mkdir(sec_path)
            copy_xyzs(xyz_dict, sec_path)

    for idx, key in enumerate(tqdm(list(react_dict.keys()))):
        stautus = react_dict[key]
        ori_r, ori_p = stautus['reactant'], stautus['product']
        r_opt_n, p_opt_n = stautus['reactant_no_equal'], \
            stautus['product_no_equal']
        r_opt_eq, p_opt_eq = stautus['reactant_equal'],  \
            stautus['product_equal']
        count = len(stautus['info']['job_idx'])

        if freq == None:
            pass
        else:
            if count < freq:
                continue

        path = joinPaths(group_path, str(idx + 1))
        mkdir(path)

        collects = [ori_r, ori_p, 
            r_opt_n, p_opt_n,     
            r_opt_eq, p_opt_eq]
        names = ['_Re_Ori', '_Pr_Ori', '_Re_Opt_No_Eq', 
                 '_Pr_Opt_No_Eq', '_Re_Opt_Eq', '_Pr_Opt_Eq', ]
        
        for (col, name) in zip(collects, names):
            built_dir(path, name, col)
        

def react_dict_2_dataframe(react_dict):
    concise_data = []
    path_data = []
    for (r_smi, p_smi), status in react_dict.items():
        concise_entry = {
            'r_smi': r_smi,
            'p_smi': p_smi
        }
        detailed_entry = {
            'r_smi': r_smi,
            'p_smi': p_smi,
            'spin_chrg': status['info']['spin_chrg'],
            'job_idx': status['info']['job_idx']
        }

        for key in ['reactant_no_equal', 'product_no_equal', 
                    'reactant_equal', 'product_equal',
                    'reactant', 'product', 'info']:
            if key == 'info':              
                count = len(status[key]['job_idx']) * 3
                concise_entry['count'] = count
                concise_entry['spin_chrg'] = Counter(
                    status[key]['spin_chrg']).most_common()[0][0]

            else:
                concise_entry[key] = [
                    (smi, idx + 1, len(paths)) for idx, (smi, paths) 
                        in enumerate(status[key].items())]

                paths_as_str = [
                    smi + ': ' + ', '.join([str(path) for path in paths]) 
                        for smi, paths in status[key].items()]
                detailed_entry[key] = ' | '.join(paths_as_str)
            
        concise_data.append(concise_entry)
        path_data.append(detailed_entry)

    df_detailed = pd.DataFrame(path_data)
    df_concise = pd.DataFrame(concise_data)

    return df_detailed, df_concise


def sort_react_dict_by_count(react_dict):
    react_items = list(react_dict.items())
    sorted_react_items = sorted(react_items, 
        key=lambda item: len(item[1]['info']['job_idx']), reverse=True)
    sorted_react_dict = {item[0]: item[1] for item in sorted_react_items}
    return sorted_react_dict


def intermediate_preprocess(step, freq):
    react_dict, group_path = integrate_xtb_opt('xtbopt', step)
    sort_react_dict = sort_react_dict_by_count(react_dict)
    df_detailed, df_concise = react_dict_2_dataframe(sort_react_dict)
    df_detailed.to_csv(
        joinPaths(group_path, 'react_detailed_data.csv'), index=False)
    df_concise.to_csv(
        joinPaths(group_path, 'react_concise_data.csv'), index=False)
    operate_xtb_opt(sort_react_dict, group_path, freq)

#########################################

def commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--step', type=int, default=None)
    parser.add_argument(
        '--freq', type=int, default=None)
    parser.add_argument(
        '--unite', action='store_true')
    parser.add_argument(
        '--sort', action='store_true')
    
    parser.add_argument(
        '--ts', action='store_true')
    parser.add_argument(
        '--neb', action='store_true')
    parser.add_argument(
        '--path', action='store_true')
    parser.add_argument(
        '--xtbopt', action='store_true')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    step = commandline().step
    sort = commandline().sort
    freq = commandline().freq
    unite = commandline().unite

    if commandline().ts:
        ts_preprocess('ts', step, unite, sort)
    if commandline().neb:
        ts_preprocess('neb', step, unite, sort)
    if commandline().path:
        ts_preprocess('path', step, unite, sort)
    if commandline().xtbopt:
        intermediate_preprocess(step, freq)
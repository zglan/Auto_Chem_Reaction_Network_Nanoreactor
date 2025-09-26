#!/usr/bin/env python3

import os
import re
import shutil

from pathlib import Path

#================================#
#          I/O Toolkit           #
#================================#

def mkdir(path):
    """Build a folder if it doesn't exist"""
    path = Path(path)
    if Path.exists(path) == False:
        Path.mkdir(path, parents=True)


def rmf(path):
    """Remove a file"""
    path = Path(path)
    if Path.exists(path) == True:
        path.unlink()


def rmdir(path):
    """Remove a folder"""
    path = Path(path)
    if Path.exists(path) == True:
        shutil.rmtree(path)


def copyFile(src, dst):
    """Copy a file"""
    try:
        shutil.copy(src, dst)
    except:
        print('IOError')
    

def copyFolder(src, dst):
    """Copy a folder"""
    try:
        shutil.copytree(src, dst)
    except:
        pass

def coverFolder(src, dst):
    """Cover a folder"""
    if os.path.exists(dst):
        shutil.rmtree(dst)
    try:
        shutil.copytree(src,dst)
    except:
        pass


def moveFolderList(src, dst):
    """Move a folder"""
    try:
        shutil.move(src, dst)
    except:
        pass


def loadFile2List(path):
    """Load a file content to list"""
    with open(path) as file:
        return file.readlines()


def writeFile(path, mode, string):
    """Write a string content to file"""
    with open(path, mode) as file:
        file.write(string)


def getLocalFolderList():
    """Get local folders with numbers"""
    local = Path('.')
    folder_list = sorted(
        [child for child in local.iterdir()
            if child.is_dir() and str(child).isdigit()],
        key = lambda l : int(re.findall('\d+', str(l))[0]))
    return folder_list


def getSubFolderList(path):
    """Get all sub folders"""
    folder_list = sorted(
        [child for child in path.iterdir()
            if child.is_dir()],
        key = lambda l : int(re.findall('\d+', str(l))[0]))
    return folder_list


def getSelectedFolderList(path):
    """Get selected folders with numbers"""
    path = Path(path)
    folder_list = sorted(
        [child for child in path.iterdir()
            if child.is_dir()],
        key = lambda l : int(re.findall('\d+', str(l))[0]))
    return folder_list


def extendFolderList(folder_list, name):
    """Extend a folder list to a file or a child folder list
        [file1, file2, ...] or [childfolder1, ...]
    """
    file_list = [f for f in 
        [folder.joinpath(name) for folder in folder_list]
            if f.exists()]
    return file_list


def getNameFileList(folder_list, name):
    """Get files with specified name in selected list of folder"""
    file_list = [
        [file for file in folder.glob('**/%s' % (name)) 
            if file.exists()]
                for folder in folder_list]
    return file_list


def getSuffixFileList(folder_list, suffix):
    """Get files with specified suffix in selected list of folder
        [[file1, file2, ...], [file3, file4], ...]
    """
    file_list = [
        [file for file in folder.glob('**/*.%s' % (suffix)) 
            if file.exists()]
                for folder in folder_list]
    return file_list


def getSuffixNumsFileList(folder_list, suffix, num):
    """Get files with specified suffix in selected list of folder
        the number of file will not take  
        more than the specified num in each folder 
    """
    file_list = []
    for folder in folder_list:
        tmp_list = sorted(
            [child for child in folder.iterdir() 
                if child.is_dir()]
                    )[:num]
        for file in tmp_list:
            file_list.extend(
                sorted(file.glob('*.%s' % (suffix))))
    return file_list


#========================================================#
#                  Network Analysis                      #
#========================================================#
import networkx as nx

def getConnection(edges_list):  
    G = nx.Graph()
    G.add_edges_from(edges_list) # list of edge tuples
    G_con = [
        [sorted(list(
            G.subgraph(c).nodes)), list(G.subgraph(c).edges)]
             for c in nx.connected_components(G)]
    return G_con

def diffNetworks(G1, G2):
    """
        Return: nx.Graph
    """
    G_new = nx.Graph()
    G_d1 = nx.difference(G1, G2)
    G_d2 = nx.difference(G2, G1)
    G_d = nx.compose(G_d1, G_d2)
    G_new.add_edges_from(G_d.edges)
    return G_new


#========================================================#
#                    Format Transform                    #
#========================================================#
import numpy as np

def atom_coord2array(atom_coord, idx_list=None):
    atom = []
    coord = []
    if idx_list != None:
        for idx, line in enumerate(atom_coord):
            if idx in idx_list:
                coord.append(line[-1])
                atom.append(line[0])
        atom = np.array(atom)
        coord = np.array(coord)
    elif idx_list == None:
        atom = np.array([line[0] for line in atom_coord])
        coord = np.array([line[-1] for line in atom_coord])
    return atom, coord

def atom_coord2str(atom_coord, idx_list=None):
    str = ''
    for idx, line in enumerate(atom_coord):
        if idx_list != None and idx in idx_list:
            str += '%s  %.14f  %.14f  %.14f\n' % (
                line[0], line[1][0], line[1][1], line[1][2])
        elif idx_list == None:
            str += '%s  %.14f  %.14f  %.14f\n' % (
                line[0], line[1][0], line[1][1], line[1][2])
    return str

def array2atom_coord(atom_array, coord_array):
    atom_coord = []
    for idx, atom in enumerate(atom_array):
        atom_coord.append((atom, coord_array[idx].tolist()))
    return atom_coord

#===============================================#
#                     Else                      #
#===============================================#

def limitList(arr, up, down):
    arr = arr[arr <= up]
    arr = arr[arr >= down]
    return arr

def coordStr2fxyz(xyz, coord_str, atom_nums):
    with open('%s' % xyz, 'w') as f:
        f.write('%s\n\n%s' % (atom_nums, coord_str))

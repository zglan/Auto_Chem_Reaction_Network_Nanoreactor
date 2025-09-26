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
        the number of file will not take  more than the specified num in each folder 
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

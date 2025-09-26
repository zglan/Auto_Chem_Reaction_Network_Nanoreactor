#!/usr/bin/env python3

import os
import argparse
import operator

import numpy as np 

from link_judge import linkJudge
from ioToolkit import rmf
from logger import MyLog



myLogger = MyLog().logger

def commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--inputfile', 
        help='Input gaussion log file', required=True
    )
    parser.add_argument(
        '-t', '--type',
        help='Set gaussion job type', required=True
    )
    args = parser.parse_args()
    return args


def listToStr(l):
    s = ''
    for i in l:
        s += '%s  %8f  %8f  %8f' % (i[0], i[1], i[2], i[3]) + '\n'
    return s

def writeXyz(name, coord):
    num = str(len(coord)) + '\n\n'
    coord = listToStr(coord)
    with open('%s.xyz' % (name), 'w') as xyz:
        xyz.write(num + coord)

def xyz2coord(f):
    with open(f) as xyz:
        coord = ''.join(xyz.readlines()[2:])
    return coord

def readOrcaInp(f):
    with open(f) as inp:
        line = inp.readlines()[-1].split('')
        print(line)



class SetGaussionJob:
    
    def __init__(self, **kwargs):
        
        necessary = ['jobname','jobtype', 'coord']
        default = {
            'state': 'sin',
            'reactiontype': 'radical',
            'nproc': 4, 'mem': 4000,
            'base': 'SVP', 'qmmethod': 'B3LYP',
            'maxcyc': 300, 'maxstep': 10, 
            'maxpoints': 300, 'stepsize': 10}
            # when maxcyc is set large than 100, it will be set 100
        addition = ['elenums', 'charge', 'multiplicity']

        # addition = ['cartesian', 'em=GD3BJ'] # IRC: FormBX had a problem

        for arg in default:
            kwargs.setdefault(arg, default[arg])
        
        for arg in necessary:
            kwargs[arg] = kwargs[arg]
        
        for arg in addition:
            kwargs.setdefault(arg, None)

        self.__dict__.update(kwargs)

    def setChrgSpin(self):
        if self.charge != None and self.multiplicity != None:
            chrg = self.charge
            spin = self.multiplicity

        elif self.state == 'sin' and self.reactiontype == 'ion':
            spin = 1
            if self.charge == None:
                chrg = None 
            elif 0.3 <= self.charge <= 0.65:
                chrg = None
            elif -0.3 <= self.charge <= 0.3:
                chrg = 0
            elif 0.65 <= self.charge <= 1.4:
                chrg = 1
            elif 1.4 <= self.charge <= 2.4:
                chrg = 2
            elif 2.4 <= self.charge <= 3.4:
                chrg = 3
            elif self.charge >= 3.4:
                chrg = None
            elif -0.65 <= self.charge <= -0.3:
                chrg = None
            elif -1.4 <= self.charge <= -0.65:
                chrg = -1
            elif -2.4 <= self.charge <= -1.4:
                chrg = -2
            elif -3.4 <= self.charge <= -2.4:
                chrg = -3
            elif self.charge <= -3.4:
                chrg = None
        
        elif self.state == 'sin' and self.reactiontype == 'radical':
            if self.state == 'sin' and self.elenums % 2 == 0:
                chrg, spin = 0, 1
            else:
                chrg, spin = 0, 2

        elif self.state == 'tri':
            chrg, spin = 0, 3

        return str(chrg), str(spin)


    def setTittle(self):
        if self.jobtype == 'ts':
            tittle = '#opt=(%s,calcfc,recalc=5,noeigen,maxcyc=%s) freq %s/%s em=GD3BJ\n\n' % (
                self.jobtype, self.maxcyc, self.qmmethod, self.base)

        if self.jobtype == 'opt':
            tittle = '#opt=(tight,maxcyc=%s,maxstep=%s) %s/%s em=GD3BJ\n\n' % (
                self.maxcyc, self.maxstep, self.qmmethod, self.base)

        if self.jobtype == 'irc':
            tittle = '#irc=(LQA,calcfc,noeigen,maxpoints=%s,stepsize=%s) %s/%s em=GD3BJ\n\n' % (
                self.maxcyc, self.stepsize, self.qmmethod, self.base)
        
        if self.jobtype == 'singlepoint':
            tittle = '#%s/%s em=GD3BJ\n\n' % (
                self.qmmethod, self.base)
            
        return tittle

    def setJob(self):
        with open('%s.gjf' % (self.jobname), 'w') as gjf:
            gjf.write(
                '%%nproc=%s\n%%mem=%sMB\n'
                % (self.nproc, self.mem))
                #'%%chk=%s\n%%nproc=%s\n%%mem=%sMB\n'
                #% (self.jobname.split('/')[-1], self.nproc, self.mem))

            gjf.write(
                self.setTittle())
            
            gjf.write(
                '%s-%s\n\n' 
                % (self.jobtype, self.jobname.split('/')[-1]))

            gjf.write(
                '%s\n%s\n' % (
                ' '.join(self.setChrgSpin()), self.coord))


class ReadGaussionJob:
    
    def __init__(self, jobname):
        self.jobname = jobname
        self.log = self.getLog()
        self.error = self.checkError()
        if self.error == None:
            self.chrg, self.spin, self.atomnums, self.atomlist = self.getAtomInfo()
            self.coord = self.extractCoord()
        else:
            myLogger.error('[ %s ] Exists  %s  ' % (self.jobname, self.error))
            # exit()
    
    def getLog(self):
        with open(self.jobname) as log:
            return log.readlines()
    
    def checkError(self):
        error = None
        if len(self.log) <= 3:
            error = 'G16WorkNoStarts'
        for line in self.log:
            if 'Wanted an integer as input.' in line:
                error = 'ErrorChrgSpinSet'
            if 'The combination of multiplicity' and 'electrons is impossible.' in line:
                error = 'ErrorChrgSpinSet'
            if 'Convergence criterion not met.' in line:
                error = 'ErrorConvergence'
        return error
    
    def getTime(self):
        for index, line in enumerate(self.log):
            if 'Job cpu time' in line:
                ct = line.split()
                ct = np.array([int(ct[3]), int(ct[5]), int(ct[7]), float(ct[9])])
                print(ct)
            if 'Elapsed time' in line:
                wt = line.split()
                wt = np.array([int(wt[2]), int(wt[4]), int(wt[6]), float(wt[8])])
                print(wt)
        return ct, wt

    def getAtomInfo(self):
        atomnums = 0
        atomlist = []
        for index, line in enumerate(self.log):
            if 'Charge' and 'Multiplicity' in line and atomnums == 0:
                chrg = int(line.split()[2])
                spin = int(line.split()[-1])
                judge = True
                while judge:
                    if self.log[index + 1 + atomnums] == ' \n':
                        judge = False
                    else:
                        atomlist.append(self.log[index + 1 + atomnums].split()[0])
                    atomnums += 1
        return chrg, spin, atomnums, atomlist

    def extractCoord(self):
        coords = []
        for index, line in enumerate(self.log):
            if len(self.atomlist) == 0:
                break
            if 'Input orientation:' in line:
                coord = [
                    [self.atomlist[idx]] + [float(x) for x in acoord.strip().split()[-3: ]]
                        for idx, acoord in enumerate(
                            self.log[index + 5: index + 4 + self.atomnums])
                    ]
                coords.append(coord)
        return coords

    def extractStandardOrientation(self):
        self.chrg, self.spin, self.atomnums, self.atomlist = self.getAtomInfo()
        standardOri = []
        for index, line in enumerate(self.log):
            if len(self.atomlist) == 0:
                break
            if 'Standard orientation:' in line:
                coord = [
                    [self.atomlist[idx]] + [float(x) for x in acoord.strip().split()[-3: ]]
                        for idx, acoord in enumerate(
                            self.log[index + 5: index + 4 + self.atomnums])
                    ]
                standardOri.append(coord)
        return standardOri


class ReadOptJob(ReadGaussionJob): 
    
    def extractEnergy(self):
        energy = []
        for index, line in enumerate(self.log):
            if 'SCF Done' in line:
                energy.extend([
                    float(line.split()[4]) * 627.510])
        try:
            energy = energy[-1]
        except:
            pass
        return energy
    
    def extractTerminType(self):
        type = None
        for line in self.log:
            if 'Optimization completed.' in line:
                type = 'Normal'
        return type
    
    def extractAtomChrg(self):
        for index, line in enumerate(self.log):
            if 'Mulliken charges:' in line:
                atomchrg = [float(i.split()[-1]) 
                    for i in self.log[index + 2: index + 1 + self.atomnums]]
        return atomchrg
 
    def main(self):
        # if self.error == None and self.extractTerminType() == 'Normal':
            energy = self.extractEnergy()
            coords = self.extractStandardOrientation()[-1]
            # atomchrg  = self.extractAtomChrg()
            con = linkJudge(coords)
            
            for idx, i in enumerate(con):
                n_coord = []
                for j in i:
                    n_coord.append(coords[j])
                jobname = str(self.jobname).split('.')[0] + '-' + str(idx)
                SetGaussionJob(
                    jobname=jobname,
                    jobtype='singlepoint',
                    coord=listToStr(n_coord),
                    charge=99,
                    multiplicity=1).setJob()
                os.system('g16 ' + jobname)
                standard = ReadGaussionJob(jobname + '.log').extractStandardOrientation()[0]
                writeXyz(jobname, standard)
                rmf(jobname + '.log')
                rmf(jobname + '.gjf')


class ReadTsJob(ReadGaussionJob): 
    
    def extractEnergy(self):
        energy = []
        for index, line in enumerate(self.log):
            if 'SCF Done' in line:
                energy.extend([
                    float(line.split()[4]) * 627.510])
        try:
            energy = energy[-1]
        except:
            pass
        return energy

    def extractFreq(self):
        freq = []
        freqError = True
        for index, line in enumerate(self.log):
            if '1 imaginary frequencies (negative Signs)' in line and line.strip().split()[0] != 1:
                freqError = False
            if 'Frequencies --' in line:
                freq.extend(
                    [float(f) for f in line.split()[-3:]])
        return freq, freqError
    
    def extractTerminType(self):
        type = 'Normal'
        for line in self.log:
            if 'Normal termination' in line:
                pass
            
            if 'Error termination' in line:
                if 'FormBX had a problem' in line:
                    pass
        return type
            
    def setIrcJob(self):
        SetGaussionJob(
            nproc=8,
            jobname='%s-irc' % (self.jobname.split('.')[0]),
            jobtype='irc',
            coord=listToStr(
                self.extractCoord()[-1]),
            charge=self.chrg,
            multiplicity=self.spin
        ).setJob()
            
    def main(self):
        if self.error == None:
            freq, freqError = self.extractFreq()
            energy = self.extractEnergy()
            if freqError == False:
                self.setIrcJob()


class ReadIRCJob(ReadGaussionJob):
    
    def extractPointCoord(self):
        pointCoord = []
        for index, line in enumerate(self.log):
            if 'Input orientation' in line:
                # [(atom1, [x1, y1, z1]),
                #  (atom2, [x2, y2, z2]) ]
                coord = [
                    (self.atomlist[int(line.split()[0]) - 1], 
                        list(map(float, line.split()[-3:])))
                    for line in 
                        self.log[index + 5: index + 4 + self.atomnums]]
                    
                judge = True
                n = 0
                while judge:
                    if 'Point Number:' and 'Path Number:' in self.log[index + n]:
                        point = int(self.log[index + n].strip().split()[2])
                        path = int(self.log[index + n].strip().split()[-1])
                        judge = False
                        pointCoord.append((path, point, coord))
                    if 'Job cpu time' in self.log[index + n]:
                        judge = False
                    n += 1
        return pointCoord

    def extractPointEnergy(self):
        pointEnergy = []
        for index, line in enumerate(self.log):
            # extract energy
            if 'SCF Done:' in line:
                energy = float(line.strip().split()[4]) * 627.510 # hatree to kcal/mol
                judge = True
                n = 0
                # extract point and path
                while judge:                
                    if 'Point Number:' and 'Path Number:' in self.log[index + n]:
                        point = int(self.log[index + n].strip().split()[2])
                        path = int(self.log[index + n].strip().split()[-1])
                        judge = False
                    
                    n += 1
                pointEnergy.append((path, point, energy))
        return pointEnergy
        
    def extractTerminType(self):
        type = 'Normal'
        for line in self.log:
            if 'Reaction path calculation complete.' in line:
                pass
            if 'Problem with the distance matrix.' in line:
                type = 'Error'
        return type
    
    def getPath(self):
        if self.error == None and self.extractTerminType() == 'Normal':
            path = operator.itemgetter(0)
            pointEnergy = self.extractPointEnergy()
            pointCoord = self.extractPointCoord()

            first, second = [], []
            for onePoint in pointEnergy:
                if path(onePoint) == 1:
                    first.append(onePoint)
                if path(onePoint) == 2:
                    second.append(onePoint)
            energyPath = (first[0], first[-1], second[-1])
            
            first, second = [], []
            for onePoint in pointCoord:
                if path(onePoint) == 1:
                    first.append(onePoint)
                if path(onePoint) == 2:
                    second.append(onePoint)
            coordPath = (first[0], first[-1], second[-1])
            energyGap = abs(energyPath[1][-1] - energyPath[-1][-1])
        else:
            energyGap, energyPath, coordPath = None, None, None
        return energyGap, energyPath, coordPath
    
    def main(self):
        energyGap, energyPath, coordPath = self.getPath()
        print(energyPath[0][-1] - energyPath[1][-1], energyPath[0][-1] - energyPath[-1][-1])
            # print(coordPath[1][-1])
            # reacCoords = listToStr(coordPath[1][-1])
            # prodCoords = listToStr(coordPath[-1][-1])
            # SetGaussionJob(
            #     jobname=self.jobname.split('.')[0] + '-r-opt', 
            #     jobtype='opt',
            #     state='sin',
            #     charge=self.chrg,
            #     multiplicity=self.spin,
            #     coord=reacCoords).setJob()
            # SetGaussionJob(
            #     jobname=self.jobname.split('.')[0] + '-p-opt', 
            #     jobtype='opt',
            #     state='sin',
            #     charge=self.chrg,
            #     multiplicity=self.spin,
            #     coord=prodCoords).setJob()


def main():
    log = commandline().inputfile
    if commandline().type == 'ts':
        ReadTsJob(log).main()
    if commandline().type == 'irc':
        ReadIRCJob(log).main()
    if commandline().type == 'opt':
        ReadOptJob(log).main()
    if commandline().type == 'time':
        ReadGaussionJob(log).getTime()



if __name__ == '__main__':
    main()
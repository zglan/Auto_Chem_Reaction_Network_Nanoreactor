#!/usr/bin/env python3

import argparse
import operator

from _logger import MyLog
from _setting_share import SharedSetting


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


class G16_Set(SharedSetting):

    def __init__(self, gjf, job_type, coord, chrg=0, spin=1):
        super().__init__()
        self.gjf = gjf
        self.job_type = job_type
        self.coord = coord
        self.chrg = chrg
        self.spin = spin

    def setGjf(self):
        qm = f'''{self.qm_config["method"]}/{self.qm_config["g16_base"]} {self.disp_config["g16"]}'''

        if self.job_type == 'opt':
            addit = self.opt_config["g16_addit2"]
            outline = f'''#opt=({addit}maxcyc={self.opt_config["maxcyc"]},maxstep={self.opt_config["maxstep"]}) {qm}'''

        elif self.job_type == 'ts':
            addit = self.ts_config["g16_addit"]
            outline = f'''#opt=({addit}{self.job_type},recalc={self.ts_config['recalc']},maxcyc={self.ts_config['maxcyc']}) freq {qm}'''

        elif self.job_type == 'irc':
            addit = self.irc["g16_addit"]
            outline = f'''#irc=({addit}maxpoints={self.irc_config["maxpoints"]},stepsize={self.irc_config["stepsize"]}) {qm}'''
        
        comp = f'''%nproc={self.comp_config["nproc"]}\n%mem={self.comp_config["core"]}MB'''
        sub_line = f'''{self.job_type}-{self.gjf.split('/')[-1]}'''
        chrgspin = f'''{self.chrg} {self.spin}'''
        
        gjf = f'''{comp}\n{outline}\n\n{sub_line}\n\n{chrgspin}\n{self.coord}\n'''
        with open('%s' % (self.gjf), 'w') as f:
            f.write(gjf)

class G16_Read:
    
    def __init__(self, jobname):
        self.jobname = jobname
        self.log = self.getLog()
        self.error = self.checkError()
        if self.error == None:
            self.chrg, self.spin,\
                self.atomnums, self.atomlist = self.getAtomInfo()
            self.coord = self.extractCoord()
        else:
            myLogger.error('[ %s ] Exists  %s  ' % (self.jobname, self.error))
    
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
            if 'The combination of multiplicity' and \
            'electrons is impossible.' in line:
                error = 'ErrorChrgSpinSet'
            if 'Convergence criterion not met.' in line:
                error = 'ErrorConvergence'
        return error
    
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
                coord = [(self.atomlist[idx], 
                    [float(x) for x in acoord.strip().split()[-3: ]])
                        for idx, acoord in enumerate(
                            self.log[index + 5: index + 4 + self.atomnums])]
                coords.append(coord)
        return coords

    def extractStandardOrient(self):
        self.chrg, self.spin, self.atomnums, self.atomlist = self.getAtomInfo()
        standardOri = []
        for index, line in enumerate(self.log):
            if len(self.atomlist) == 0:
                break
            if 'Standard orientation:' in line:
                coord = [(self.atomlist[idx], 
                    [float(x) for x in acoord.strip().split()[-3: ]])
                        for idx, acoord in enumerate(
                            self.log[index + 5: index + 4 + self.atomnums])
                    ]
                standardOri.append(coord)
        return standardOri


class G16_Ropt(G16_Read): 
    
    def extractEnergy(self):
        energy = []
        for line in self.log:
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
 


class G16_Rts(G16_Read): 
    
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
        for line in self.log:
            if '1 imaginary frequencies (negative Signs)'\
                  in line and line.strip().split()[0] != 1:
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
        from _toolkit import atom_coord2str
        G16_Set(
            gjf='%s-irc' % (self.jobname.split('.')[0]),
            job_type='irc',
            coord=atom_coord2str(self.extractCoord()[-1]),
            chrg=self.chrg,
            spin=self.spin).setGjf()
            
    def ts2irc(self):
        if self.error == None:
            freq, freqError = self.extractFreq()
            energy = self.extractEnergy()
            if freqError == False:
                self.setIrcJob()


class G16_Rirc(G16_Read):
    
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
    
    def setOptJob(self, coords):
        G16_Set(
            gjf='%s-opt' % (self.jobname.split('.')[0]),
            job_type='opt',
            coord=coords,
            chrg=self.chrg,
            spin=self.spin).setGjf()

    def irc2opt(self):
        energyGap, energyPath, coordPath = self.getPath()
        if energyGap > 2:
            print('%.1f kcal/mol | %.1f kcal/mol' 
                  % (energyPath[0][-1] - energyPath[1][-1], energyPath[0][-1] - energyPath[-1][-1]))
            from _toolkit import atom_coord2str
            reacCoords = atom_coord2str(coordPath[1][-1])
            prodCoords = atom_coord2str(coordPath[-1][-1])
            self.setOptJob(reacCoords)
            self.setOptJob(prodCoords)


def main():
    log = commandline().inputfile
    if commandline().type == 'ts':
        G16_Rts(log).ts2irc()
    if commandline().type == 'irc':
        G16_Rirc(log).irc2opt()
    if commandline().type == 'opt':
        pass



if __name__ == '__main__':
    main()
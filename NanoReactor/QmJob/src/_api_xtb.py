import os

from pathlib import Path

spinp = '''
$chrg %s
$spin %s
$scc 
   maxiterations=800
   temp=6000.0000000000000
   broydamp=.4000000000000000
$write
   wiberg=true
   charges=true
$end'''

pathinp= '''
$chrg %s
$spin %s
$scc 
   maxiterations=800
   temp=6000.0000000000000
   broydamp=.4000000000000000
$path
   nrun=1
   npoint=50
   anopt=10
   kpush=0.003
   kpull=-0.015
   ppull=0.05
   alp=0.6
$end
'''

class Xtb_Set:
    
    def __init__(self, xyz, jobType='sp', chrg=0, spin=0):
        self.jobType = jobType
        self.chrg = chrg
        self.spin = spin

        self.xyz = xyz
        if jobType == 'path':
            self.inp = str(Path(self.xyz).parent) + '/path.inp'
            self.log = str(Path(self.xyz).parent) + '/path.log'
        else:
            self.inp = xyz.split('.')[0] + '.inp'
            self.log = xyz.split('.')[0] + '.log'


        self.atom_nums = self.getNums()
    
    def getNums(self):
        with open(self.xyz) as f:
            atom_nums = len(f.readlines()) - 2
            return atom_nums
    
    def makeInp(self):
        if self.jobType == 'sp':
            content = spinp % (self.chrg, self.spin)
        if self.jobType == 'path':
            content = pathinp % (self.chrg, self.spin)
        with open('%s' % (self.inp), 'w') as inp:
            inp.write(content)

    def runXTB(self):
        self.makeInp()
        if not os.path.exists(self.inp):
            raise  OSError
        
        if self.jobType == 'sp':
            os.system('xtb %s --input %s > %s' %  (self.xyz, self.inp, self.log))
        elif self.jobType == 'path':
            os.system('xtb start.xyz --path end.xyz --input %s  > %s' % (self.xyz, self.inp, self.log))

    def getChrg_charges(self):
        f_chrg = str(Path(self.xyz).parent) + '/charges'
        with open(f_chrg) as f:
            chrg_list = [float(line.split()[0]) for line in f.readlines()]
        nf_chrg = self.xyz.split('.')[0] + '.charges'
        os.rename(f_chrg, nf_chrg)
        return chrg_list

    def getChrg_log(self):
        with open(self.log) as f:
            lines = f.readlines()
            for idx, line in enumerate(lines):
                if  '#   Z          covCN         q      C6AA      α(0)' in line:
                    chrg_list  = [float(i.split()[4]) for i in lines[idx + 1: idx + 1 + self.atom_nums]]
                    return chrg_list

    def getFragmentChrg(self, idx_list):
        self.runXTB()
        chrg_list = self.getChrg_log()
        chrg = 0
        for idx in idx_list:
            chrg += chrg_list[idx]
        return chrg
    
    def getBondOrder_wbo(self):
        f_wbo = str(Path(self.xyz).parent) + '/wbo'
        with open(f_wbo) as f:
            wbo_list = [[
                (int(line.split()[0]) - 1, int(line.split()[1]) - 1), 
                    float(line.split()[-1])] 
                for line in f.readlines()]
        
        nf_wbo = self.xyz.split('.')[0] + '.wbo'
        os.rename(f_wbo, nf_wbo)
        return wbo_list
   
    def getFragmentInfo(self, idx_list):
        self.runXTB()
        fwbo_list = self.getBondOrder_wbo()
        fchrg_list = self.getChrg_charges()
       
        fgt_wbo_list, fgt_chrg_list = [], []
        fgt_chrg = 0
        for idx in idx_list:
            fgt_chrg += fchrg_list[idx]
            fgt_chrg_list.append((fchrg_list[idx], idx))
            for wbo in fwbo_list:
                if idx in wbo[0] and wbo not in fgt_wbo_list:
                    fgt_wbo_list.append(wbo)
        fgt_wbo_list.sort()
        return fgt_wbo_list, fgt_chrg_list, fgt_chrg

    def getFODPop_log(self):
        with open(self.log) as f:
            lines = f.readlines()
            for idx, line in enumerate(lines):
                if  'Loewdin FODpop      n(s)   n(p)   n(d)' in line:
                    chrg_list  = [float(i.split()[1]) for i in lines[idx + 1: idx + 1 + self.atom_nums]]
                    return chrg_list
    
    def getFragmentFODPop(self, idx_list):
        self.runXTB()
        pop_list = self.getPop_log()
        pop = 0
        for idx in idx_list:
            pop += pop_list[idx]
        return pop
    
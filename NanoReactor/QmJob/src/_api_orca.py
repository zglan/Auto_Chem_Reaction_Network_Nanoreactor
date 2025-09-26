from pathlib import Path

from ._logger import MyLog
from ._setting_share import SharedSetting

myLogger = MyLog().logger

class Orca_Set(SharedSetting):

    def __init__(self, inp, job_type, xyz, chrg=0, spin=1):
        super().__init__()
        self.inp = inp
        self.job_type = job_type
        self.xyz = xyz
        self.chrg = chrg
        self.spin = spin

    def setInp(self):      
        if self.job_type == 'neb':
            sxyz = self.xyz[0]
            outline = 'CI-NEB freq'
            details = f'''%%neb\n    NEB_END_XYZFILE "{sxyz}"\nend'''
            self.xyz = self.xyz[1]

        elif self.job_type == 'ts':
            outline = 'RIJCOSX optTS freq noautostart nopop miniprint'
            details = f'''%%geom\n    Calc_Hess true\n    Recalc_Hess {self.ts_config['recalc']}\n    MaxIter {self.ts_config['maxcyc']}\nend'''

        elif self.job_type == 'irc':
            outline = 'RIJCOSX freq noautostart nopop miniprint'
            details = f'''%%irc\n    MaxIter {self.irc_config['maxpoints']}\nend'''
        
        comp = f'''%%PAL NPROCS {self.comp_config["nproc"]} END \n%%maxcore {self.comp_config["core"]}'''
        xyz_file = f'''* XYZfile {self.xyz} {self.chrg} {self.spin}'''
        inp = f'''! {self.qm_config["method"]} {self.disp_config["orca"]} {self.qm_config["orca_base"]} {outline}\n{comp}\n{details}\n{xyz_file}'''

        with open('%s' % (self.inp), 'w') as f:
            f.write(inp)


class Orca_Read:

    def __init__(self, folder, name, job_type):
        self.name = name
        self.job_type = job_type
        self.folder = Path(folder)
        self.out = self.folder.joinpath('%s.out' % (name))
        self.tmp_list = [
            file for file in self.folder.glob('**/*.tmp') 
                if file.exists()]
       
    def get_eleE_Freq(self):
        freqS_list, freqE_list = [], []
        eleE_list = []
        
        with open(self.out) as out:
            lines = out.readlines()
            for idx, line in enumerate(lines):
                if 'VIBRATIONAL FREQUENCIES' in line:
                    freqS_list.append(idx + 5)
                elif 'NORMAL MODES' in line:
                    freqE_list.append(idx - 4)
                elif 'Electronic energy' in line:
                    eleE_list.append(float(line.split()[-2]) * 627.510)

        fin_img_freq = []
        fin_freq = lines[freqS_list[-1]: freqE_list[-1]]
        for l in fin_freq:
            if 'imaginary mode' in l:
                fin_img_freq.append(float(l.split()[1]))     

        fin_eleE = eleE_list[0]
        
        return fin_eleE, fin_img_freq
    
    def get_ePath(self):
        idx_list, ePath = [], []
        
        with open(self.out) as out:
            lines = out.readlines()
            for idx, line in enumerate(lines):
                if 'IRC PATH SUMMARY' in line:
                    idx_list.append(idx + 5)
                elif '<= TS' in line:
                    idx_list.append(idx)
                elif 'Timings for individual modules' in line:
                    idx_list.append(idx - 2)
        
        for idx in idx_list:
            ePath.append(float(lines[idx].split()[2]))
        
        return ePath
        
    def readTsOut(self):
        eleE, freq = 0, 0
        error = False
        if len(self.tmp_list) > 0:
            myLogger.error(
                '  Orca %s Job [ %s ] Failed  '
                  % (self.job_type, self.out))
        else:
            eleE, freq = self.get_eleE_Freq()
            if len(freq) > 1:
                myLogger.error(
                    '  Orca %s Job [ %s ] exists more than 1 imaginary mode'
                      % (self.job_type, self.out))
            else:
                freq = freq[0]
                error = True

        return eleE, freq, error

    def readIrcOut(self):
        ePath = []
        error = False
        if len(self.tmp_list) > 0:
            myLogger.error(
                '  Orca %s Job [ %s ] Failed  '
                  % (self.job_type, self.out))
        else:
            ePath = self.get_ePath()
            error = True
        return ePath, error

import json, os
import numpy as np
from src.tool.constant import kb_au, au2per_ps, onefs, amu2au
from src.nanoreactor.wall_potential import Nanoreactor
from src.mtd.metadynamics import MTD
from src.its.its_main import ITS
from typing import Dict, Any
from ase.io import read


class SimulationSettings:
    def __init__(self, config_path: str = "./job.json"):
        self.config = self._load_config(config_path)
        
        self.md = self.config.get("md", {})
        self.dt = self.md.get("dt", 1.0) * onefs                            # fs in au
        self.time = self.md.get("time", 5) * onefs * 1e3                    # ps in au
        self.T_target = self.md.get("Temperature")

        self.thermostat = self.md.get('thermostat', {})
        self.calculator_type = self.md.get('calculator_type', "GFN2-xTB")
        if self.thermostat:             
            print(self.thermostat)
            if self.thermostat.get('type') == 1:
                print("Langevin")
                self.gamma = self.thermostat.get('gamma', 7) / au2per_ps                  # ps-1
                self.c1 = np.exp(-self.gamma * self.dt / 2)
                self.c2 = np.sqrt((1 - self.c1 ** 2) * kb_au * self.T_target)
            elif self.thermostat.get('type') == 2:
                print("Berendsen") 
                self.tau = self.thermostat.get('tau', 200) * self.md.get("dt")

        self.mtd = self.config.get("mtd", {})
        if self.mtd:
            self.cvs_type = self.mtd.get('cvs_type', "RMSD")
            self.atom_ids = self.mtd.get('atom_ids', None)
            self.k_mtd = self.mtd.get('k_mtd', 0.2)
            self.sigma_mtd = self.mtd.get('sigma_mtd', 0.7)
            self.dump = self.mtd.get('dump', 5000)                  # Dump a reference structure after every 5000 steps.
            self.domtd = MTD(cvs_type=self.cvs_type,
                             atom_ids=self.atom_ids,
                             kpush=self.k_mtd,
                             sigma=self.sigma_mtd,
                             dump=self.dump)

        self.nano = self.config.get("nano", {})
        if self.nano:
            self.k_wall = self.nano.get('k_wall', 0.019)
            self.beta_wall = self.nano.get('beta_wall', 10.)
            self.R_wall = self.nano.get('R_wall')
            self.donano = Nanoreactor(k_wall=self.k_wall,
                                      beta_wall=self.beta_wall,
                                      R=self.R_wall)
        
        self.its = self.config.get("its", {})
        if self.its:
            self.bin_temp = self.its.get('bin_temp', 10)
            self.temp_low = self.its.get('temp_low', 300)
            self.temp_high = self.its.get('temp_high', 650)
            if self.temp_low > self.temp_high:
                raise ValueError("The lowest temperature Can not larger than highest temperature.")
            self.doits = ITS(T_target=self.T_target,
                             bin_temp=self.bin_temp,
                             temp_low=self.temp_low,
                             temp_high=self.temp_high)
        

        

    def _load_config(self, path: str) -> Dict[str, Any]:
        try:
            with open(path, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            raise RuntimeError(f"config file {path} is not exist!!")

    def get_module_config(self, module_name: str) -> Dict:
        return getattr(self, module_name, {})


class System:
    def __init__(self, xyz_file: str = "./stru.xyz", vel_file: str = "./vel.ini"):

        self.T_target = SimulationSettings().T_target
        if os.path.exists(xyz_file):
            self.Atom = read(f"{xyz_file}")
            xyz = self.Atom.get_positions()
            self.N = len(xyz)
            self.xyz = xyz / 5.2917720859e-1  # xyz: a.u.
            self.mass = np.reshape(self.Atom.get_masses() * amu2au, (self.N, 1))  # mass: a.u.
        else:
            raise FileExistsError(f"Flie {xyz_file} is not exist!!! Please Check the file path of {xyz_file}.")

        if os.path.exists(vel_file):
            self.v = np.loadtxt(f"{vel_file}")
        else:
            self.v = np.random.normal(loc=0, scale=np.sqrt(kb_au * self.T_target / self.mass), size=(self.N, 3))

            np.savetxt(f"{vel_file}", self.v)
             
        self.force = None      

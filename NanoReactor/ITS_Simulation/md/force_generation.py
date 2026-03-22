from ase.units import Bohr, Hartree
import numpy as np
from xtb.ase.calculator import XTB 


def force_xtb(system):
    """xtb-python checked!"""
    # print("Use GFN2-xTB")
    system.Atom.set_positions(system.xyz * Bohr)
    system.Atom.calc = XTB(method="GFN2-xTB", accuracy=1.0, electronic_temperature=300, max_iterations=5000)
    U = system.Atom.get_potential_energy() / Hartree
    force = system.Atom.get_forces() * Bohr / Hartree
    # print(f"force: {force}")
    return U, force

def force_DPA_3_F(system):
    """Not Checked!"""
    from deepmd.pt.utils.ase_calc import DPCalculator
    system.Atom.set_positions(system.xyz * Bohr)
    calc = DPCalculator(model="DPA-3-F@DFT.pt", device='cuda')
    system.Atom.calc = calc
    U = system.Atom.get_potential_energy() / Hartree
    force = system.Atom.get_forces() * Bohr / Hartree
    return U, force

def force_DPA_3_DF(system):
    """Not Checked!"""
    from deepmd_xtb import DP_xTB
    system.Atom.set_positions(system.xyz * Bohr)
    calc = DP_xTB(model="DPA-3-DF@DFT.pt", device='cuda')
    system.Atom.calc = calc
    U = system.Atom.get_potential_energy() / Hartree
    force = system.Atom.get_forces() * Bohr / Hartree
    return U, force

def force_calculation(settings, system, nstep):

    calculator_map = {
            "GFN2-xTB": force_xtb,
            "DPA-3-F": force_DPA_3_F,
            "DPA-3-DF": force_DPA_3_DF
        }

    if settings.calculator_type not in calculator_map:
            raise ValueError(f"Invalid calculator type: {settings.calculator_type}. "
                             f"Valid options: {list(calculator_map.keys())}")

    force_funct = calculator_map[settings.calculator_type]
    U, force_ini = force_funct(system)

    force = force_ini
    U_system = U
        
    if settings.mtd:
        U_mtd, grad_Umtd = settings.domtd.mtd_cal(system, nstep)
        force = force - grad_Umtd
        # print(f"gradient of mtd: {grad_Umtd}")
        U_system = U_system + U_mtd
        settings.domtd.Umtd_list.append(U_mtd)

    if settings.nano:
        U_nano, grad_nano = settings.donano.nano_cal(system)
        force = force -grad_nano
        # print(f"gradient of nanoreactor: {grad_nano}")
        U_system = U_system + U_nano
        settings.donano.Uwall_list.append(U_nano)

    if settings.its:
        if nstep == 0:
            settings.doits.pe_initial = U
        else:
            if settings.doits.pe_initial > U:
                ck_new = np.exp(-settings.doits.beta * (U - settings.doits.pe_initial))
                settings.doits.nk = settings.doits.nk * ck_new
                settings.doits.pe_initial = U
        settings.doits.pe_ref_list.append(settings.doits.pe_initial)
        S, settings.doits.lnPk = settings.doits.its_force(U=U_system, pe_initial=settings.doits.pe_initial, lnPk=settings.doits.lnPk)
        settings.doits.S_list.append(S)
        # print(f"pe initial: {self.pe_initial}")
        # print(f"Pk from force generation:\n {self.ITS.Pk}")
        # print(f"U from force generation(do its):\n {U}")
        # print("Enhance Factor:", S)
        force = S * force

    return force
from settings import SimulationSettings, System
from src.md.integrator import get_integrator
from src.md.force_generation import force_calculation
from src.temperature_control.thermostat import get_thermostat
from src.tool.constant import kb_au
from ase.io import write
from ase.units import Bohr
import os
import matplotlib.pyplot as plt
import numpy as np


def main():

    settings = SimulationSettings()
    md_config = settings.get_module_config("md")
    system = System()

    integrator_type = md_config.get("integrator", 1)
    dt = settings.dt
    total_time = settings.time
    steps = int(total_time / dt)

    system.force = force_calculation(settings, system, nstep=0)
    integrator = get_integrator(integrator_type)
    
    thermostat_config = settings.thermostat

    T_list = []
    if thermostat_config:
        thermostat_type = settings.thermostat.get('type')
        thermostat = get_thermostat(thermostat_type)
        if thermostat_type == 1:
            for _ in range(1, steps):
                thermostat(settings, system)
                integrator(settings, system, dt, _)
                thermostat(settings, system)
                write('md.xyz', system.Atom, format='xyz', append=(_ > 0))

                T = 0.5 * np.sum(system.mass * system.v ** 2) / (1.5 * kb_au * system.N)
                T_list.append(T)

                # with open("coord.xyz", "a") as f:
                #     f.write(str(system.xyz) + '\n')
                
                # with open("velocity.dat", "a") as g:
                #     g.write(str(system.v) + '\n')

                # # Update nk
                if settings.its:
                    if _ % settings.doits.update_step == 0 and _ != 0 and os.path.exists('{dataPath}/nk.dat') is not True:
                        settings.doits.its_update()
                        print(_)

            np.savetxt("T.dat", T_list)
            plt.plot(T_list)
            plt.savefig("T~t.png")
            plt.close()

            if settings.its:
                np.savetxt("lnnk.dat", settings.doits.lnnk_list)
                np.savetxt("S.dat", settings.doits.S_list)
                np.savetxt("U_ref.dat", settings.doits.pe_ref_list)
                plt.plot(settings.doits.lnnk_list)
                plt.savefig("lnnk.png")
                plt.close()
                plt.plot(settings.doits.S_list)
                plt.savefig("Enhanced_Factor.png")
                plt.close()

            if settings.nano:
                np.savetxt(fname="Uwall.dat", X=settings.donano.Uwall_list)
            
            if settings.mtd:
                np.savetxt(fname="Umtd.dat", X=settings.domtd.Umtd_list)
                np.savetxt(fname="cv_ref.dat", X=settings.domtd.cv_ref_list)
                for ii in range(len(settings.domtd.x_ref_list)):
                    np.savetxt(fname=f"{ii}_ref.dat", X=settings.domtd.x_ref_list[ii])

        elif thermostat_type == 2:
            for _ in range(steps):
                integrator(settings, system, dt, _)
                thermostat(settings, system, dt)
                print(_)
    else:
        for _ in range(steps):
            integrator(settings, system, dt, _)
            print(_)
    

if __name__ == "__main__":
    main()

import numpy as np

def Langevin(settings, system):
    r"""
    Langevin Dynamics:
        Reference:
            Bussi, G.; Parrinello, M. Accurate sampling using Langevin dynamics. Phys. Rev. E 2007, 75, 056707.

        How to Update?
            c1 = exp(-\gamma * dt / 2)
            c2 = \sqrt((1 - c1 ** 2) * mass * kb * T_target)

            p = c1 * p + c2 * np.random.standard_normal(size=(self.N, 3))

            xyz = xyz + p * dt / mass - 0.5 * grad * dt ** 2 / mass
            p = p - (grad + grad_new) * dt / 2

            p = c1 * p + c2 * np.random.standard_normal(size=(self.N, 3))
    """
    # print("Langevin Thermostat is ON.")
    system.v = settings.c1 * system.v + settings.c2 * np.random.standard_normal(size=system.v.shape) / np.sqrt(system.mass)
    # system.v = settings.c1 * system.v + settings.c2 * np.random.standard_normal() / np.sqrt(system.mass)        # For HO test

def Berendsen(settings, system, dt):
    T = 0.5 * np.sum(system.mass * system.v ** 2) / 1.5
    factor = np.sqrt(1. + (settings.T_target - T) * dt / (settings.tau * T))
    system.v = factor * system.v

def get_thermostat(thermostat_type: int):
    thermostats = {
        1: Langevin,
        2: Berendsen
    }
    return thermostats.get(thermostat_type, Langevin)
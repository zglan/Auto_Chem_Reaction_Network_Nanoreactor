from .force_generation import force_calculation
from ase.units import Bohr

def velocity_verlet(settings, system, dt, nstep):
    system.xyz += dt * system.v + 0.5 * dt**2 * system.force / system.mass
    system.Atom.set_positions(system.xyz * Bohr)
    f_new = force_calculation(settings, system, nstep)
    system.v += 0.5 * dt * (system.force + f_new) / system.mass
    system.force = f_new

def leapfrog(settings, system, dt, steps):
    exit()

def get_integrator(integrator_type: int):
    integrators = {
        1: velocity_verlet,
        2: leapfrog
    }
    return integrators.get(integrator_type, velocity_verlet)
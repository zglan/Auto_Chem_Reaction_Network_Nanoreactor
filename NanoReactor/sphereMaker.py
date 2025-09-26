#!/usr/bin/env python3

import os
import math
import random
import argparse

import numpy as np
from rdkit import Chem



def mole_getmass(f):
    """AMU mass""" 
    mass = 0
    with open(f, 'r') as xyz:   
        for idx, atom in enumerate(xyz):   
            mass += Chem.Atom(atom.split()[0]).GetMass() if idx > 1 else 0 
    return mass


def system_getradii(mass, density, length_unit='A'):
    """
    Get system radii for the input mass and density.

    Parameters
    ----------
    mass: unit in AMU.
    density: unit in g/ml.
    length unit: A or bohr.

    Returns
    -------
    system radii: unit in A or bohr.

    """
    NA = 6.02e23
    BohrToA = 0.52917720859
    radii = ((mass / NA / density) / 3*4 / 3.14)**(1 / 3) / 1e-8
    radii = radii if length_unit == 'A' else radii / BohrToA
    return round(radii, 1)



def packmol_make(density, nums, xyzfiles, name, center, dis=2):
    mole_mass = np.array([mole_getmass(f) for f in xyzfiles])
    mass = np.sum(np.array(nums)*mole_mass)
    radii = system_getradii(mass, density)
    print('*** Sphere Radii is %.2f A.' %(radii))
    print('*** Sphere Radii is %.2f Bohr.' %(radii / 0.52918))
    with open('%s.inp'%(name), 'w') as inp:
        inp.write(
            'tolerance %s\n' % (dis) +
            'seed -1\n' +
            'randominitialpoint\n' +
            'filetype xyz\n' +
            'output %s.xyz\n'% (name))
        for num, file in zip(nums, xyzfiles):
            if center == False:
                inp.write(
                    'structure %s\n'%(file) + 
                    '  number %d\n'%(num) + 
                    '  inside sphere %.1f %.1f %.1f %.1f\n'%(0, 0, 0, radii) + 
                    'end structure\n')
            else:
                if file not in center:
                    inp.write(
                        'structure %s\n'%(file) + 
                        '  number %d\n'%(num) + 
                        '  inside sphere %.1f %.1f %.1f %.1f\n'%(0, 0, 0, radii) + 
                        'end structure\n')
                else:
                    inp.write(
                        'structure %s\n'%(file) + 
                        '  number %d\n'%(num) + 
                        '  center\n' +
                        '  fixed 0. 0. 0. %f %f %f\n' % (
                            math.radians(random.randint(0, 180)), 
                            math.radians(random.randint(0, 180)), 
                            math.radians(random.randint(0, 180))) + 
                        'end structure\n') 
                    
    try:
        os.system('packmol < %s > packmol.log'%('%s.inp'%(name)))
        print('=== %s.xyz is built' %(name))
    except:
        raise Exception('Can not use packmol, failed to generate')


def xtb_mtd(radii, xyzf, t, T, dt, k):
    with open('input.inp', 'w') as inp:
        inp.write(
        '$md\n' +
        '   time=%d\n'%(t) +
        '   step=%.1f\n'%(dt) + 
        '   temp=%d\n'%(T) + 
        '   dump=1\n' + 
        '   shake=0\n' + 
        '$end\n' + 
        '$metadyn\n' + 
        '   save=10\n' + 
        '   kpush=%.2f\n'%(k) + 
        '   alp=0.7\n' + 
        '$end\n' + 
        '$wall\n' + 
        '   potential=logfermi\n' + 
        '   beta=10.0\n' + 
        '   temp=6000.0\n' + 
        '   sphere: %.1f, all\n'%(radii) + 
        '$end\n')

    try:
        os.system('xtb --md --i %f %f'%(inp, xyzf))
    except:
        raise Exception


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument(
    '-i', '--inputfile', nargs='*',
    help='Input mole file, e.g. water.xyz OH.xyz', required=True)
parser.add_argument(
    '-n', '--nums', nargs='*',
    help='Input mole nums, e.g. 10 1', type=int, required=True)
parser.add_argument(
    '-d', '--density',
    help='Set system density, unit in g/ml, e.g. 10', type=float, required=True)
parser.add_argument(
    '--names',
    help='Set Output file name', nargs='*', required=True)
parser.add_argument(
    '--center',
    help='Set Central Molecules', nargs='*', required=False, default=False)
parser.add_argument(
    '--dis',
    help='Set Mol distance', type=float, default=2)  
args = parser.parse_args()


inputfile = args.inputfile
nums = args.nums
density = args.density
names = args.names
dis = args.dis
center = args.center


for name in names:
    packmol_make(density, nums, inputfile, name, center, dis)
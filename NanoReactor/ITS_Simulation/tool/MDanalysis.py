#!/share/home/lijie/anaconda3/bin/python

import numpy as np
import argparse
import matplotlib.pyplot as plt

"""
Only .xyz file is supported !
"""

def read_xyz(filename="position.xyz"):
    with open(f"{filename}", "r") as f:
        """
        Natom: Number of Atoms, type: int
        atom_list: Name of Atoms, type: list(str)
        xyz: xyz of Atoms, type: array(float)
        """
        lines = f.readlines()
        Natom = int(lines[0])
        
        xyz = []
        for i, _ in enumerate(lines):
            # print(i)
            if (i - 2) % (Natom + 2) == 0:
                # print(lines[i: i+Natom])
                temp = []
                for line in lines[i: i+Natom]:
                    temp.append(np.array(line.split()[1:], dtype=float))
                xyz.append(np.array(temp, dtype=float))
        xyz = np.array(xyz)
        return xyz

def cal_bond(xyz: np.array, 
             atom_1: int, 
             atom_2: int) -> np.array:
    """
    xyz: (shot_num, atom_num, dim_num)
    atom_1: int
    atom_2: int
    """
    A = xyz[:, atom_1, :]
    B = xyz[:, atom_2, :]
    AB = A - B
    # print(AB.shape)
    distance = np.linalg.norm(AB, axis=1)
    # print(distance)
    return distance

def cal_angle(xyz: np.array, 
             atom_1: int, 
             atom_2: int,
             atom_3: int) -> np.array:
    A = xyz[:, atom_1, :]
    B = xyz[:, atom_2, :]
    C = xyz[:, atom_3, :]
    AB = A - B
    AC = A - C
    AB_norm = np.linalg.norm(AB, axis=1)
    AC_norm = np.linalg.norm(AC, axis=1)
    AB_dot_AC = np.einsum("ij,ij->i", AB, AC)
    angle = np.degrees(np.arccos(AB_dot_AC) / (AB_norm * AC_norm))
    return angle

def cal_dihedral(xyz: np.array, 
                atom_1: int, 
                atom_2: int,
                atom_3: int,
                atom_4: int):
    A = xyz[:, atom_1, :]
    B = xyz[:, atom_2, :]
    C = xyz[:, atom_3, :]
    D = xyz[:, atom_4, :]
    AB = B - A
    BC = C - B
    CD = D - C
    normal1 = np.cross(AB, BC)
    normal1 = normal1 / np.linalg.norm(normal1, axis=1, keepdims=True)
    normal2 = np.cross(BC, CD)
    normal2 = normal2 / np.linalg.norm(normal2, axis=1, keepdims=True)
    dihedral = np.einsum("ij, ij->i", normal1, normal2)

    normal_cross = np.cross(normal1, normal2)
    sign = np.sign(np.einsum("ij, ij->i", normal_cross, BC))
    print(sign)
    dihedral = dihedral * sign
    return dihedral

def evolution(filename, data):
    plt.plot(data)
    plt.title(f'Data Visualization: Evolution of {filename}')
    plt.xlabel('The Number of Step(N)')
    plt.ylabel(f'The Value of {filename}')
    plt.grid(True)
    plt.savefig(f"{filename}-evolution", dpi=150, bbox_inches='tight')
    plt.close()

def distribution(filename, data, bin_num=500):
    data_max, data_min = data.max(), data.min()
    bin_length = (data_max - data_min) / bin_num
    hist, bin_edges = np.histogram(data, bin_num)
    hist = hist / (len(data) * bin_length)

    np.savetxt(
        f"{filename}-prob_density.dat",          
        np.column_stack((bin_edges[:-1], hist)), 
        fmt='%.6f %.6f'                          
    )

    plt.title(f'Data Visualization: Distribution of {filename}')
    plt.xlabel(f'The Value of {filename}')
    plt.ylabel(f'The Probability Density of {filename}')
    plt.bar(bin_edges[:-1], hist, width=bin_length)
    plt.grid(True)
    plt.savefig(f"{filename}-distribution", dpi=150, bbox_inches='tight')
    plt.close()

def main():
    # Parameter Provide:
    parser = argparse.ArgumentParser(description="Calculation of bond, angle, dihedral and others for MD trajectory file.")
    parser.add_argument("--input", help="Trajectory Filename(eg: position.xyz)")
    parser.add_argument("-b", "--bond", nargs=2, type=int, metavar=("A1", "A2"), action="append",
                      help="Bond Calculation (Atomic Ids: 0, 1, ...)")
    parser.add_argument("-a", "--angle", nargs=3, type=int, metavar=("A1", "A2", "A3"), action="append",
                      help="Angle Calculation (Atomic Ids: 0, 1, ...), Notes that A1 is the center atom. eg: angle H-O-H, \
                           Assume atom order is O H H. For angle H-O-H, -a 0 1 2 or -a 0 2 1. For angle H-H-O, -a 1 2 0 or -a 2 1 0")
    parser.add_argument("-d", "--dihedral", nargs=4, type=int, metavar=("A1", "A2", "A3", "A4"), action="append",
                      help="Dihedral Calculation (Atomic Ids: 0, 1, ...)")
    # parser.add_argument("--draw",  nargs='*', type=str, action="append",
    #                   help="Draw type, Evolution(--draw t) or Distribution(--draw p) or both(--draw t p)")
    args = parser.parse_args()

    # Coordinate Reading:
    xyz = read_xyz(filename=args.input)
    # print(xyz.shape)

    # Job Process:
    if args.bond: 
        # print(args.bond)
        for bond in args.bond:
            # print(bond)
            outfile = f"bond_{bond[0]}_{bond[1]}"
            distance = cal_bond(xyz, bond[0], bond[1])
            np.savetxt(fname=outfile, X=distance)
            evolution(outfile, distance)
            distribution(outfile, distance)

    if args.angle: 
        # print(args.angle)
        for angle in args.angle:
            # print(angle)
            outfile = f"angle_{angle[0]}_{angle[1]}_{angle[2]}"
            angle = cal_angle(xyz, angle[0], angle[1], angle[2])
            np.savetxt(fname=outfile, X=angle)
            evolution(outfile, angle)
            distribution(outfile, angle)

    if args.dihedral: 
        # print(args.dihedral)
        for dihedral in args.dihedral:
            # print(dihedral)
            outfile = f"dihedral_{dihedral[0]}_{dihedral[1]}_{dihedral[2]}_{dihedral[3]}"
            dihedral = cal_dihedral(xyz, dihedral[0], dihedral[1], dihedral[2], dihedral[3])
            np.savetxt(fname=outfile, X=dihedral)
            evolution(outfile, dihedral)
            distribution(outfile, dihedral)
    


if __name__ == '__main__':
    main()

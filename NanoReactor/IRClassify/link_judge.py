import numpy as np
import networkx as nx


# Covalent radii from Cordero et al. 'Covalent radii revisited' Dalton Transactions 2008, 2832-2838.
Radii = [0.31, 0.28, # H and He
         1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, # First row elements
         0.00, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, # Second row elements
         1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, # Second row elements
         2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.61, 1.52, 1.50,
         1.24, 1.32, 1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16, # Third row elements, K through Kr
         2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42,
         1.39, 1.45, 1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.40, # Fourth row elements, Rb through Xe
         2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98,
         1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87, # Fifth row elements, s and f blocks
         1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36,
         1.36, 1.32, 1.45, 1.46, 1.48, 1.40, 1.50, 1.50, # Fifth row elements, d and p blocks
         2.60, 2.21, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.69] # Sixth row elements

Elements = ['H','He',
            'Li','Be','B','C','N','O','F','Ne',
            'Na','Mg','Al','Si','P','S','Cl','Ar',
            'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
            'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',
            'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
            'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',
            'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt']


def linkJudge(coord, muti=1.1):
    
    def coord2RadiiMatrix(coord, muti):
        radii_matr = np.zeros((length, length))
        radii_list = []
        for c in coord:
            for idx, e in enumerate(Elements):
                if c[0] == e:
                    radii_list.append(Radii[idx])

        for i in range(length):
            for j in range(i + 1, length):
                radii_matr[i, j] = (radii_list[i] + radii_list[j]) * muti
        
        return radii_matr

    def coord2distMatrix(coord):
        dist_matr = np.zeros((length, length))
        for i in range(length):
            for j in range(i + 1, length):
                atomi = np.array(coord[i][1:])
                atomj = np.array(coord[j][1:])
                dist_matr[i, j] = np.linalg.norm(atomi - atomj)
        return dist_matr
    
    length = len(coord)
    link_list = []
    dist_matr = coord2distMatrix(coord)
    radii_matr = coord2RadiiMatrix(coord, muti)

    for i in range(length):
        for j in range(i + 1, length):
            if dist_matr[i, j] <= radii_matr[i, j]:
                link_list.append((i, j))

    G = nx.Graph()
    G.add_edges_from(link_list)
    G.add_nodes_from(list(range(length)))
    con = [sorted(list(G.subgraph(c).nodes)) for c in nx.connected_components(G)]

    return con, G.edges
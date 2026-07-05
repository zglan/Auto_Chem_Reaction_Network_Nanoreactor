import numpy as np
# import matplotlib.pyplot as plt

class Nanoreactor:
    def __init__(self,
                 R: float = 10.,
                 k_wall: float = 0.019,       # a.u.
                 beta_wall: float = 10.):           # a.u.
        self.R = R
        self.k_wall = k_wall
        self.beta_wall = beta_wall

        self.Uwall_list = []
    
    def mass_center(self, system):
        r"""
        :param xyz:  array, (N, 3)
        :param mass:  array, (N, 1)
        :return: center of mass:  array, (3, )
        """
        # x = np.reshape(system.xyz[:, 0], (len(system.mass), 1))
        # y = np.reshape(system.xyz[:, 1], (len(system.mass), 1))
        # z = np.reshape(system.xyz[:, 2], (len(system.mass), 1))
        # x_mass = np.sum(system.mass * x) / np.sum(system.mass)
        # y_mass = np.sum(system.mass * y) / np.sum(system.mass)
        # z_mass = np.sum(system.mass * z) / np.sum(system.mass)
        # return np.array([x_mass, y_mass, z_mass])
        return np.array([0., 0., 0.])           

    def distance(self, xyz, xyz_mass):
        """
        :param xyz_mass:  array, (N, 3)
        :param xyz:  (3, )
        :return: distance, array, (N, 1)
        """
        r = []
        for i in xyz:
            temp = np.linalg.norm(i - xyz_mass)
            r.append(temp)
        r = np.array(r)
        r = np.reshape(r, (len(r), 1))
        return r

    def wall_pe(self, r: float):
        """
        R: Radius of Sphere, float
        r: Distance between Mass-Center and Every Atom, np.array

        in a.u.
        """
        temp = np.log(1 + np.exp(-self.beta_wall * (self.R - r)))
        E_wall = self.k_wall * np.sum(temp)
        return E_wall

    def nano_cal(self, system):
        r"""
        E_wall = k * \sum_{i}^{N}{log(1 + exp(-beta * (R - ri)))}
        dE_wall / dqi = (k * log(1 + exp(-beta * (R - ri))))'
                      = (k / (1 + exp(-beta * (R - ri)))) * exp(-beta * (R - ri)) * beta * (qi - q_mass) / ri
        :return:
        """
        xyz_mass = self.mass_center(system)
        r = self.distance(xyz=system.xyz, xyz_mass=xyz_mass)
        # wall pe:
        temp = np.log(1 + np.exp(-self.beta_wall * (self.R - r)))
        E_wall = self.k_wall * np.sum(temp)

        # wall grad:
        displacement = (system.xyz - xyz_mass) / r
        grad_wall = self.k_wall / (1 + np.exp(-self.beta_wall * (self.R - r)))
        grad_wall = grad_wall * np.exp(-self.beta_wall * (self.R - r))
        grad_wall = grad_wall * self.beta_wall
        grad_wall = grad_wall * displacement
        return E_wall, grad_wall
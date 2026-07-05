import numpy as np
from .cvs_cal import CVs_Calculation

class MTD:
    def __init__(self,
                 cvs_type: str,                  # eg: 'bond'
                 atom_ids: list = None,          # eg: [0, 4]
                 kpush: float = 0.2,             # a.u.
                 sigma: float = 0.7,             # a.u.
                 kappa: float = 0.03,
                 dump: int = 500):               # a.u.
        self.cvs_type = cvs_type
        self.atom_ids = atom_ids
        self.kpush = kpush
        self.sigma = sigma
        self.kappa = kappa
        self.dump = dump

        self.counter = []
        self.cv_ref_list = []
        self.x_ref_list = []
        self.Umtd_list = []

    
    def factor_cal(self, 
                   i: int,
                   xyz: np.array):
        """
        Input:
            i: ith step of simulation.
            xyz: (N, D), np.array.
        Output:
            self.counter: (I, ), np.array
        """
        if i % self.dump == 0:
            if i != 0 or not self.x_ref_list:
                self.counter = list(self.counter)
                self.counter.append(0)
                self.x_ref_list.append(xyz.copy())
                print(f"x_ref: {self.x_ref_list}")
                
        
        self.counter = np.array(self.counter)
        self.counter = self.counter + np.ones_like(self.counter)
        return 2 / (1 + np.exp(-self.kappa * self.counter)) - 1
    

    """CVs Calculation Control"""
    @staticmethod
    def cv_cal(self,
                system,
                cvs_type: str,
                atom_ids: list):
        # print(cvs_type)
        if cvs_type in ['Bond', 'bond', 'b', 'B']:
            # print(cvs_type)
            cv = CVs_Calculation.cal_bond(system.xyz, atom_ids[0], atom_ids[1])

        elif cvs_type in ['Angle', 'angle', 'a', 'A']:
            # print(cvs_type)
            cv = CVs_Calculation.cal_angle(system.xyz, atom_ids[0], atom_ids[1], atom_ids[2])
        
        elif cvs_type in ['Dihedral', 'dihedral', 'd', 'D']:
            # print(cvs_type)
            cv = CVs_Calculation.cal_dihedral(system.xyz, atom_ids[0], atom_ids[1], atom_ids[2], atom_ids[3])
        
        elif cvs_type in ['RMSD', 'rmsd', 'r', 'R']:
            # print(cvs_type)
            x_ref = np.array(self.x_ref_list)
            cv = CVs_Calculation.cal_rmsd(system.xyz, x_ref)
        return cv
    
    @staticmethod
    def grad_cv_cal(self,
                    system,
                    cvs_type: str,
                    atom_ids: list):
        if cvs_type in ['Bond', 'bond', 'b', 'B']:
            # print(cvs_type)
            grad_cv = CVs_Calculation.grad_bond_cal(system.xyz, atom_ids[0], atom_ids[1])
        elif cvs_type in ['Angle', 'angle', 'a', 'A']:
            # print(cvs_type)
            pass
        elif cvs_type in ['Dihedral', 'dihedral', 'd', 'D']:
            # print(cvs_type)
            pass
        elif cvs_type in ['RMSD', 'rmsd', 'r', 'R']:
            x_ref = np.array(self.x_ref_list)             # May wrong?
            grad_cv = CVs_Calculation.grad_rmsd_cal(system.xyz, x_ref)
        return grad_cv


    def mtd_cal(self,
                system,
                i: int) -> np.array:
        if self.cvs_type in ['RMSD', 'rmsd', 'r', 'R']:
            """
            k exp (-rmsd ** 2 / (2 * sigma ** 2)) * (x_shift - y_rot) / (N * sigma ** 2)
            """
            r"""
            U_mtd = \sum_{i=1}^{n} {k_i exp(-\alpha * CV ** 2)}
            where:
                CV: RMSD, 
                    CV = \sqrt{
                            \sum_{j=1}^N (r_j - r_j^{ref, i}) ** 2 / N
                            }
                    calculate by quaternion algorithm of Coutsias.
                    (Journal of Computational Chemistry (2004), 25 (15), 1849-1857CODEN: JCCHDD; ISSN:0192-8651. (John Wiley & Sons, Inc.))
                r: coordinate of current structure
                r_ref: coordinate of reference structure
                N: number of atoms
                i: ith reference structure
            
            Instantaneous addition of a biasing potential in eq 2 can cause instabilities in the MD with large time steps, a factor for U_mtd is needed:
                f_dump = \dfrac{2}{1 + exp(- \kappa * k)}
                where:
                    \kappa = 0.03
                    k: MD step counter, eg: when new U_mtd dump, k = 0, and independent.
            """
            """
            Input:
                x: (N, D), N: atom number; D: dimension
                x_ref: (I, N, D), I: reference structure number
                counter: (I, )
            Output:
                U_mtd: float
            """
            if i % self.dump == 0:
                self.x_ref_list.append(system.xyz)

            rmsd_val = MTD.cv_cal(self,
                        system=system,
                        cvs_type=self.cvs_type,
                        atom_ids=self.atom_ids)
            # print(rmsd_val)
            grad_rmsd = MTD.grad_cv_cal(self,
                                system=system,
                                cvs_type=self.cvs_type,
                                atom_ids=self.atom_ids)
            
            U_mtd = np.sum(self.kpush * np.exp(-rmsd_val ** 2 / (2 * self.sigma ** 2)))
            self.Umtd_list.append(U_mtd)

            factor = np.exp(-rmsd_val ** 2 / (2 * self.sigma ** 2)) / (self.sigma ** 2)
            factor = np.expand_dims(factor, axis=-1)
            factor = np.expand_dims(factor, axis=-1)
            grad_Umtd = np.sum(factor * grad_rmsd, axis=0)

        else:
            # print(f"CVs type: {self.cvs_type}")
            cv = MTD.cv_cal(self,
                        system=system,
                        cvs_type=self.cvs_type,
                        atom_ids=self.atom_ids)
            # print(f"cv value: {cv}")
            if i % self.dump == 0:
                self.cv_ref_list.append(cv)
            # print(f"cv_ref value: {self.cv_ref_list}")
            grad_cv = MTD.grad_cv_cal(self,
                                system=system,
                                cvs_type=self.cvs_type,
                                atom_ids=self.atom_ids)
            
            cv_ref = np.array(self.cv_ref_list)
            # print(f"grad cv: {grad_cv}")
            # print(f"k value: {self.kpush}")
            # print(f"exp (delta cv ** 2 / (2 * sigma ** 2)): {np.exp(-(cv - cv_ref) ** 2 / (2 * self.sigma ** 2))}")
            U_mtd = np.sum(self.kpush * np.exp(-(cv - cv_ref) ** 2 / (2 * self.sigma ** 2)))
            self.Umtd_list.append(U_mtd)

            factor = np.exp(-(cv - cv_ref) ** 2 / (2 * self.sigma ** 2)) * \
                        (-(cv - cv_ref)) / (self.sigma ** 2)
            # print(f"dUmtd / ds:  {factor}")
            factor = np.expand_dims(factor, axis=-1)
            factor = np.expand_dims(factor, axis=-1)

            grad_cv = np.expand_dims(grad_cv, axis=0)

            # grad_Umtd = np.sum(self.kpush * np.exp(-(cv - cv_ref) ** 2 / (2 * self.sigma ** 2)) * \
            #               (-(cv - cv_ref) / (self.sigma ** 2)) * grad_cv)
            grad_Umtd = np.sum(factor * grad_cv, axis=0)
        # print(f"U_mtd: {U_mtd}")
        # print(f"grad_mtd:\n {grad_Umtd}")
        print("-----------------------------------------------------------------------------------------")
        return U_mtd, grad_Umtd

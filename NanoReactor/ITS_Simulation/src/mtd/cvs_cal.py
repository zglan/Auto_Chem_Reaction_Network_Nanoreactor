import numpy as np


class CVs_Calculation:

    r"""
    Reference:
        Coutsias, E. A.; Seok, C.; Dill, K. A. Using Quaternions to Calculate RMSD. 
        J. Comput. Chem. 2004, 25, 1849-1857,  DOI: 10.1002/jcc.20110

    How to implement this code?
        1. calculate barycenter x0 and generate x_shift by x - x0
        2. calculate correlation matrix R_matrix
        3. calculate F_matrix by R_matrix
        4. obtain eigenvalue \lambda_max and its corresponding eigenvector q(quarternion vector)
        5. calculate best-fit RMSD by q, norms of x_shift and y_shift
        6. obtain rotation matrix U from \lambda_max and q
        7. calculate the gradient of RMSD by U, RMSD
    """
    """RMSD Calculation"""

    @staticmethod
    def barycenter_cal(x: np.array):
        r"""
        Barycenter:
                q_center = \sum_i^{N} {q_i} / N

        x: (I, N, D) array
        ::returns::
        x_center: (I, D, ) array
        """
        x_center = np.mean(x, axis=1)
        return x_center

    @staticmethod
    def coordinate_shift(x: np.array, x_center: np.array):
        r"""
        x: (I, N, D) array
        x_center: (I, D, ) array
        ::returns::
        x_shift: (I, N, D) array
        """
        return x - x_center[:, np.newaxis, :]      # x_center[:, np.newaxis, :]: (I, D, ) --> (I, 1, D)

    @staticmethod
    def makeR(x: np.array, x_ref: np.array):
        r"""
        Correlation Matrix R_matrix:
                R_ij = \sum_{k=1}^N {x_ik y_jk} i, j = 1, 2, 3
        x: (1, N, D) array
        x_ref: (I, N, D) array
        ::returns::
        R_matirx: (I, D, D) array
        """
        return np.einsum("mni,lnj-> lij", x, x_ref)

    @staticmethod
    def makeF(R_matrix: np.array):
        r"""
        F_matrix:
        
        R_matirx: (D, D) array
        ::returns::
        F_matrix: (4, 4) array. Hermitian Matrix.
        """
        """
        May have some mistakes because of R_matrix. For the difference of index provided by fortran and python.
        """
        I = R_matrix.shape[0]
        F_matrix = np.zeros((I, 4, 4))

        F_matrix[:, 0, 0] = R_matrix[:, 0, 0] + R_matrix[:, 1, 1] + R_matrix[:, 2, 2]
        F_matrix[:, 1, 0] = F_matrix[:, 0, 1] = R_matrix[:, 1, 2] - R_matrix[:, 2, 1]
        F_matrix[:, 2, 0] = F_matrix[:, 0, 2] = R_matrix[:, 2, 0] - R_matrix[:, 0, 2]
        F_matrix[:, 3, 0] = F_matrix[:, 0, 3] = R_matrix[:, 0, 1] - R_matrix[:, 1, 0]

        F_matrix[:, 1, 1] = R_matrix[:, 0, 0] - R_matrix[:, 1, 1] - R_matrix[:, 2, 2]
        F_matrix[:, 1, 2] = F_matrix[:, 2, 1] = R_matrix[:, 0, 1] + R_matrix[:, 1, 0]
        F_matrix[:, 1, 3] = F_matrix[:, 3, 1] = R_matrix[:, 0, 2] + R_matrix[:, 2, 0]

        F_matrix[:, 2, 2] = -R_matrix[:, 0, 0] + R_matrix[:, 1, 1] - R_matrix[:, 2, 2]
        F_matrix[:, 2, 3] = F_matrix[:, 3, 2] = R_matrix[:, 1, 2] + R_matrix[:, 2, 1]

        F_matrix[:, 3, 3] = -R_matrix[:, 0, 0] - R_matrix[:, 1, 1] + R_matrix[:, 2, 2]
        return F_matrix

    @staticmethod
    def eigen_cal(F_matrix: np.array):
        r"""
            Calculate eigenvalues and eigenvectors of F_matrix, 
            and take the maximum eigenvalue \lambda_max and the corresponding eigenvector q.
        
        F_matrix: (I, 4, 4) array. Hermitian Matrix.
        ::returns::
        eigenvalues: (I, ) array
        eigenvectors: (I, 4) array
        """
        val, vector = np.linalg.eigh(F_matrix)
        max_idx = np.argmax(val, axis=1)
        lambda_max = np.take_along_axis(val, max_idx[:, np.newaxis], axis=1)
        lambda_max = lambda_max.squeeze(axis=1)
        q = np.take_along_axis(vector, max_idx[:, np.newaxis, np.newaxis], axis=2)
        q = q.squeeze(axis=2)
        return lambda_max, q

    @staticmethod
    def quarternion_rmsd(x_shift: np.array, 
                        x_ref_shift: np.array,
                        lambda_max: float):
        r"""
        RMSD:
            \sqrt{
            [\sum_{k=1}^N {(|x_shift,k|^2 + |x_ref_shift,k|^2)} - 2 * \lambda_max] / N
            }
        
        x_shift: (1, N, D) array
        x_ref_shift: (I, N, D) array
        lambda_max: eigenvalues of F_matrix, (I, 1, ) array
        ::returns::
        rmsd: (I, )
        """
        N = x_shift.shape[1]
        x_norm = np.sum(x_shift * x_shift, axis=(1, 2))
        y_norm = np.sum(x_ref_shift * x_ref_shift, axis=(1, 2))
        value = np.maximum(0, (x_norm + y_norm - 2 * lambda_max) / N)
        return np.sqrt(value)

    @staticmethod
    def rotation_matrix(q):
        r"""
        q: eigenvectors of F_matrix, (I, 4) array
        ::return::
        U: rotation matrix, (I, 3, 3) array
        """
        q0, q1, q2, q3 = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
        U = np.array([
            [
                q0**2 + q1**2 - q2**2 - q3**2,
                2*(q1*q2 - q0*q3),
                2*(q1*q3 + q0*q2)
            ],
            [
                2*(q1*q2 + q0*q3),
                q0**2 - q1**2 + q2**2 - q3**2,
                2*(q2*q3 - q0*q1)
            ],
            [
                2*(q1*q3 - q0*q2),
                2*(q2*q3 + q0*q1),
                q0**2 - q1**2 - q2**2 + q3**2
            ]
        ])
        U = np.moveaxis(U, -1, 0)
        return U

    @staticmethod
    def grad_RMSD_cal(U: np.array, 
                x_shift: np.array, 
                x_ref_shift: np.array, 
                rmsd: float):
        r"""
            Calculate the gradient of best-fit rmsd.
        U: rotation matrix, (I, 3, 3) array
        x_shift: (1, N, D) array
        x_ref_shift: (I, N, D) array
        rmsd: (I, )
        ::return::
        grad: (I, N, D) array
        """
        N = x_shift.shape[1]
        y_rotated = np.einsum("ijk,ink->ijn", U, x_ref_shift)
        # grad = (x_shift - np.transpose(y_rotated, (0, 2, 1))) / (rmsd[:, np.newaxis, np.newaxis] * N)  # (N, 3)
        # return grad
        return (x_shift - np.transpose(y_rotated, (0, 2, 1))) / N

    @staticmethod
    def rmsd_main(coord1, coord2):
        """
        coord1: (1, N, D)
        coord2: (I, N, D), I >= 1
        ::return::
        rmsd: (I, )
        grad_rmsd: (I, N, D)
        """
        x1_center = CVs_Calculation.barycenter_cal(coord1)
        x2_center = CVs_Calculation.barycenter_cal(coord2)

        x1_shift = CVs_Calculation.coordinate_shift(coord1, x1_center)
        x2_shift = CVs_Calculation.coordinate_shift(coord2, x2_center)

        R_matrix = CVs_Calculation.makeR(x1_shift, x2_shift)
        F_matrix = CVs_Calculation.makeF(R_matrix)

        lambda_max, q = CVs_Calculation.eigen_cal(F_matrix)
        U = CVs_Calculation.rotation_matrix(q)

        rmsd = CVs_Calculation.quarternion_rmsd(x1_shift, x2_shift, lambda_max)
        grad_rmsd = CVs_Calculation.grad_RMSD_cal(U, x1_shift, x2_shift, rmsd)
        return rmsd, grad_rmsd

    """CVs Calculation"""
    @staticmethod
    def cal_bond(xyz: np.array, 
             atom_1: int, 
             atom_2: int) -> np.array:
        """
        xyz: (shot_num, atom_num, dim_num)
        atom_1: int
        atom_2: int
        """
        A = xyz[atom_1, :]
        B = xyz[atom_2, :]
        AB = A - B
        # print(AB.shape)
        distance = np.linalg.norm(AB, axis=0)
        # print(distance)
        return distance
    
    @staticmethod
    def cal_angle(xyz: np.array, 
                atom_1: int, 
                atom_2: int,
                atom_3: int) -> np.array:
        A = xyz[atom_1, :]
        B = xyz[atom_2, :]
        C = xyz[atom_3, :]
        AB = A - B
        AC = A - C
        AB_norm = np.linalg.norm(AB, axis=0)
        AC_norm = np.linalg.norm(AC, axis=0)
        AB_dot_AC = np.einsum("ij,ij->i", AB, AC)
        angle = np.degrees(np.arccos(AB_dot_AC) / (AB_norm * AC_norm))
        return angle

    @staticmethod
    def cal_dihedral(xyz: np.array, 
                    atom_1: int, 
                    atom_2: int,
                    atom_3: int,
                    atom_4: int):
        A = xyz[atom_1, :]
        B = xyz[atom_2, :]
        C = xyz[atom_3, :]
        D = xyz[atom_4, :]
        AB = B - A
        BC = C - B
        CD = D - C
        normal1 = np.cross(AB, BC)
        normal1 = normal1 / np.linalg.norm(normal1, axis=0, keepdims=True)
        normal2 = np.cross(BC, CD)
        normal2 = normal2 / np.linalg.norm(normal2, axis=0, keepdims=True)
        dihedral = np.einsum("ij, ij->i", normal1, normal2)

        normal_cross = np.cross(normal1, normal2)
        sign = np.sign(np.einsum("ij, ij->i", normal_cross, BC))
        # print(sign)
        dihedral = dihedral * sign
        return dihedral


    @staticmethod    
    def cal_rmsd(xyz: np.array,
                 xyz_ref: np.array):
        """
        xyz: (atom_num, dim_num)
        xyz_ref: (atom_num, dim_num)
        """
        if len(xyz.shape) == 1:
            raise ValueError("The shape of xyz is wrong!")

        if len(xyz.shape) == 2:
            """xyz: (N, D)"""
            xyz = xyz[np.newaxis, :, :]

        rmsd_val, grad_rmsd = CVs_Calculation.rmsd_main(coord1=xyz,
                               coord2=xyz_ref)
        return rmsd_val
        

    """Gradient of CVs Calculation"""
    @staticmethod
    def grad_bond_cal(xyz: np.array, 
             atom_1: int, 
             atom_2: int) -> np.array:
        """
        xyz: (atom_num, dim_num)
        atom_1: int
        atom_2: int
        """
        print(f"xyz: {xyz}")
        A = xyz[atom_1, :]
        B = xyz[atom_2, :]
        AB = A - B
        print(f"delta q: {AB}")
        distance = np.linalg.norm(AB, axis=0)
        # print(distance)
        temp = (A - B) / distance
        print(f"deduce: {temp}")
        grad_bond = np.zeros_like(xyz)
        grad_bond[atom_1] = temp
        grad_bond[atom_2] = -temp
        return grad_bond
    
    @staticmethod
    def grad_angle_cal(xyz: np.array, 
             atom_1: int, 
             atom_2: int,
             atom_3: int) -> np.array:
        pass

    @staticmethod
    def grad_dihedral(xyz: np.array, 
                    atom_1: int, 
                    atom_2: int,
                    atom_3: int,
                    atom_4: int) -> np.array:
        pass


    @staticmethod
    def grad_rmsd_cal(xyz: np.array,
                      xyz_ref: np.array):
        if len(xyz.shape) == 1:
            raise ValueError("The shape of xyz is wrong!")

        if len(xyz.shape) == 2:
            """xyz: (N, D)"""
            xyz = xyz[np.newaxis, :, :]
        
        """
        xyz: (1, N, D)
        xyz_ref: (I, N, D)
        """
        rmsd_val, grad_rmsd = CVs_Calculation.rmsd_main(coord1=xyz,
                                        coord2=xyz_ref)
        return grad_rmsd
    

    """CVs Calculation Control"""
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
            x_ref = np.array(self.x_ref)
            cv = CVs_Calculation.cal_rmsd(system.xyz, x_ref)
        return cv

    def grad_cv_cal(self,
                    system,
                    cvs_type: str,
                    atom_ids: list):
        # if cvs_type == 'Bond' or 'bond' or 'b' or 'B':
        if cvs_type in ['Bond', 'bond', 'b', 'B']:
            # print(cvs_type)
            grad_cv = CVs_Calculation.grad_bond_cal(system.xyz, atom_ids[0], atom_ids[1])
        # elif cvs_type == 'Angle' or 'angle' or 'a' or 'A':
        elif cvs_type in ['Angle', 'angle', 'a', 'A']:
            # print(cvs_type)
            pass
        # elif cvs_type == 'Dihedral' or 'dihedral' or 'd' or 'D':
        elif cvs_type in ['Dihedral', 'dihedral', 'd', 'D']:
            # print(cvs_type)
            pass
        # elif cvs_type == 'RMSD' or 'rmsd' or 'R':
        elif cvs_type in ['RMSD', 'rmsd', 'r', 'R']:
            x_ref = np.array(self.x_ref)             # May wrong?
            grad_cv = CVs_Calculation.grad_rmsd_cal(system.xyz, x_ref)
        return grad_cv

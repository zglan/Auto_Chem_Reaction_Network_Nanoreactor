import numpy as np
from ..tool.constant import kb_au
from scipy.special import logsumexp


class ITS:

    r"""
    Why pe_initial:
            exp(-beta * U) / Z
        where Z = \int_{-inf}^{inf} exp(-beta * U), Set U_system = U - pe_initial, then:

            exp(-beta * U_system) * C / \int_{-inf}^{inf} exp(-beta * U_system) * C
        where C = exp(-beta * pe_initial)
        
        Then:
            exp(-beta * U) / Z = C * exp(-beta * U_system) / (C * Z_system)
                               = exp(-beta * U_system) / Z_system
        More Stability For MD Simulation.

    For ITS:
            exp(-beta * U_eff) = \sum_{nk * exp(-beta_k * U)}
            exp(-beta * U_eff) / Z
        where Z = \int_{-inf}^{inf} exp(-beta * U_eff) = \int_{-inf}^{inf} \sum_{nk * exp(-beta_k * U)}, Set U_system = U - pe_initial, then:

            \sum_{nk * exp(-beta_k * U)} = \sum_{nk * exp(-beta_k * U_system) * exp(-beta_k * pe_initial)}
                                         = \sum_{nk * ck * exp(-beta_k * U_system)}
            Z = \int_{-inf}^{inf} \sum_{nk * exp(-beta_k * U)}
              = \int_{-inf}^{inf} \sum_{nk * ck * exp(-beta_k * U_system)}

        if set nk_reduced = nk * ck, or actually nk_reduced is nk, then:
            exp(-beta * U_eff) = \sum_{nk_reduced * exp(-beta_k * U_system)}
            Z = \int_{-inf}^{inf} \sum_{nk_reduced * exp(-beta_k * U_system)}
        Basically, the same.
    
    If there are many minima, how to deal with it?
        U1, U2. Where U1 > U2 and pe_initial = U1. Set U_system = U - U1, then:
            exp(-beta * U_eff) = \sum_{nk_reduced * exp(-beta_k * U_system)}
            Z = \int_{-inf}^{inf} \sum_{nk_reduced * exp(-beta_k * U_system)}
        If minima needed to be changed, set U_system_new = U - U2, then:
            exp(-beta * U_eff) = \sum_{nk_reduced_new * exp(-beta_k * U_system_new)}
            Z = \int_{-inf}^{inf} \sum_{nk_reduced_new * exp(-beta_k * U_system_new)}
        Initially, we set nk = nk_reduced = nk_real * ck, and nk_reduced_new = nk_real * ck_new, the relationship between nk_reduced
    and nk_reduced_new is given as follow:
            nk_reduced_new = nk_real * ck_new = nk_reduced * ck_new / ck
    """

    def __init__(self,
                 T_target: float,
                 temp_low: float,
                 temp_high: float,
                 bin_temp: int,
                 update_step: int = 500):
        self.kb = 1.380649e-23 / 4.3597447222071e-18   # a.u. / K
        self.T_target = T_target
        self.beta_target = 1 / (self.kb * self.T_target)

        self.T_list = np.linspace(temp_low, temp_high, num=bin_temp)
        self.kbT_list = self.T_list * self.kb
        self.beta = np.reciprocal(self.kbT_list)

        # self.nk = np.ones_like(self.beta) * 1e-2
        # self.nk = np.array([2 ** (-k-1) for k in range(len(self.beta))])
        self.nk = np.array([np.sqrt(self.T_list[0] / self.T_list[k]) for k in range(len(self.beta))])

        self.lnnk = np.log(self.nk)
        self.lnPk = np.full(len(self.beta), -np.inf)
        self.lnWk = np.full(len(self.beta)-1, -np.inf)

        self.pe_initial = 0.0      # Obtain by structure optimization or the first point.
        self.update_step = update_step
        self.bias = 1.0
        self.ln_bias = np.log(self.bias)

        self.lnnk_list = []
        self.pe_ref_list = []
        self.S_list = []

    def its_update(self):
        # print(self.lnPk)
        # print(f"Pk: {np.exp(self.lnPk)}")
        self.lnPk_norm = self.lnPk - logsumexp(self.lnPk)
        # exit()
        # print(f"lnPk_norm: {self.lnPk_norm }")
        # print(np.exp(self.lnPk_norm))
        # exit()
        lnmk = self.lnnk[:-1] - self.lnnk[1:]
        # print(f"lnmk: {lnmk}")
        # print(np.exp(lnmk))
        # exit()
        ln_w = 0.5 * (self.lnPk_norm[1:] + self.lnPk_norm[:-1])
        # print(np.exp(ln_w))
        # exit()
        ln_Pk_norm_ratio = self.lnPk_norm[1:] - self.lnPk_norm[:-1]
        # print(np.exp(ln_Pk_norm_ratio))
        # exit()
        lnw_norm = np.logaddexp(self.lnWk, ln_w)
        # print(np.exp(lnw_norm))
        # exit()
        lnwbias = self.ln_bias + ln_w
        # print(np.exp(lnwbias))
        # exit()
        lnwbias_norm = np.logaddexp(self.lnWk, lnwbias)
        # print(np.exp(lnwbias_norm))
        # exit()
        lnmk_new = lnmk + np.logaddexp(lnwbias + ln_Pk_norm_ratio, self.lnWk) - lnwbias_norm
        # print(np.exp(lnmk_new))
        # exit()
        lnnk_new = [0.]
        for i in range(len(self.lnnk)-1):
            lnnk_temp = lnnk_new[i] - lnmk_new[i]
            # lnnk_temp = lnnk_new[i] - lnmk_new[i] + np.log(self.T_list[i] / self.T_target)    # Test
            # print(np.exp(lnnk_temp))
            # exit()
            lnnk_new.append(lnnk_temp)
        lnnk_new = np.array(lnnk_new)
        # print(f"lnnk new: ",lnnk_new)

        self.lnWk = lnw_norm
        self.lnnk = lnnk_new
        self.lnPk = np.full_like(self.beta, -np.inf)

        self.lnnk_list.append(self.lnnk)
        

    def its_force(self, U, pe_initial, lnPk):
        """U: potential energy"""
        U = U - pe_initial
        upper = logsumexp(self.lnnk + np.log(self.beta) - self.beta * U)
        # print(f"upper: {np.sum(np.exp(self.lnnk) * self.beta * np.exp(- self.beta * U))}")
        # print(f"upper: {np.exp(upper)}")
        # exit()
        lower = logsumexp(self.lnnk - self.beta * U)
        # print(f"lower: {np.sum(np.exp(self.lnnk) * np.exp(- self.beta * U))}")
        # print(f"lower: {np.exp(lower)}")
        # exit()
        S = np.exp(upper - lower) / self.beta_target
        print(f"enhanced factor: {S}")
        
        lnPk = np.logaddexp(lnPk, self.lnnk - self.beta * U)
        # print(np.exp(lnPk))
        # print(f"lnPk: {lnPk}")
        # Pk = np.zeros(len(self.beta))
        # Pk = Pk + self.nk * np.exp(-self.beta * U)
        # print(Pk)
        # print(self.nk * np.exp(-self.beta * U))
        return S, lnPk

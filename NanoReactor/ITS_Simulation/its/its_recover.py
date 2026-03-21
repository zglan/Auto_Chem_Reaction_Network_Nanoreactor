import numpy as np

class ITS_recover:
    
    @ staticmethod
    def its_only(nk, beta_k, beta_target, U):
        """
        ITS only
        
        Probability Density of U: rou = exp(-beta * U) / Q, where Q is partition function.
        According to ITS:
            rou = exp(-beta * U) / Q
                = c0 * exp(-beta * U_its) / Q
        So:
            c0 = exp(-beta * U) * exp(beta * U_its)
               = sum_{k}{nk * exp(beta_k * U)} * exp(-beta * U)
               = 1 / sum_{k}{nk * exp((beta - beta_k) * U)}
        """
        weight = np.sum(nk * np.exp((beta_target - beta_k) * U))
        c0 = 1 / weight
        return c0

    @ staticmethod
    def its_bias(nk, beta_k, beta_target, U, U_mtd, U_nano):
        """
        ITS with additional potential like V_mtd, V_wall and so on.
        
        Probability Density of U: rou = exp(-beta * U) / Q, where Q is partition function.
        According to ITS:
            rou = exp(-beta * U) / Q
                = c0 * exp(-beta * U_its) / Q
        Because:
            U_system = U + U_mtd + U_nano
            U_additon = U_mtd + U_nano 
            U_its = -ln{sum_{k}{nk * exp(-beta_k * U_system)}} / beta
        Then:
            c0 = exp(-beta * U) * exp(beta * U_its)
               = sum_{k}{nk * exp(beta_k * U_system)} * exp(-beta * U)
               = exp(-beta * U_addition) / sum_{k}{nk * exp((beta - beta_k) * U_system)}
        """
        U_system = U + U_mtd + U_nano
        U_addition = U_mtd + U_nano
        weight = np.sum(nk * np.exp((beta_target - beta_k) * U_system))
        c0 = np.exp(beta_target * U_addition) / weight
        return c0
    

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path

class ExactRiemann:
    """ Class to solve the exact 1D Riemann problem """
    def __init__(self, rho_L, u_L, p_L, rho_R, u_R, p_R, gamma, casename):
        # Casename for plotting and printing
        self.casename = casename
        # Input parameters
        self.rho_L = rho_L
        self.u_L = u_L
        self.p_L = p_L
        self.rho_R = rho_R
        self.u_R = u_R
        self.p_R = p_R
        self.gamma = gamma
        self.rtol = 1e-6
        self.a_L = np.sqrt(self.gamma * self.p_L / self.rho_L)
        self.a_R = np.sqrt(self.gamma * self.p_R / self.rho_R)

        # Solution in the star region
        self.p_star = 0.0
        self.rho_star = 0.0
        self.rhoL_star = 0.0
        self.rhoR_star = 0.0
    
    def solve_pstar(self, initial_guess='mean'):
        """ Find the pressure in the star region """
        maxits = 10
        # Pressure at step k and k + 1
        if initial_guess == 'mean':
            p_k = 0.5 * (self.p_L + self.p_R)
        elif initial_guess == 'TR':
            p_k = ((self.a_L + self.a_R - 0.5 * (self.gamma - 1) * (self.u_R - self.u_L)) 
                / (self.a_L / self.p_L**((self.gamma - 1) / 2 / self.gamma) 
                    + self.a_L / self.p_L**((self.gamma - 1) / 2 / self.gamma)))**(2 * self.gamma / (self.gamma - 1))
        it = 0
        err = np.inf
        while (err > self.rtol):
            p_km1 = p_k
            p_k -= pressure_function(p_k, self.rho_L, self.u_L, self.p_L, 
                    self.rho_R, self.u_R, self.p_R, self.gamma) \
                / dpressure_function(p_k, self.rho_L, self.u_L, self.p_L, 
                    self.rho_R, self.u_R, self.p_R, self.gamma)
            err = 2 * abs(p_k - p_km1) / (p_k + p_km1)
            if it > maxits:
                print('Maximum number of iterations reached')
                return p_k
        return p_k
    
    def solve_riemann(self, initial_guess='mean'):
        """ Solve the Riemann problem by providing p_star, u_star, rhoL_star, rhoR_star """
        self.p_star = self.solve_pstar(initial_guess)
        self.u_star = 0.5 * (self.u_L + self.u_R 
            + f_K(self.p_star, self.rho_R, self.p_R, self.gamma) 
            - f_K(self.p_star, self.rho_L, self.p_L, self.gamma))
        self.rhoL_star = rho_star(self.rho_L, self.p_star / self.p_L, self.gamma)
        self.rhoR_star = rho_star(self.rho_R, self.p_star / self.p_R, self.gamma)
    
    def __str__(self):
        result_str = f'{self.casename}:\n'
        result_str += f'{self.rho_L:10.3e} | {self.u_L:10.3e} | {self.p_L:10.3e}\n'
        result_str += f'{self.p_star:10.3e} | {self.u_star:10.3e} | {self.rhoL_star:10.3e} | {self.rhoR_star:10.3e}\n'
        result_str += f'{self.rho_R:10.3e} | {self.u_R:10.3e} | {self.p_R:10.3e}\n'
        return result_str

def f_K(p, rho_K, p_K, gamma):
    """ Return the implicit function to find the pressure (Toro Chaper 4) """
    A_K = 2 / (gamma + 1) / rho_K
    B_K = (gamma - 1) / (gamma + 1) * p_K
    ss_K = np.sqrt(gamma * p_K / rho_K)
    return np.where(p > p_K, (p - p_K) * np.sqrt(A_K / (p + B_K)),
                2 * ss_K / (gamma - 1) * ((p / p_K)**((gamma - 1) / 2 / gamma) - 1))

def df_K(p, rho_K, p_K, gamma):
    """ Return the derivative of f_K """
    A_K = 2 / (gamma + 1) / rho_K
    B_K = (gamma - 1) / (gamma + 1) * p_K
    ss_K = np.sqrt(gamma * p_K / rho_K)
    return np.where(p > p_K, np.sqrt(A_K / (B_K + p)) * (1 - (p - p_K) / 2 / (B_K + p)),
                (p / p_K)**(-(gamma + 1) / 2 / gamma) / rho_K / ss_K)

def pressure_function(p, rho_L, u_L, p_L, rho_R, u_R, p_R, gamma):
    """ Return the implicit pressure function """
    return f_K(p, rho_L, p_L, gamma) + f_K(p, rho_R, p_R, gamma) + u_R - u_L

def dpressure_function(p, rho_L, u_L, p_L, rho_R, u_R, p_R, gamma):
    """ Return the implicit pressure function """
    return df_K(p, rho_L, p_L, gamma) + df_K(p, rho_R, p_R, gamma)

def rho_rarefaction(rho_K, press_ratio, gamma):
    """ Rho star as a function of 
    shock parameters from isentropic law """
    return rho_K * press_ratio**(1 / gamma)

def rho_shock(rho_K, press_ratio, gamma):
    """ Rho star as a function of
    shock parameters from Rankine-Hugoniot relations """
    frac_g = (gamma - 1) / (gamma + 1)
    return rho_K * (press_ratio + frac_g) / (frac_g * press_ratio + 1)

def rho_star(rho_K, press_ratio, gamma):
    """ Depending on the value of the pressure 
    inside the star region the density is computed """
    return np.where(press_ratio > 1, rho_shock(rho_K, press_ratio, gamma), 
        rho_rarefaction(rho_K, press_ratio, gamma))

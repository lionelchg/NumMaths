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

        # Gamma related variables
        self.gp12 = (self.gamma + 1) / 2
        self.gm12 = (self.gamma - 1) / 2
        self.gp1ogm1 = (self.gamma + 1) / (self.gamma - 1)
        self.gp12og = self.gp12 / self.gamma
        self.gm12og = self.gm12 / self.gamma

        # Solution in the star region
        self.p_star = 0.0
        self.rho_star = 0.0
        self.rhoL_star = 0.0
        self.rhoR_star = 0.0
        self.SL = None
        self.SR = None
    
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
        elif isinstance(initial_guess, float):
            p_k = initial_guess
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
    
    def plot_pfunction(self, figname):
        """ Plot the pressure function for debugging """
        p = np.linspace(0.8 * self.p_L, 1.2 * self.p_R, 201)
        fig, ax = plt.subplots()
        ax.plot(p, pressure_function(p, self.rho_L, self.u_L, self.p_L, 
            self.rho_R, self.u_R, self.p_R, self.gamma))
        ax.plot(p, np.zeros_like(p), color='k')
        ax.grid(True)
        ax.set_xlabel('$p$')
        ax.set_ylabel('$f(p)$')
        ax.legend()
        ax.spines['top'].set_color('none')
        ax.spines['right'].set_color('none')
        fig.savefig(figname, bbox_inches='tight')
    
    def solve_riemann(self, initial_guess='mean'):
        """ Solve the Riemann problem by providing p_star, u_star, rhoL_star, rhoR_star """
        self.p_star = self.solve_pstar(initial_guess)
        self.u_star = 0.5 * (self.u_L + self.u_R 
            + f_K(self.p_star, self.rho_R, self.p_R, self.gamma) 
            - f_K(self.p_star, self.rho_L, self.p_L, self.gamma))
        self.rhoL_star = rho_star(self.rho_L, self.p_star / self.p_L, self.gamma)
        self.rhoR_star = rho_star(self.rho_R, self.p_star / self.p_R, self.gamma)
        self.aL_star = np.sqrt(self.gamma * self.p_star / self.rhoL_star)
        self.aR_star = np.sqrt(self.gamma * self.p_star / self.rhoR_star)

        # Compute shock speeds or head and tail speeds if necessary
        if self.p_star > self.p_L:
            self.SL = self.u_L - self.a_L * np.sqrt(self.gp12og * self.p_star / self.p_L 
                        + self.gm12og)
        else:
            # Rarefaction wave, compute the speeds
            self.aL_star = self.a_L * (self.p_star / self.p_L)**((self.gamma - 1) / 2 / self.gamma)
            self.S_HL = self.u_L - self.a_L
            self.S_TL = self.u_star - self.aL_star

        if self.p_star > self.p_R:
            self.SR = self.u_R + self.a_R * np.sqrt(self.gp12og * self.p_star / self.p_R 
                        + self.gm12og)
        else:
            # Rarefaction wave, compute the speeds
            self.aR_star = self.a_R * (self.p_star / self.p_R)**((self.gamma - 1) / 2 / self.gamma)
            self.S_HR = self.u_star + self.aR_star
            self.S_TR = self.u_R + self.a_R
            
    def __str__(self):
        result_str = f'{self.casename}:\n'
        if self.p_star > self.p_L:
            result_str += f'Left Shock:  S_L = {self.SL:.3e} m/s\n'
        else:
            result_str += f'Left RW:     S_HL = {self.S_HL:.3e} m/s - S_TL = {self.S_TL:.3e} m/s\n'
        if self.p_star > self.p_R:
            result_str += f'Right Shock: S_R = {self.SR:.3e} m/s\n'
        else:
            result_str += f'Right RW:    S_TR = {self.S_TR:.3e} m/s - S_HR = {self.S_HR:.3e} m/s\n'
        result_str += f"         | {'rho':^10s} | {'u':^10s} | {'p':^10s} | {'a':^10s}\n"
        result_str += '-'*60 + '\n'
        result_str += f'W_L:     | {self.rho_L:10.3e} | {self.u_L:10.3e} | {self.p_L:10.3e} | {self.a_L:10.3e}\n'
        result_str += f'W_Lstar: | {self.rhoL_star:10.3e} | {self.u_star:10.3e} | {self.p_star:10.3e} | {self.aL_star:10.3e}\n'
        result_str += f'W_Rstar: | {self.rhoR_star:10.3e} | {self.u_star:10.3e} | {self.p_star:10.3e} | {self.aR_star:10.3e}\n'
        result_str += f'W_R:     | {self.rho_R:10.3e} | {self.u_R:10.3e} | {self.p_R:10.3e} | {self.a_R:10.3e}\n'

        return result_str
    
    def construct_sol(self, x, time):
        """ Construct the solution of Riemann problem at time t and for vector x """
        x_left = x[x / time <= self.u_star]
        left_sol = np.zeros((len(x_left), 3))
        if self.p_star > self.p_L:
            # Shock
            left_sol[:, 0] = np.where(x_left / time < self.SL, self.rho_L, self.rhoL_star)
            left_sol[:, 1] = np.where(x_left / time < self.SL, self.u_L, self.u_star)
            left_sol[:, 2] = np.where(x_left / time < self.SL, self.p_L, self.p_star)
        else:
            # Filter the vectors into three parts
            x_L = x_left[x_left / time < self.S_HL]
            x_rwave = x_left[(x_left / time >= self.S_HL) & (x_left / time <= self.S_TL)]
            x_Lstar = x_left[x_left / time > self.S_TL]

            # Compute the fan
            W_Lfan = np.zeros((len(x_rwave), 3))
            # Density
            W_Lfan[:, 0] = self.rho_L * (2 / (self.gamma + 1) 
                + (self.gamma - 1) / (self.gamma + 1) / self.a_L 
                * (self.u_L - x_rwave / time))**(2 / (self.gamma - 1))
            # Speed
            W_Lfan[:, 1] = 2 / (self.gamma + 1) * (self.a_L 
                + (self.gamma - 1) / 2 * self.u_L + x_rwave / time)
            # Pressure
            W_Lfan[:, 2] = self.p_L * (2 / (self.gamma + 1) 
                + (self.gamma - 1) / (self.gamma + 1) / self.a_L 
                * (self.u_L - x_rwave / time))**(2 * self.gamma / (self.gamma - 1))

            # Reconstruct solution
            left_sol[:, 0] = np.concatenate((self.rho_L * np.ones_like(x_L), W_Lfan[:, 0], self.rhoL_star * np.ones_like(x_Lstar)))
            left_sol[:, 1] = np.concatenate((self.u_L * np.ones_like(x_L), W_Lfan[:, 1], self.u_star * np.ones_like(x_Lstar)))
            left_sol[:, 2] = np.concatenate((self.p_L * np.ones_like(x_L), W_Lfan[:, 2], self.p_star * np.ones_like(x_Lstar)))
        
        x_right = x[x / time > self.u_star]
        right_sol = np.zeros((len(x_right), 3))
        if self.p_star > self.p_R:
            # Shock
            right_sol[:, 0] = np.where(x_right / time > self.SR, self.rho_R, self.rhoR_star)
            right_sol[:, 1] = np.where(x_right / time > self.SR, self.u_R, self.u_star)
            right_sol[:, 2] = np.where(x_right / time > self.SR, self.p_R, self.p_star)
        else:
            # Filter the vectors into three parts
            x_R = x_right[x_right / time < self.S_HR]
            x_rwave = x_right[(x_right / time >= self.S_HR) & (x_right / time <= self.S_TR)]
            x_Rstar = x_right[x_right / time > self.S_TR]

            # Compute the fan
            W_Rfan = np.zeros((len(x_rwave), 3))
            # Density
            W_Rfan[:, 0] = self.rho_R * (2 / (self.gamma + 1) 
                - (self.gamma - 1) / (self.gamma + 1) / self.a_R 
                * (self.u_R - x_rwave / time))**(2 / (self.gamma - 1))
            # Speed
            W_Rfan[:, 1] = 2 / (self.gamma + 1) * (- self.a_R 
                + (self.gamma - 1) / 2 * self.u_R + x_rwave / time)
            # Pressure
            W_Rfan[:, 2] = self.p_R * (2 / (self.gamma + 1) 
                - (self.gamma - 1) / (self.gamma + 1) / self.a_R 
                * (self.u_R - x_rwave / time))**(2 * self.gamma / (self.gamma - 1))

            # Reconstruct solution
            right_sol[:, 0] = np.concatenate((self.rhoR_star * np.ones_like(x_R), W_Rfan[:, 0], self.rho_R * np.ones_like(x_Rstar)))
            right_sol[:, 1] = np.concatenate((self.u_star * np.ones_like(x_R), W_Rfan[:, 1], self.u_R * np.ones_like(x_Rstar)))
            right_sol[:, 2] = np.concatenate((self.p_star * np.ones_like(x_R), W_Rfan[:, 2], self.p_R * np.ones_like(x_Rstar)))
        
        sol = np.concatenate((left_sol, right_sol), axis=0)
        return sol
    
    
    def plot_solution(self, x, x0, t, figname):
        """ Plot the solution for vector x at instant t """
        fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(8, 8))
        sol = self.construct_sol(x - x0, t)
        internal_energy = sol[:, 2] / sol[:, 0] / (self.gamma - 1)
        axes[0][0].plot(x, sol[:, 0], 'k')
        ax_prop(axes[0][0], r'$x$ [m]', r'$\rho$ [kg/m$^3$]')
        axes[0][1].plot(x, sol[:, 1], 'k')
        ax_prop(axes[0][1], r'$x$ [m]', r'$u$ [m/s]')
        axes[1][0].plot(x, sol[:, 2], 'k')
        ax_prop(axes[1][0], r'$x$ [m]', r'$p$ [Pa]')
        axes[1][1].plot(x, internal_energy, 'k')
        ax_prop(axes[1][1], r'$x$ [m]', r'$e$ [J]')
        fig.tight_layout()
        fig.savefig(figname, bbox_inches='tight')
    
    def print_solution(self, x, x0, t, filename):
        """ Print the solution to filename """
        sol = self.construct_sol(x - x0, t)
        fp = open(filename, 'w')
        fp.write(f"{'x':^10s} {'rho':^10s} {'u':^10s} {'p':^10s}")
        for i in range(len(x)):
            fp.write(f'{x[i]:10.3e} {sol[i, 0]:10.3e} {sol[i, 1]:10.3e} {sol[i, 2]:10.3e}\n')
        fp.close()

def ax_prop(ax, xlabel, ylabel):
    ax.grid(True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

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

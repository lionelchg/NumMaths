import numpy as np
import yaml
import argparse
import scipy.constants as co
import matplotlib.pyplot as plt

from pathlib import Path

class Grid:
    def __init__(self, cfg):
        
        # Initialize grid:
        self.npts = cfg['grid']['npts']
        self.x_min, self.x_max = cfg['grid']['x_min'], cfg['grid']['x_max']
        self.x = np.linspace(self.x_min, self.x_max, self.npts, endpoint=True)
        self.dx = self.compute_dx(self.npts, self.x)
        
        # Initialize mixture / species
        self.spec_names = cfg['mixture']['spec_names']
        self.sp_ind = {sp: i for i, sp in enumerate(self.spec_names)}
        self.nsp = len(self.spec_names)
        
        # Species particle density, its residuals and its gradient
        self.w = np.zeros((self.nsp, self.npts))
        self.dw = np.zeros((self.nsp, self.npts))
        self.gradw = np.zeros((self.nsp, self.npts))
        # Species fluxes
        self.flux_w = np.zeros((self.nsp, self.npts))
        # Species mobility coefficient
        self.mu = np.zeros((self.nsp, self.npts))
        # Species diffusion coefficient
        self.D = np.zeros((self.nsp, self.npts))
        # Species charge
        self.charge = self.compute_charge()
        
        # Physical properties:
        self.T = cfg['gas']['T']
        self.P = cfg['gas']['P']
        self.N = self.P / ( self.T * co.k )
        
        # initialize species profiles
        if 'init' in cfg:
            self.init_species(cfg)
        
        # Initialize some variables:
        # Electric field
        self.E = np.zeros(self.npts)
        
    @staticmethod
    def compute_dx(npts, x):
        dx = np.zeros(npts)
        dx[0] = 0.5 * (x[1] - x[0])
        dx[-1] = 0.5 * (x[-1] - x[-2])
        for n in range(1, npts - 1):
            xi = 0.5 * (x[n - 1] + x[n])
            xj = 0.5 * (x[n] + x[n + 1])
            dx[n] = xj - xi
            
        return dx
        
    def init_species(self, cfg):
        """ Initialize species with gaussian form """
        for i in range(self.nsp):
            if self.spec_names[i] in cfg['init']:
                self.init_single_species(self.spec_names[i], **cfg['init'][self.spec_names[i]])
        
    def init_single_species(self, name, x0, A, sigma, A0 = 0.0):
        """ Initialize a single species based on its name """
        self.w[self.sp_ind[name], :] = self.gaussian(self.x, x0, A, sigma, A0) 
    
    @staticmethod
    def gaussian(x, x0, A, sigma_x, A0 = 0.0):
        """ Gaussian with amplitude A and variance sigma_x with background A0 """
        return A0 + A * np.exp(- (x - x0)**2 / sigma_x**2)
    
    def compute_charge(self):
        """ Compute the charge of species (+/- 1) """
        charge = np.zeros(self.nsp)
        for n in range(self.nsp):
            if self.spec_names[n] == 'Electron':
                charge[n] = - 1.0
            elif self.spec_names[n][-1] == '+':
                charge[n] = 1.0
            elif self.spec_names[n][-1] == '-':
                charge[n] = - 1.0
        return charge
    
    def compute_gradient(self):
        """ Compute gradient of w """
        for s in range(self.nsp):
            self.gradw[s, 0] = (self.w[s, 1] - self.w[s, 0]) / (self.x[1] - self.x[0])
            self.gradw[s, -1] = (self.w[s, -1] - self.w[s, -2]) / (self.x[-1] - self.x[-2])
            self.gradw[s, 1:-1] = (self.w[s, 2:] - self.w[s, 0:-2]) / (self.x[2:] - self.x[0:-2])
        
    def compute_ambipolar(self):
        """ Compute the ambipolar electric field """
        
        # Num / Denom for ambipolar formula
        A, B = np.zeros(self.npts), np.zeros(self.npts)
        
        for n in range(self.npts):
            for i in range(self.nsp):
                A[n] += self.charge[i] * self.D[i, n] * self.gradw[i, n]
                B[n] += self.mu[i, n] * self.w[i, n]
                
        self.E = A / B
        
    def compute_thermo(self):
        """ Compute mobility and diffusion coefficients """
        
        for n in range(self.nsp):
            if self.spec_names[n] == "Electron":
                self.mu[n, :] = self.electron_mobility(self.npts, self.E, self.N)
                self.D[n, :] = self.electron_diffusion(self.npts, self.E, self.mu[n, :], self.N)
            elif self.spec_names[n][-1] == "+":
                self.mu[n, :] = self.positive_ion_mobility(self.N)
                self.D[n, :] = self.einstein_mu2D(self.mu[n, :], self.T)
            elif self.spec_names[n][-1] == "-":
                self.mu[n, :] = self.negative_ion_mobility(self.npts, self.E, self.N)
                self.D[n, :] = self.einstein_mu2D(self.mu[n, :], self.T)
            
    @staticmethod
    def electron_mobility(npts, E, N):
        """ Electron mobility from Morrow and Lowke """
        
        mu = np.zeros(npts)
        for n in range(npts):
            Eprim = np.abs(E[n]) * 1.0e-2   # Units [V.cm-1]
            E_N = np.abs(E[n]) / N * 1.0e4  # Units [V.cm2]
            if ( E_N > 2.0e-15 ):
                mu[n] = (7.4e21 * E_N + 7.1e6) / Eprim
            elif ( E_N > 1e-16 and E_N <= 2.0e-15 ):
                mu[n] = (1.03e22 * E_N + 1.3e6) / Eprim
            elif ( E_N > 2.6e-17 and E_N <= 1e-16 ):
                mu[n] = (7.2973e21 * E_N + 1.63e6) / Eprim
            elif ( E_N > 1e-19 ):
                mu[n] = (6.87e22 * E_N + 3.38e4) / Eprim
            else: # E_N < 0.01 Td As in Tholin page 12
                mu[n] = (6.87e22 * 1e-19 + 3.38e4) / ( (1e-19 * N / 1e4) * 1e-2 ) # Constant mobility at low electric field.
            
            # [m2/V/s]    [cm2/V/s]
            mu[n] = mu[n] * 1.0e-4 # m2.V-1.s-1
        
        return mu
    
    @staticmethod
    def electron_diffusion(npts, E, mu, N):
        """ Electron diffusion from Morrow and Lowke """
        
        D = np.zeros(npts)
        for n in range(npts):
            E_N = np.abs(E[n]) / N * 1e4

            if ( E_N > 1e-19 ): # > 0.01 Td
                D[n] = (0.3341e9 * E_N**0.54069) * (mu[n] * 1e4) * 1e-4   # Mobility in cm2.V-1.s-1
            else: # Same limit as in Tholin PhD page 13 
                D[n] = (0.3341e9 * 1e-19**0.54069) * (mu[n] * 1e4) * 1e-4 # Mobility is computed before diffusion coeff 
                                                                          # Criterium for low E field already applied 
        return D
    
    @staticmethod
    def positive_ion_mobility(N):
        """ Positive ion mobility from Morrow and Lowke """
        return 2.34e-4 * (N / 2.45e25)

    @staticmethod
    def negative_ion_mobility(npts, E, N):
        """ Negative ion mobility from Morrow and Lowke """
        mu = np.zeros(npts)
        for n in range(npts):
            E_N = np.abs(E[n]) / N * 1e4
            if ( E_N > 5e-16 ):
                mu[n] = 2.7e-4 * (N/ 2.45e25)
            else:
                mu[n] = 1.86e-4 * (N / 2.45e25)
        return mu
    
    @staticmethod
    def einstein_mu2D(mu, T):
        """ Compute diffusion coefficient using mobility coefficient according to Einstein relation """
        return mu * T * co.k / co.e


def main():
    """ Wrapper around parser """
    args = argparse.ArgumentParser(description='Single case run')
    args.add_argument('-c', '--config', type=str,
                        help='Config filename', required=True)
    args = args.parse_args()

    with open(args.config, 'r') as yaml_stream:
        cfg = yaml.safe_load(yaml_stream)
        
    # Figures directory
    fig_dir = Path('figures')
    fig_dir.mkdir(parents=True, exist_ok=True)
    
    grid = Grid(cfg)
    
    
    

if __name__ == '__main__':
    main()


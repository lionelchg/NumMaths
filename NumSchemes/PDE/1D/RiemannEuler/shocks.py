import numpy as np
import scipy.constants as co

def density_ratio(rel_mach, gamma):
    """ Density ratio as a function of relative mach number """
    return (gamma + 1) * rel_mach**2 / ((gamma - 1) * rel_mach**2 + 2)

def pressure_ratio(rel_mach, gamma):
    """ Pressure ratio as a function of relative mach number """
    return (2 * gamma * rel_mach**2 - (gamma - 1)) / (gamma + 1)

def rel_mach(press_ratio, gamma):
    """ Relative mach number for a 1-wave (+ sqrt) 3-wave (- sqrt) shock """
    return np.sqrt((gamma + 1) / 2 / gamma * press_ratio + (gamma - 1) / 2 / gamma)

if __name__ == '__main__':
    # Density, speed and pressure of right state
    rhoR = 1.225
    uR = 0.0
    pR = 1.01325e5
    gamma = 1.4

    # Mach number
    M3 = 3

    # Compute right sound speed and mach numbers
    aR = np.sqrt(gamma * pR / rhoR)
    S3 = M3 * aR
    MR = uR / aR
    M_RS = MR - M3

    # Deduce density and pressure
    rho_star = density_ratio(M_RS, gamma) * rhoR
    p_star = pressure_ratio(M_RS, gamma) * pR
    u_star = S3 * (1 - rhoR / rho_star) + rhoR / rho_star * uR
    a_star = np.sqrt(gamma * p_star / rho_star)

    # Print quantities
    print(f'Input data:')
    print(f'rhoR = {rhoR:.3e} kg/m3 - uR = {uR:.3e} m/s - pR = {pR:.3e} Pa - M3 = {M3:.1f}')
    print(f'aR = {aR:.2e} m/s - S3 = {S3:.2e} m/s - M_R = {MR:.1f} - (M_R - M_S) = {M_RS:.1f}')
    print(f'Star state:')
    print(f'rho_star = {rho_star:.3e} kg/m3 - u_star = {u_star:.3e} m/s - p_star = {p_star:.3e} Pa - a_star = {a_star:.3e} m/s')

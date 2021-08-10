import numpy as np
import scipy.constants as co
import matplotlib.pyplot as plt
from pathlib import Path

def density_ratio(rel_mach, gamma):
    """ Density ratio as a function of relative mach number """
    return (gamma + 1) * rel_mach**2 / ((gamma - 1) * rel_mach**2 + 2)

def pressure_ratio(rel_mach, gamma):
    """ Pressure ratio as a function of relative mach number """
    return (2 * gamma * rel_mach**2 - (gamma - 1)) / (gamma + 1)

def inside_speed(S3, rhoR, rho_star, uR):
    """ Return the star region speed """
    return S3 * (1 - rhoR / rho_star) + rhoR / rho_star * uR

def rel_mach(press_ratio, gamma):
    """ Relative mach number for a 1-wave (+ sqrt) 3-wave (- sqrt) shock """
    return np.sqrt((gamma + 1) / 2 / gamma * press_ratio + (gamma - 1) / 2 / gamma)

def star_from_right(rhoR, uR, pR, MS, gamma, lprint=False):
    """ Properties inside star region for given right state and shock mach number """
    # Compute right sound speed and mach numbers
    aR = np.sqrt(gamma * pR / rhoR)
    S3 = MS * aR
    MR = uR / aR
    M_RS = MR - MS

    # Deduce density and pressure
    rho_star = density_ratio(M_RS, gamma) * rhoR
    p_star = pressure_ratio(M_RS, gamma) * pR
    u_star = inside_speed(S3, rhoR, rho_star, uR)
    a_star = np.sqrt(gamma * p_star / rho_star)

    if lprint:
        # Print quantities
        print(f'Input data:')
        print(f'rhoR = {rhoR:.3e} kg/m3 - uR = {uR:.3e} m/s - pR = {pR:.3e} Pa - MS = {MS:.1f}')
        print(f'aR = {aR:.2e} m/s - S3 = {S3:.2e} m/s - M_R = {MR:.1f} - (M_R - M_S) = {M_RS:.1f}')
        print(f'Star state:')
        print(f'rho_star = {rho_star:.3e} kg/m3 - u_star = {u_star:.3e} m/s - p_star = {p_star:.3e} Pa - a_star = {a_star:.3e} m/s')

    return rho_star, p_star, u_star, a_star

def ax_prop(ax, title):
    ax.grid(True)
    ax.set_title(title)

if __name__ == '__main__':
    # Density, speed and pressure of right state
    rhoR = 1.225
    uR = 0.0
    pR = 1.01325e5
    gamma = 1.4

    # Compute star region physical values
    star_from_right(rhoR, uR, pR, 3.0, gamma, lprint=True)

    # Plotting tendencies for fixed right state and varying shock mach number
    MS = np.linspace(1, 5, 201)
    rho_star, p_star, u_star, a_star = star_from_right(rhoR, uR, pR, MS, gamma)

    fig_dir = Path('figures')
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(8, 8), sharex=True)
    # rho_star
    axes[0][0].plot(MS, rho_star)
    ax_prop(axes[0][0], r'$\rho_*$ [kg/m$^3$]')
    # p_star
    axes[0][1].plot(MS, p_star / 1e5)
    ax_prop(axes[0][1], r'$p_*$ [bar]')
    # u_star
    axes[1][0].plot(MS, u_star)
    ax_prop(axes[1][0], r'$u_*$ [m/s]')
    axes[1][0].set_xlabel(r'$M_S$')
    # a_star
    axes[1][1].plot(MS, a_star)
    ax_prop(axes[1][1], r'$a_*$ [m/s]')
    axes[1][1].set_xlabel(r'$M_S$')
    fig.savefig(fig_dir / 'shock_params', bbox_inches='tight')

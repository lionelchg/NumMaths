import os
import numpy as np
import cmath
import matplotlib.pyplot as plt
from pathlib import Path

def ampl_factor(phi, cfl, fourier, scheme):
    """ Return the amplfication factor for advection diffusion scheme """
    if scheme == 'FOU_CD':
        return 1 - cfl * (1 - np.exp(-1j * phi)) + 2 * fourier * (np.cos(phi) - 1)
        # return 1 - cfl * (1 - np.exp(-1j * phi)) + 4 * fourier * np.sin(phi)**2
    elif scheme == 'CD_CD':
        return 1 - 1j * cfl * np.sin(phi) + 2 * fourier * (np.cos(phi) - 1)
        # return 1 - 1j * cfl * np.sin(phi) + 4 * fourier * np.sin(phi)**2

def plot_G(scheme, cfls, fouriers, fig_dir):
    """ Plot the diffusion and dispersion errors from the amplification factor """
    phi = np.linspace(0, np.pi, 300)
    phi_deg = phi * 180 / np.pi    
    fig, axes = plt.subplots(ncols=3, figsize=(10, 4))
    for (cfl, fourier) in zip(cfls, fouriers):
        df_err, dp_err = errors(ampl_factor(phi, cfl, fourier, scheme), phi, cfl, fourier)
        axes[0].plot(phi_deg, df_err, label=f'CFL = {cfl:.2f} - Fourier = {fourier:.2f}')
        axes[1].plot(phi_deg, dp_err, label=f'CFL = {cfl:.2f} - Fourier = {fourier:.2f}')
        axes[2].plot(phi_deg, dp_err * phi_deg, label=f'CFL = {cfl:.2f} - Fourier = {fourier:.2f}')
    ax_prop_G(axes[0], r'$\varepsilon_D$')
    ax_prop_G(axes[1], r'$\varepsilon_\phi$')
    ax_prop_G(axes[2], r'$\phi_\mathrm{num}$')
    axes[2].legend()
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    fig.suptitle(f'{scheme} Spectral Analysis')
    fig.savefig(fig_dir / f'errors_{scheme}', bbox_inches='tight')

def ax_prop_G(ax, ylabel, ylim=None):
    """ Ax properties of plot_G """
    ax.grid(True)
    ax.set_xlabel(r'$\phi$')
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=0)
    ax.set_xlim([0, 180])

def errors(G_num, phi, sigma, beta):
    """ Computation of diffusion and dispersion error for a constant advection
    speed problem """
    diff_err = abs(G_num) / np.exp(- beta * phi**2)
    disp_err = np.zeros_like(phi)
    disp_err[0] = 1
    disp_err[1:] = np.array([- cmath.phase(G_num[i]) / sigma / phi[i] for i in range(1, len(phi))])
    return diff_err, disp_err

if __name__ == '__main__':
    fig_dir = Path('figures')
    fig_dir.mkdir(parents=True, exist_ok=True)
    plot_G("FOU_CD", [0.1, 0.1, 0.1], [1e-2, 5e-2, 0.2], fig_dir)
    plot_G("CD_CD", [0.1, 0.1, 0.1], [1e-2, 5e-2, 0.2], fig_dir)
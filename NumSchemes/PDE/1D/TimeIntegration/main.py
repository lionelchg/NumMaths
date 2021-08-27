import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

def ax_prop(ax, xlabel, ylabel, title):
    ax.grid(True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()

def fourier_symbol(phi, sigma, scheme):
    """ Fourier symbol of the specified scheme at cfl *sigma* """
    if scheme == 'FOU':
        return - sigma * (1 - np.exp(-1j * phi))
    elif scheme == 'SOU':
        return - sigma / 2 * (3 - 4 * np.exp(-1j * phi) + np.exp(-1j * 2 * phi))
    elif scheme == 'Fromm':
        return - sigma / 4 * (np.exp(1j * phi) + 3 - 5 * np.exp(-1j * phi) + np.exp(-1j * 2 * phi))
    elif scheme == 'Quick':
        return - sigma / 8 * (3 * np.exp(1j * phi) + 3 - 7 * np.exp(-1j * phi) + np.exp(-1j * 2 * phi))
    elif scheme == 'Third':
        return - sigma / 6 * (2 * np.exp(1j * phi) + 3 - 6 * np.exp(-1j * phi) + np.exp(-1j * 2 * phi))

if __name__ == '__main__':
    # Default params of matplotlib overwritten
    mpl.rc('lines', linewidth=2)

    # Figures directory
    fig_dir = Path('figures')
    fig_dir.mkdir(parents=True, exist_ok=True)

    # List of schemes
    schemes = ['FOU', 'SOU', 'Fromm', 'Quick', 'Third']
    cfls = [0.5, 1.0]
    phi = np.linspace(-np.pi, np.pi, 201)
    for icfl, cfl in enumerate(cfls):
        fig, ax = plt.subplots()
        for scheme in schemes:
            Omega_dt = fourier_symbol(phi, cfl, scheme)
            ax.plot(Omega_dt.real, Omega_dt.imag, label=scheme)
        ax_prop(ax, 'Re($\Omega \Delta t$)', 'Im($\Omega \Delta t$)', f'CFL = {cfl:.2f}')
        fig.savefig(fig_dir / f'cfl_{icfl:d}', bbox_inches='tight')
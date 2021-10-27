from rungekutta import ax_plot_rk_stab
from schemes import fourier_symbol, ax_prop
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from cycler import cycler

default_cycler = (cycler(color=['lightcoral', 'firebrick', 'darkred']) +
                  cycler(linestyle=['-', '--', ':']))

plt.rc('lines', linewidth=1.8)
plt.rc('axes', prop_cycle=default_cycler)

if __name__ == '__main__':
    fig_dir = Path('figures/schemes_rk')
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Plot third order schemes cfls 0.5/1.0 - RK
    schemes = ['SOU', 'Fromm', 'Quick', 'Third']
    for scheme in schemes:
        cfls = [0.5, 1.0, 1.5]
        phi = np.linspace(-np.pi, np.pi, 201)
        fig, ax = plt.subplots()
        for icfl, cfl in enumerate(cfls):
            Omega_dt = fourier_symbol(phi, cfl, scheme)
            ax.plot(Omega_dt.real, Omega_dt.imag, label=f'CFL = {cfl:.1f}')
        lambda_dt_re = np.linspace(-3.5, 1.0, 201)
        lambda_dt_im = np.linspace(-3.0, 3.0, 201)
        X, Y = np.meshgrid(lambda_dt_re, lambda_dt_im)
        ax_plot_rk_stab(ax, X, Y)
        ax.legend()
        fig.savefig(fig_dir / f'{scheme}_order_rks', bbox_inches='tight')
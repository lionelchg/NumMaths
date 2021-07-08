import numpy as np
import matplotlib.pyplot as plt
from test_funcs import gaussian, step, packet_wave
import sp_analysis
from sp_analysis import ampl_factor


def plot_sim(x_th, x, x0, u_gauss, u_step, u_2pw, u_4pw, schemes, figtitle, figname):
    """ Plot the results of the finite difference schemes """
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(9, 9), sharex=True)

    axes = axes.reshape(-1)

    axes[0].plot(x_th, gaussian(x_th, x0, 0.3))
    for i_scheme, scheme in enumerate(schemes):
        axes[0].plot(x, u_gauss[i_scheme, :], '.', ms=4, label=scheme)
    ax_prop_sim(axes[0])

    axes[1].plot(x_th, step(x_th, x0))
    for i_scheme, scheme in enumerate(schemes):
        axes[1].plot(x, u_step[i_scheme, :], '.', ms=4, label=scheme)
    ax_prop_sim(axes[1])

    axes[2].plot(x_th, packet_wave(x_th, x0, 0.5))
    for i_scheme, scheme in enumerate(schemes):
        axes[2].plot(x, u_2pw[i_scheme, :], '.', ms=4, label=scheme)
    ax_prop_sim(axes[2])

    axes[3].plot(x_th, packet_wave(x_th, x0, 0.25))
    for i_scheme, scheme in enumerate(schemes):
        axes[3].plot(x, u_4pw[i_scheme, :], '--', ms=4, label=scheme)
    ax_prop_sim(axes[3])

    fig.suptitle(figtitle)
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    fig.savefig(figname, bbox_inches='tight')
    plt.close(fig)

def ax_prop_sim(ax):
    ax.legend()
    ax.grid(True)

def plot_G(scheme, cfls, fig_dir):
    """ Plot the diffusion and dispersion errors from the amplification factor """
    phi = np.linspace(0, np.pi, 300)
    phi_deg = phi * 180 / np.pi    
    fig, axes = plt.subplots(ncols=2, figsize=(8, 5))
    for cfl in cfls:
        # df_err, dp_err = sp_analysis.errors(getattr(sp_analysis, f'G_{scheme}'), phi, cfl)
        df_err, dp_err = sp_analysis.errors(ampl_factor(phi, cfl, scheme), phi, cfl)
        axes[0].plot(phi_deg, df_err, label=f'CFL = {cfl:.2f}')
        axes[1].plot(phi_deg, dp_err, label=f'CFL = {cfl:.2f}')
    ax_prop_G(axes[0], r'$\varepsilon_D$')
    ax_prop_G(axes[1], r'$\varepsilon_\phi$')
    fig.suptitle(f'{scheme} Spectral Analysis')
    fig.savefig(fig_dir / f'errors_{scheme}', bbox_inches='tight')

def ax_prop_G(ax, ylabel, ylim=None):
    """ Ax properties of plot_G """
    ax.legend()
    ax.grid(True)
    ax.set_xlabel(r'$\phi$')
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=0)
    ax.set_xlim([0, 180])

def plot_cvg(nnxs, schemes, functions, errors, figtitle, figname):
    """ Plot convergence of the solution to the exact one for schemes and functions 
    given """
    fig, axes = plt.subplots(ncols=len(functions), figsize=(5 * len(functions), 8))
    for i_func, function in enumerate(functions):
        axes[i_func].set_title(function)
        for i_scheme, scheme in enumerate(schemes):
            axes[i_func].plot(nnxs, errors[:, i_scheme, i_func], label=f'{scheme}', lw=1.5)
        ax_prop_cvg(axes[i_func])
    fig.suptitle(figtitle)
    fig.savefig(figname, bbox_inches='tight')
    plt.close()

def ax_prop_cvg(ax):
    """ Ax properties of plot_cvg """
    ax.grid(True)
    ax.set_xlabel('$n_x$')
    ax.set_ylabel(r'$\varepsilon_1$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
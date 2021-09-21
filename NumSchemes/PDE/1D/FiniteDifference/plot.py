import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from test_funcs import gaussian, step, packet_wave
import sp_analysis
from sp_analysis import ampl_factor

def plot_sim(x_th, x, x0, u_gauss, u_step, u_2pw, u_4pw, schemes, figtitle, figname):
    """ Plot the results of the finite difference schemes """
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(9, 9), sharex=True)

    axes = axes.reshape(-1)

    axes[0].plot(x_th, gaussian(x_th, x0, 0.15))
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

def plot_G(scheme, cfls, fig_dir):
    """ Plot the diffusion and dispersion errors from the amplification factor """
    phi = np.linspace(0, np.pi, 300)
    phi_deg = phi * 180 / np.pi    
    fig, axes = plt.subplots(ncols=3, figsize=(10, 4))
    for cfl in cfls:
        df_err, dp_err = sp_analysis.errors(ampl_factor(phi, cfl, scheme), phi, cfl)
        axes[0].plot(phi_deg, df_err, label=f'CFL = {cfl:.2f}')
        axes[1].plot(phi_deg, dp_err, label=f'CFL = {cfl:.2f}')
        axes[2].plot(phi_deg, dp_err * phi_deg, label=f'CFL = {cfl:.2f}')
    ax_prop_G(axes[0], r'$\varepsilon_D$')
    ax_prop_G(axes[1], r'$\varepsilon_\phi$')
    ax_prop_G(axes[2], r'$\phi_\mathrm{num}$')
    axes[2].legend()
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    fig.suptitle(f'{scheme} Spectral Analysis')
    fig.savefig(fig_dir / f'errors_1D_{scheme}', bbox_inches='tight')

def ax_prop_G(ax, ylabel, ylim=None):
    """ Ax properties of plot_G """
    ax.grid(True)
    ax.set_xlabel(r'$\phi$')
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=0)
    ax.set_xlim([0, 180])

def round_up(number: float, decimals:int=1) -> float:
    """ Rounding of a number

    :param number: A number to round
    :type number: float
    :param decimals: Number of decimals defaults to 1
    :type decimals: int, optional
    :return: The rounded number
    :rtype: float
    """

    if number != 0.0:
        power = int(np.log10(number))
        digit = number / 10**power * 10**decimals
        return np.ceil(digit) * 10**(power - decimals)
    else:
        return 1.0

def plot_ax_scalar(fig, ax, X, Y, field, title):
    """ Plot a 2D field on mesh X and Y with contourf and contour. Automatic
    handling of maximum values for the colorbar with an up rounding to a certain
    number of decimals. """
    # If not specified find the maximum value and round this up to 1 decimal
    max_value = round_up(np.max(np.abs(field)), decimals=1)

    # Depending on the scale (log is typically for streamers) the treatment
    # is not the same
    field_ticks = np.linspace(0, max_value, 5)
    levels = np.linspace(0, max_value, 101)

    cs1 = ax.contourf(X, Y, field, levels, cmap='Blues')

    # Contours
    clevels = np.linspace(np.min(field), np.max(field), 6)[1:-1]
    cs2 = ax.contour(X, Y, field, levels=clevels, colors='k', linewidths=0.9)
    ax.clabel(cs2, fmt='%.2f')
    cs2 = ax.contour(X, Y, field, levels=[1.0], colors='y', linewidths=0.9)
    ax.clabel(cs2, fmt='%.2f')

    # Adjust the size of the colorbar
    xmax, ymax = np.max(X), np.max(Y)
    fraction_cbar = 0.1
    aspect = 7
    # Set the colorbar in scientific notation
    sfmt = mpl.ticker.ScalarFormatter(useMathText=True) 
    sfmt.set_powerlimits((0, 0))
    fig.colorbar(cs1, ax=ax, pad=0.05, fraction=fraction_cbar, aspect=aspect, 
        ticks=field_ticks, format=sfmt)

    # Apply same formatting to x and y axis with scientific notation
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

    ax.set_title(title)

def plot_G_2D(scheme, fig_dir):
    """ Plot the diffusion and dispersion errors from the amplification factor
    as 2D contourf and contours """
    phi_1d = np.linspace(0, np.pi, 101)
    cfl_1d = np.linspace(0.01, 1, 51)
    phi, cfl = np.meshgrid(phi_1d, cfl_1d)
    phi_deg = phi * 180 / np.pi    
    df_err = np.zeros_like(phi)
    dp_err = np.zeros_like(phi)

    # Compute the diffusion and dispersion error
    for i in range(len(cfl_1d)):
        df_err[i, :], dp_err[i, :] = sp_analysis.errors(ampl_factor(phi_1d, 
                    cfl_1d[i], scheme), phi_1d, cfl_1d[i])
    
    # 2D plotting
    fig, axes = plt.subplots(ncols=2, figsize=(10, 4))
    plot_ax_scalar(fig, axes[0], phi_deg, cfl, df_err, "Diffusion Error")
    plot_ax_scalar(fig, axes[1], phi_deg, cfl, dp_err, "Dispersion Error")
    fig.savefig(fig_dir / f'errors_2D_{scheme}', bbox_inches='tight')
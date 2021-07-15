import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.patches import Circle
import cmath
from plot import plot_ax_scalar

def G_LW(cfl, cos_alpha, sin_alpha, phi_X, phi_Y):
    """ 2D Lax-Wendroff amplification factor """
    return (1 - 1j * cfl * (cos_alpha * np.sin(phi_X) + sin_alpha * np.sin(phi_Y))
        + cfl**2 * (cos_alpha**2 * (np.cos(phi_X) - 1) 
            + sin_alpha**2 * (np.cos(phi_Y) - 1) 
            + cos_alpha * sin_alpha * (np.cos(phi_X + phi_Y) - np.cos(phi_X - phi_Y))))

def plot_G_2D(scheme, cfl, figname):
    """ Plot 2D amplification factor """
    fig, axes = plt.subplots(ncols=2, figsize=(9, 5))
    # Compute amplification factor, diffusion and dispersion error
    tmp_G = eval(f'G_{scheme}(cfl, cos_alpha, sin_alpha, phi_X, phi_Y)')
    mod_G = abs(tmp_G)
    phi_G = np.array([[- cmath.phase(tmp_G[j, i]) for i in range(nnphi)] for j in range(nnphi)])
    disp_err = phi_G / cfl / (phi_X * cos_alpha + phi_Y * sin_alpha)

    # Plot quantities
    plot_ax_scalar(fig, axes[0], phi_X, phi_Y, mod_G, 
        f'CFL = {cfl:.2f}', cmap='Blues', max_value=1.2)
    plot_ax_scalar(fig, axes[1], phi_X, phi_Y, disp_err, 
        f'CFL = {cfl:.2f}', cmap='Blues', max_value=1.1)

    # Print circles for delimitation
    axes[0].add_patch(Circle((0, 0), np.pi, fill=False, color='yellow'))
    axes[1].add_patch(Circle((0, 0), np.pi, fill=False, color='yellow'))
    fig.savefig(figname, bbox_inches='tight')

if __name__ == '__main__':
    # Figures directory
    fig_dir = Path('figures')
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Schemes
    schemes = ['LW']

    # Creation of vectors and grid for 2D plotting
    cfls = [0.1, 0.3, 0.5, 0.7, 0.9]
    nnphi = 300
    phi_x = np.linspace(-np.pi, np.pi, nnphi)
    phi_y = np.linspace(-np.pi, np.pi, nnphi)
    phi_X, phi_Y = np.meshgrid(phi_x, phi_y)
    phi_norm = np.sqrt(phi_X**2 + phi_Y**2)
    cos_alpha = phi_X / phi_norm
    sin_alpha = phi_Y / phi_norm

    # Figure drawing
    for scheme in schemes:
        scheme_dir = fig_dir / scheme
        scheme_dir.mkdir(parents=True, exist_ok=True)
        for i_cfl, cfl in enumerate(cfls):
            plot_G_2D('LW', cfl, scheme_dir / f'cfl_{i_cfl}')
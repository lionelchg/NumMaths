import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def axis_center(ax):
    # Move left y-axis and bottim x-axis to centre, passing through (0,0)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')

    # Eliminate upper and right axes
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    # Show ticks in the left and lower axes only
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

def plot_stability_rk(X, Y, figtitle, figname):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax_plot_rk_stab(ax, X, Y)
    ax.set_title(figtitle)
    fig.savefig(figname, bbox_inches='tight')

def ax_plot_rk_stab(ax, X, Y):
    Omega_dt = X + 1j * Y
    Z_rk1 = np.abs(1 + Omega_dt)
    ax.contour(X, Y, Z_rk1, levels=[1.0], colors='lightsteelblue', linestyles='solid')
    ax.text(-1, 0.75, 'RK1')
    Z_rk2 = np.abs(1 + Omega_dt * (1 + 0.5 * Omega_dt))
    ax.contour(X, Y, Z_rk2, levels=[1.0], colors='royalblue', linestyles='dashed')
    ax.text(-0.9, 1.4, 'RK2')
    Z_rk3 = np.abs(1 + Omega_dt * (1 + 0.5 * Omega_dt * (1 + 1 / 3 * Omega_dt)))
    ax.contour(X, Y, Z_rk3, levels=[1.0], colors='mediumblue', linestyles='dashdot')
    ax.text(-1.1, 2.1, 'RK3')
    Z_rk3 = np.abs(1 + Omega_dt * (1 + 0.5 * Omega_dt * (1 + 0.5 * Omega_dt)))
    ax.contour(X, Y, Z_rk3, levels=[1.0], colors='mediumblue', linestyles='solid')
    ax.text(0.3, 1.5, 'RK3-Stab')
    Z_rk4 = np.abs(1 + Omega_dt * (1 + 0.5 * Omega_dt * (1 + 1 / 3 * Omega_dt * (1 + 1 / 4 * Omega_dt))))
    ax.contour(X, Y, Z_rk4, levels=[1.0], colors='darkblue', linestyles='dotted')
    ax.text(-0.55, 2.65, 'RK4')
    ax.set_aspect('equal')
    ax.set_xlabel(r'Re$(\Omega \Delta t)$')
    ax.set_ylabel(r'Im$(\Omega \Delta t)$')
    ax.grid(True)

if __name__ == '__main__':
    fig_dir = Path('figures/stability/')
    fig_dir.mkdir(parents=True, exist_ok=True)

    # RK1
    lambda_dt_re = np.linspace(-3.5, 1.0, 201)
    lambda_dt_im = np.linspace(-3.0, 3.0, 201)
    X, Y = np.meshgrid(lambda_dt_re, lambda_dt_im)
    plot_stability_rk(X, Y, 'Runge-Kutta Stability limits', fig_dir / 'runge_kutta_stab')


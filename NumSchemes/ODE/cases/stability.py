import numpy as np
import matplotlib.pyplot as plt

from odesolver.utils import create_dir

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

def plot_stability(X, Y, Z, figtitle, figname):
    fig, ax = plt.subplots()
    ax_plot_contour(ax, X, Y, Z)
    ax.set_title(figtitle)
    fig.savefig(figname, bbox_inches='tight')

def ax_plot_contour(ax, X, Y, Z):
    cs = ax.contour(X, Y, Z, levels=[0.5, 0.75, 1.0, 2.0])
    ax.clabel(cs, inline=True, fontsize=10)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\lambda_r \Delta t$')
    ax.set_ylabel(r'$\lambda_i \Delta t$')

if __name__ == '__main__':
    fig_dir = 'figures/stability/'
    create_dir(fig_dir)
    
    # Backward euler
    lambda_dt_re = np.linspace(-1.5, 3.5, 201)
    lambda_dt_im = np.linspace(-2.5, 2.5, 201)
    X, Y = np.meshgrid(lambda_dt_re, lambda_dt_im)
    Z = np.abs(1 / (1 - (X + 1j * Y)))
    plot_stability(X, Y, Z, 'Backward Euler', fig_dir + 'backward_euler')

    # Second order BDF
    alpha_1, alpha_2, beta_0 = - 4 / 3, 1 / 3, 2 / 3
    lambda_dt_re = np.linspace(-3, 5, 401)
    lambda_dt_im = np.linspace(-4, 4, 401)
    X, Y = np.meshgrid(lambda_dt_re, lambda_dt_im)
    root1 = (- alpha_1 - np.sqrt(alpha_1**2 - 4 * alpha_2 * (1 - (X + 1j * Y) * beta_0))) \
                / (2 * (1 - (X + 1j * Y) * beta_0))
    root2 = (- alpha_1 + np.sqrt(alpha_1**2 - 4 * alpha_2 * (1 - (X + 1j * Y) * beta_0))) \
                / (2 * (1 - (X + 1j * Y) * beta_0))
    g = np.maximum(np.abs(root1), np.abs(root2))
    plot_stability(X, Y, g, 'BDF2', fig_dir + 'bdf2')

    # Third order BDF
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

from cycler import cycler

default_cycler = (cycler(color=['darkblue', 'mediumblue', 'royalblue']) +
                  cycler(linestyle=['-', '--', ':']))

plt.rc('lines', linewidth=1.8)
plt.rc('axes', prop_cycle=default_cycler)

def van_leer(r):
    return (r + np.abs(r)) / (1 + r)

def min_mod(r):
    return np.where(r >= 0, np.minimum(r, np.ones_like(r)), np.zeros_like(r))

def superbee(r):
    zeros = np.zeros_like(r)
    ones = np.ones_like(r)
    return np.maximum.reduce([zeros, np.minimum(2 * r, ones), np.minimum(r, 2 * ones)])

def beta_lim(r, beta):
    zeros = np.zeros_like(r)
    ones = np.ones_like(r)
    return np.maximum.reduce([zeros, np.minimum(beta * r, ones), np.minimum(r, beta * ones)])

def osher(r):
    zeros = np.zeros_like(r)
    ones = np.ones_like(r)
    return np.maximum(zeros, np.minimum(r, 2 * ones))

def alpha_lim(r, alpha):
    zeros = np.zeros_like(r)
    ones = np.ones_like(r)
    return np.maximum(zeros,
        np.minimum.reduce([2 * r, alpha * r + (1 - alpha), 2 * ones]))

def ax_prop(ax, title):
    ax.grid(True)
    ax.set_title(title)
    ax.set_xlim([0, 3])
    ax.set_ylim([0, 2.5])
    ax.set_xlabel('$r$')
    ax.set_ylabel('$\Psi(r)$')

if __name__ == '__main__':
    fig_dir = Path('figures/limiters')
    fig_dir.mkdir(parents=True, exist_ok=True)

    r = np.linspace(0, 3, 301)
    fig, axes = plt.subplots(nrows=3, ncols=2,
        figsize=(10, 10), sharex=True, sharey=True)
    axes = axes.reshape(-1)
    axes[0].plot(r, van_leer(r))
    ax_prop(axes[0], 'Van Leer limiter')
    axes[1].plot(r, min_mod(r))
    ax_prop(axes[1], 'Min-mod limiter')
    axes[2].plot(r, superbee(r))
    ax_prop(axes[2], 'Superbee limiter')
    axes[3].plot(r, beta_lim(r, 1.5))
    ax_prop(axes[3], r'$\beta = 1.5$-limiter')
    axes[4].plot(r, osher(r))
    ax_prop(axes[4], 'Osher limiter')
    axes[5].plot(r, alpha_lim(r, 2 / 3))
    ax_prop(axes[5], r'$\alpha = 2/3$-limiter')

    # Shaded area for stability
    stab_curve_up = np.where(r <= 1, np.minimum(2 * r, np.ones_like(r)), np.minimum(r, 2 * np.ones_like(r)))
    stab_curve_down = np.where(r <= 1, r, 1)
    for ax in axes:
        ax.fill_between(r, stab_curve_down, stab_curve_up, alpha=0.3)
    fig.savefig(fig_dir / 'limiters', bbox_inches='tight')
    plt.close(fig)

    # Specific plot for sweby limiter
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(r, beta_lim(r, 1.0), label=r'$\beta = 1$ (Min-mod)')
    ax.plot(r, beta_lim(r, 1.5), label=r'$\beta = 1.5$')
    ax.plot(r, beta_lim(r, 2.0), label=r'$\beta = 2$ (Superbee)')
    ax.fill_between(r, stab_curve_down, stab_curve_up, alpha=0.3)
    ax_prop(ax, '')
    ax.legend()
    fig.savefig(fig_dir / 'sweby_limiters.pdf', format='pdf', bbox_inches='tight')
    plt.close(fig)

    # Specific plot for vanleer limiter
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(r, van_leer(r))
    ax.fill_between(r, stab_curve_down, stab_curve_up, alpha=0.3)
    ax_prop(ax, '')
    fig.savefig(fig_dir / 'vanleer_limiter.pdf', format='pdf', bbox_inches='tight')
    plt.close(fig)
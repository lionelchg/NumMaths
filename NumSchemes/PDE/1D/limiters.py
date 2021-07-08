import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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
    ax.set_title(title)
    ax.set_xlim([0, 3])
    ax.set_ylim([0, 2.5])

if __name__ == '__main__':
    fig_dir = 'figures/'
    r = np.linspace(0, 3, 301)
    sns.set_theme()
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
    fig.savefig(f'{fig_dir}limiters', bbox_inches='tight')
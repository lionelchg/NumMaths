import numpy as np
from utils import compute_eigen
from numpy import linalg as LA
import matplotlib.pyplot as plt
from pathlib import Path

# Switch the signs of matrices so that it is - \nabla^2 \phi that is solved
# and not \nabla^2 \phi, is there a possible gain in performance having the eigenvalues
# with the same sign?

def poisson_1D_cart(n: int):
    """ Create the array of poisson dirichlet 1D problem """
    mat = np.zeros((n, n))
    mat[0, 0] = 2
    mat[0, 1] = - 2
    mat[-1, -1] = 1
    for i in range(1, n - 1):
        mat[i, i - 1] = - 1
        mat[i, i] = 2
        mat[i, i + 1] = - 1
    return mat

def poisson_1D_axi(n: int):
    """ Create the array of poisson dirichlet 1D problem """
    mat = np.zeros((n, n))
    mat[0, 0] = 4
    mat[0, 1] = - 4
    mat[-1, -1] = 1
    for i in range(1, n - 1):
        mat[i, i - 1] = - (i - 0.5) / i
        mat[i, i] = 2
        mat[i, i + 1] = - (i + 0.5) / i
    return mat

def ax_prop(ax, xlabel, ylabel):
    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

def eigen_matrices(nn: int, fig_dir: Path):
    fig, ax = plt.subplots(figsize=(9, 3), sharex=True)
    ax.grid(True)

    # Cartesian 1D
    mat_dir = poisson_1D_cart(nn)
    eigenval, eigenvect = compute_eigen(mat_dir, lprint=False)
    ax.scatter(np.real(eigenval), np.imag(eigenval), c="darkblue", label='Cartesian')

    # Cylindric 1D
    mat_dir = poisson_1D_axi(nn)
    eigenval, eigenvect = compute_eigen(mat_dir, lprint=False)
    ax.scatter(np.real(eigenval), np.imag(eigenval) + 0.3, c="darkred", label='Cylindrical')
    ax_prop(ax, r'$\mathrm{Re}(\lambda)$', r'')
    ax.set_ylim([-0.2, 0.5])
    ax.yaxis.set_ticklabels([])

    fig.savefig(fig_dir / f'eigenvalues_{nn:04d}', bbox_inches='tight')
    fig.savefig(fig_dir / f'eigenvalues_{nn:04d}.pdf', format='pdf', bbox_inches='tight')
    plt.close(fig)

def ax_prop_cond(ax, xlabel, ylabel):
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()

if __name__ == '__main__':
    # Figures directory
    fig_dir = Path('figures/1D_cart_vs_axi')
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Print info on matrices
    nns = [11, 21]
    for nn in nns:
        eigen_matrices(nn, fig_dir)

    nns = [11, 21, 41, 81, 161, 401, 801]
    norm_types = [1, 2, np.inf]
    norm_strs = ['1', '2', 'inf']
    K = dict()
    K['cart'] = dict()
    K['cyl'] = dict()
    for norm_str in norm_strs:
        K['cart'][norm_str] = np.zeros(len(nns))
        K['cyl'][norm_str] = np.zeros(len(nns))
    for inn, nn in enumerate(nns):
        mat_cart = poisson_1D_cart(nn)
        mat_cyl = poisson_1D_axi(nn)
        for norm_str, norm in zip(norm_strs, norm_types):
            K['cart'][norm_str][inn] = LA.cond(mat_cart, norm)
            K['cyl'][norm_str][inn] = LA.cond(mat_cyl, norm)

    fig, ax = plt.subplots()
    ax.plot(nns, K['cart']['2'], '--', color='darkblue', label='Cartesian')
    ax.plot(nns, K['cyl']['2'], color='darkred', label='Cylindrical')
    ax_prop_cond(ax, 'Number of nodes', 'Spectral Condition number $K_2(A)$')
    fig.savefig(fig_dir / 'condition_number', bbox_inches='tight')
    fig.savefig(fig_dir / 'condition_number.pdf', format='pdf', bbox_inches='tight')
    plt.close(fig)
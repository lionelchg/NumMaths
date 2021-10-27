import numpy as np
from utils import compute_eigen
from numpy import linalg as LA
import matplotlib.pyplot as plt
from pathlib import Path
from poisson_1D import poisson_dirichlet_1D, poisson_inner_1D, poisson_dirichlet_1D_sym, \
                    poisson_1D_axi

# Switch the signs of matrices so that it is - \nabla^2 \phi that is solved
# and not \nabla^2 \phi, is there a possible gain in performance having the eigenvalues
# with the same sign?

# def poisson_dirichlet_1D(n: int):
#     """ Create the array of poisson dirichlet 1D problem """
#     mat = np.zeros((n, n))
#     mat[0, 0] = 1.0
#     mat[-1, -1] = 1.0
#     for i in range(1, n - 1):
#         mat[i, i - 1] = - 1
#         mat[i, i] = 2
#         mat[i, i + 1] = - 1
#     return mat

# def poisson_dirichlet_1D_sym(n: int):
#     """ Create the array of poisson dirichlet 1D problem with symmetric matrix """
#     mat = np.zeros((n, n))
#     mat[0, 0] = 1.0
#     mat[-1, -1] = 1.0
#     for i in range(1, n - 1):
#         mat[i, i - 1] = - 1
#         mat[i, i] = 2
#         mat[i, i + 1] = - 1
#     mat[1, 0] = 0.0
#     mat[-2, -1] = 0.0
#     return mat

# def poisson_1D_axi(n: int):
#     """ Create the array of poisson dirichlet 1D problem """
#     mat = np.zeros((n, n))
#     mat[0, 0] = 4
#     mat[0, 1] = - 4
#     mat[-1, -1] = 1
#     for i in range(1, n - 1):
#         mat[i, i - 1] = - (i - 0.5) / i
#         mat[i, i] = 2
#         mat[i, i + 1] = - (i + 0.5) / i
#     return mat

# def poisson_inner_1D(n: int):
#     """ Create the array of inner symmetric Poisson matrix in 1D """
#     mat = np.zeros((n, n))
#     mat[0, 0], mat[0, 1] = 2, - 1
#     mat[-1, -1], mat[-1, -2] = 2, - 1
#     for i in range(1, n - 1):
#         mat[i, i - 1] = - 1
#         mat[i, i] = 2
#         mat[i, i + 1] = - 1
#     return mat

def ax_prop(ax, xlabel, ylabel):
    ax.grid(True)
    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

def eigen_matrices(nn: int, fig_dir: Path):
    fig, axes = plt.subplots(nrows=4, figsize=(9, 7), sharex=True)

    # Inner part of the matrix
    mat_dir = poisson_inner_1D(nn)
    eigenval, eigenvect = compute_eigen(mat_dir, lprint=False)
    axes[0].scatter(np.real(eigenval), np.imag(eigenval), c="darkblue", label='Inner Poisson')
    ax_prop(axes[0], '', r'$\mathrm{Im}$')

    # Dirichlet non-symmetric
    mat_dir = poisson_dirichlet_1D(nn)
    eigenval, eigenvect = compute_eigen(mat_dir, lprint=False)
    axes[1].scatter(np.real(eigenval), np.imag(eigenval), c="mediumblue", label='Dirichlet Non-Sym')
    ax_prop(axes[1], '', r'$\mathrm{Im}$')

    # Dirichlet symmetric
    mat_dir = poisson_dirichlet_1D_sym(nn)
    eigenval, eigenvect = compute_eigen(mat_dir, lprint=False)
    axes[2].scatter(np.real(eigenval), np.imag(eigenval), c="lightsteelblue", label='Dirichlet Sym')
    ax_prop(axes[2], '', r'$\mathrm{Im}$')

    # Dirichlet non-symmetric
    mat_dir = poisson_1D_axi(nn)
    eigenval, eigenvect = compute_eigen(mat_dir, lprint=False)
    axes[3].scatter(np.real(eigenval), np.imag(eigenval), c="darkred", label='Axisymmetric')
    ax_prop(axes[3], r'$\mathrm{Re}$', r'$\mathrm{Im}$')

    fig.suptitle(f'N = {nn:d}')
    fig.savefig(fig_dir / f'eigenvalues_{nn:04d}', bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    # Figures directory
    fig_dir = Path('figures/eigenvalues')
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Print info on matrices
    nns = [11, 21, 41]
    for nn in nns:
        eigen_matrices(nn, fig_dir)
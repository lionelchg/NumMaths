import numpy as np
from utils import show_matrix, info_matrix, show_eigen
from iterative import rho_method, spectal_radius

def poisson_dirichlet_1D(n: int):
    """ Create the array of poisson dirichlet 1D problem """
    mat = np.zeros((n, n))
    mat[0, 0] = 1.0
    mat[-1, -1] = 1.0
    for i in range(1, n - 1):
        mat[i, i - 1] = 1
        mat[i, i] = -2
        mat[i, i + 1] = 1
    return mat

def poisson_dirichlet_1D_sym(n: int):
    """ Create the array of poisson dirichlet 1D problem with symmetric matrix """
    mat = np.zeros((n, n))
    mat[0, 0] = 1.0
    mat[-1, -1] = 1.0
    for i in range(1, n - 1):
        mat[i, i - 1] = 1
        mat[i, i] = -2
        mat[i, i + 1] = 1
    mat[1, 0] = 0.0
    mat[-2, -1] = 0.0
    return mat

def poisson_1D_axi(n: int):
    """ Create the array of poisson dirichlet 1D problem """
    mat = np.zeros((n, n))
    mat[0, 0] = -4
    mat[0, 1] = 4
    mat[-1, -1] = 1
    for i in range(1, n - 1):
        mat[i, i - 1] = (i - 0.5) / i
        mat[i, i] = -2
        mat[i, i + 1] = (i + 0.5) / i
    return mat

if __name__ == '__main__':
    # Iterative methods
    jacobi_method = {'name': 'decomp_jacobi', 'args': {}}
    gs_method = {'name': 'decomp_gauss_seidel', 'args': {}}

    # Print info on matrices
    print('Matrices for n = 11')
    show_matrix(poisson_dirichlet_1D(11), 'Dirichlet')
    show_matrix(poisson_dirichlet_1D_sym(11), 'DirichletSym')
    show_matrix(poisson_1D_axi(11), 'Axi')
    for n in [11, 21, 51, 101]:
        print('-------------')
        print(f'n = {n:d}')
        mat = poisson_dirichlet_1D(n)
        mat_sym = poisson_dirichlet_1D_sym(n)
        mat_axi = poisson_1D_axi(n)

        # Print info on matrices
        print('Base matrix')
        info_matrix(mat)
        print('Symmetric matrix')
        info_matrix(mat_sym)
        print('Aximatrix:')
        info_matrix(mat_axi)

        # Look at the spectral radius of iteration matrices
        print(f'rho_base = {spectal_radius(mat):.2e}')
        print(f'rho_jacobi = {rho_method(mat, jacobi_method):.2e}')
        print(f'rho_gauss = {rho_method(mat, gs_method):.2e}')
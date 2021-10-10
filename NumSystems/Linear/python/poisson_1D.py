import numpy as np
from utils import show_matrix, info_matrix, show_eigen
from iterative import rho_method, spectal_radius
from gradient_method import conjugate_gradient_method

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
    n = 11
    mat_dir = poisson_dirichlet_1D(n)
    mat_dir_sym = poisson_dirichlet_1D_sym(n)
    mat_axi = poisson_1D_axi(n)

    # Solve ficitous system
    x_sol = np.ones(n)
    print('Dirichlet Symmetric')
    b = np.matmul(mat_dir_sym, x_sol)
    x_0 = np.linspace(1, 11, 11)
    conjugate_gradient_method(mat_dir_sym, x_0, b, 20, 1e-3)

    print('Dirichlet Non-Symmetric')
    b = np.matmul(mat_dir, x_sol)
    x_0 = np.linspace(1, 11, 11)
    conjugate_gradient_method(mat_dir, x_0, b, 20, 1e-3)

    print('Axi matrix')
    b = np.matmul(mat_axi, x_sol)
    x_0 = np.linspace(1, 11, 11)
    conjugate_gradient_method(mat_axi, x_0, b, 20, 1e-3)

    show_matrix(mat_dir, 'Dirichlet')
    show_matrix(mat_dir_sym, 'DirichletSym')
    show_matrix(mat_axi, 'Axi')

    # Print info about matrices for different sizes
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
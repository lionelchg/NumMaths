import numpy as np
from utils import show_matrix, info_matrix, show_eigen

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
    mat[0, 0], mat[0, 1] = -2, 1
    mat[-1, -2], mat[-1, -1] = 1, -2
    for i in range(1, n - 1):
        mat[i, i - 1] = 1
        mat[i, i] = -2
        mat[i, i + 1] = 1
    return mat

if __name__ == '__main__':
    print('Matrices for n = 11')
    show_matrix(poisson_dirichlet_1D(11))
    show_eigen(poisson_dirichlet_1D(11))
    show_matrix(poisson_dirichlet_1D_sym(9))
    for n in [11, 21, 51, 101]:
        print('-------------')
        print(f'n = {n:d}')
        mat = poisson_dirichlet_1D(n)
        mat_sym = poisson_dirichlet_1D_sym(n - 2)
        print('Base matrix')
        info_matrix(mat)
        print('Symmetric matrix')
        info_matrix(mat_sym)

    
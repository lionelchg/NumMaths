import numpy as np
from numpy import linalg as LA

def print_mat(A, matname):
    print(f'Matrix {matname}')
    n = len(A[:, 0])
    for i in range(n):
        row = ''
        for j in range(n):
            row += f'{A[i, j]:>6.2f} '
        print(row)

def show_matrix(mat):
    """ Print matrix and inverse of matrix """
    mat_inv = LA.inv(mat)
    print_mat(mat, 'A')
    print_mat(mat_inv, 'A^-1')

def info_matrix(mat):
    """ Print determinant and condition numbers of the matrix """
    det = LA.det(mat)
    print(f'|A| = {det:.3e}')
    for norm in [1, 2, np.inf]:
        print(f'K_{norm:.0f} = {LA.cond(mat, norm):.2e}')

def show_eigen(mat):
    """ Print eigenvalues and eigenvectors of a given matrix """
    w, v = LA.eig(mat)
    print('Eigenvalues:')
    print(w)
    print('Eigenvectors:')
    print(v)

def vec_str(vec):
    """ Print vector in readable format """
    print_str = ''
    for scalar in vec:
        print_str += f'{scalar:10.2e}'
    return print_str
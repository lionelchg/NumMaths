import numpy as np
from numpy import linalg

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
    mat_inv = linalg.inv(mat)
    print_mat(mat, 'A')
    print_mat(mat_inv, 'A^-1')

def info_matrix(mat):
    """ Print determinant and condition numbers of the matrix """
    det = linalg.det(mat)
    print(f'|A| = {det:.3e}')
    for norm in [1, 2, np.inf]:
        print(f'K_{norm:.0f} = {linalg.cond(mat, norm):.2e}')

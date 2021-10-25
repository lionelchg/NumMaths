import numpy as np
from numpy import linalg as LA

def print_mat(A, matname):
    print(f'Matrix {matname}')
    n = len(A[:, 0])
    for i in range(n):
        row = ''
        for j in range(n):
            # Depending on the magnitude of the entry
            # use scientific notation or not
            if abs(np.log10(abs(A[i, j]))) < 3 or A[i, j] == 0.0:
                row += f'{A[i, j]:>6.2f} '
            else:
                row += f'{A[i, j]:>6.2e} '

        print(row)

def show_matrix(mat, name_mat='A'):
    """ Print matrix and inverse of matrix """
    mat_inv = LA.inv(mat)
    print_mat(mat, name_mat)
    print_mat(mat_inv, f'{name_mat}^-1')

def info_matrix(mat):
    """ Print determinant and condition numbers of the matrix """
    det = LA.det(mat)
    print(f'|A| = {det:.3e}')
    for norm in [1, 2, np.inf]:
        print(f'K_{norm:.0f} = {LA.cond(mat, norm):.2e}')

def compute_eigen(mat, lprint=True):
    """ Print eigenvalues and eigenvectors of a given matrix """
    eigenvalues, eigenvectors = LA.eig(mat)
    if lprint:
        print('Eigenvalues:')
        print(eigenvalues)
        print('Eigenvectors:')
        print(eigenvectors)
    return eigenvalues, eigenvectors

def vec_str(vec):
    """ Print vector in readable format """
    print_str = ''
    for scalar in vec:
        print_str += f'{scalar:10.2e}'
    return print_str
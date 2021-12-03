########################################################################################################################
#                                                                                                                      #
#                                            2D Poisson solver using numpy                                             #
#                                                                                                                      #
#                                          Lionel Cheng, CERFACS, 10.03.2020                                           #
#                                                                                                                      #
########################################################################################################################
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
from utils import info_matrix, show_matrix
import yaml
import copy
import numpy as np
from scipy import sparse

from PlasmaNet.poissonsolver.poisson import PoissonLinSystem

def matrix_cart_tri(nx, nr):
    """ Build the matrix for the axisymmetric configuration. """
    diags = np.zeros((5, nx * nr))

    # Filling the diagonals, first the down neumann bc, then the dirichlet bc and finally the interior nodes
    for i in range(nx * nr):
        if 0 < i < nx - 1 or i >= (nr - 1) * nx or i % nx == 0 or i % nx == nx - 1:
            diags[0, i] = 1
            diags[1, min(i + 1, nx * nr - 1)] = 0
            diags[2, max(i - 1, 0)] = 0
            diags[3, min(i + nx, nx * nr - 1)] = 0
            diags[4, max(i - nx, 0)] = 0
        else:
            diags[0, i] = - 4
            diags[1, i + 1] = 1
            diags[2, i - 1] = 1
            diags[3, i + nx] = 1
            diags[4, i - nx] = 1

    # Creating the matrix
    return sparse.csc_matrix(
        sparse.dia_matrix((diags, [0, 1, -1, nx, -nx]), shape=(nx * nr, nx * nr)))

def matrix_cart_quad(nx, nr):
    """ Build the matrix for the axisymmetric configuration. """
    diags = np.zeros((5, nx * nr))

    # Filling the diagonals, first the down neumann bc, then the dirichlet bc and finally the interior nodes
    for i in range(nx * nr):
        if 0 < i < nx - 1 or i >= (nr - 1) * nx or i % nx == 0 or i % nx == nx - 1:
            diags[0, i] = 1
            diags[1, min(i + 1, nx * nr - 1)] = 0
            diags[2, max(i - 1, 0)] = 0
            diags[3, min(i + nx, nx * nr - 1)] = 0
            diags[4, max(i - nx, 0)] = 0
        else:
            diags[0, i] = - 2
            diags[1, i + 1 + nx] = 0.5
            diags[2, i + 1  - nx] = 0.5
            diags[3, i - 1 + nx] = 0.5
            diags[4, i - 1 - nx] = 0.5

    # Creating the matrix
    return sparse.csc_matrix(
        sparse.dia_matrix((diags, [0, 1 + nx, 1 - nx, -1 + nx, - 1 - nx]), shape=(nx * nr, nx * nr)))


if __name__ == '__main__':
    with open('poisson_ls_xy.yml', 'r') as yaml_stream:
        cfg_cart = yaml.safe_load(yaml_stream)

    nnxs = [5, 11, 21]
    for nnx in nnxs:
        print('------------\n'
                f'nnx = {nnx:d}')

        print('CDM K(A):')
        mat = matrix_cart_tri(nnx, nnx)
        info_matrix(mat.todense())

        mat_quad = matrix_cart_quad(nnx, nnx)
        info_matrix(mat_quad.todense())

        if nnx == 5:
            show_matrix(mat.todense(), 'Cart Tri')
            show_matrix(mat_quad.todense(), 'Cart Quad')

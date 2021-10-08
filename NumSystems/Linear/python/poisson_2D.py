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

from PlasmaNet.poissonsolver.poisson import PoissonLinSystem

if __name__ == '__main__':
    with open('poisson_ls_xy.yml', 'r') as yaml_stream:
        cfg_cart = yaml.safe_load(yaml_stream)
    cfg_axi = copy.deepcopy(cfg_cart)
    cfg_axi['geom'] = 'cylindrical'
    nnxs = [5, 11, 21]
    for nnx in nnxs:
        print('------------\n'
                f'nnx = {nnx:d}')
        cfg_cart['nnx'], cfg_cart['nny'] = nnx, nnx
        poisson = PoissonLinSystem(cfg_cart)

        cfg_axi['nnx'], cfg_axi['nny'] = nnx, nnx
        poisson_rx = PoissonLinSystem(cfg_axi)

        print('CDM K(A):')
        info_matrix(poisson.mat.todense())

        print('\nADM K(A):')
        info_matrix(poisson_rx.mat.todense())

        if nnx == 5:
            show_matrix(poisson.mat.todense(), 'Cart')
            show_matrix(poisson_rx.mat.todense(), 'Axi')

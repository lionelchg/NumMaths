import os
import numpy as np
import argparse
from numba import njit

from plot import plot_sim, plot_G, plot_cvg
from test_funcs import gaussian, step, packet_wave
from utils import create_dir
from fd_schemes import its_fd
from errors import L1error, L2error, Linferror
from pathlib import Path

def main(args):
    # figures directory
    fig_dir = Path(f'figures/{args.figdir}/')
    cvg_dir = fig_dir / 'cvg'
    sim_dir = fig_dir / 'sim'
    sp_dir = fig_dir / 'sp'

    # Create the directories
    cvg_dir.mkdir(parents=True, exist_ok=True)
    sim_dir.mkdir(parents=True, exist_ok=True)
    sp_dir.mkdir(parents=True, exist_ok=True)
    
    # Schemes selected
    schemes = args.schemes
    print(f"Schemes selected: {' '.join(schemes)}")

    # Number of periods
    n_periods = 2.0

    # Mesh properties
    xmin, xmax, nnx = -1, 1, 401
    Lx, ncx = xmax - xmin, nnx - 1
    x_th = np.linspace(xmin, xmax, 1001)
    dx = (xmax - xmin) / ncx

    # Center of the domain
    x0 = (xmax + xmin) / 2

    # Number of schemes to test
    n_schemes = len(schemes)

    # CFL and simulation properties
    a = 1.0
    cfls = [0.1, 0.3, 0.5, 0.7, 0.9]

    print('\n-------------------------------------------------------')
    print(f'Launching simulations at nnx = {nnx:d} and a = {a:.2f}')
    print('-------------------------------------------------------')

    # Launching all the cases
    for index, cfl in enumerate(cfls):
        dt = dx * cfl / a
        
        # initialization (number of timesteps required to do a full round)
        x = np.linspace(xmin, xmax, nnx)
        u_gauss = np.tile(gaussian(x, x0, 0.3), n_schemes).reshape(n_schemes,  nnx)
        u_step = np.tile(step(x, x0), n_schemes).reshape(n_schemes,  nnx)
        u_2pw = np.tile(packet_wave(x, x0, 0.5), n_schemes).reshape(n_schemes,  nnx)
        u_4pw = np.tile(packet_wave(x, x0, 0.25), n_schemes).reshape(n_schemes,  nnx)
        res = np.zeros_like(x)
        nt = int(n_periods * Lx / a / dt)

        print(f'CFL = {cfl:.2f} - nt = {nt:d}')

        # Iteration of the schemes
        for i_scheme, scheme in enumerate(schemes):
            its_fd(nt, res, u_gauss[i_scheme, :], cfl, scheme)
            its_fd(nt, res, u_step[i_scheme, :], cfl, scheme)
            its_fd(nt, res, u_2pw[i_scheme, :], cfl, scheme)
            its_fd(nt, res, u_4pw[i_scheme, :], cfl, scheme)

        # One plot per cfl
        plot_sim(x_th, x, x0, u_gauss, u_step, u_2pw, u_4pw, 
                schemes, f'CFL = {cfl:.2f} - dx = {dx:.2e} m - dt = {dt:.2e} s - nits = {nt:d}', 
                sim_dir / f'sim_cfl_{index}')

    # Plot the diffusion and disperson errors 
    # from the amplification factors of the schemes
    print(f'\n--> Plotting amplifications factors...')
    for scheme in schemes:
        plot_G(scheme, cfls, sp_dir)

    print('\n-------------------------------------------------------')
    print(f'Studying mesh convergence')
    print('-------------------------------------------------------')
    
    # Mesh convergence of the schemes
    functions = ['gaussian(x, x0, 0.3)', 'step(x, x0)', 'packet_wave(x, x0, 0.5)']
    nnxs = np.array([51, 101, 201, 501])
    errors = np.zeros((len(nnxs), n_schemes, len(functions)))
    n_periods_cvg = 2.0
    for index, cfl in enumerate(cfls):
        errors[:] = 0.0
        for i_mesh, nnx in enumerate(nnxs):
            ncx = nnx - 1
            dx = (xmax - xmin) / ncx
            dt = dx * cfl / a
            # initialization (number of timesteps required to do the
            # number of periods required)
            nt = int(n_periods_cvg * Lx / a / dt)
            x = np.linspace(xmin, xmax, nnx)
            res = np.zeros_like(x)

            print(f'CFL = {cfl:.2f} - nt = {nt:d}')
            
            for i_func, function in enumerate(functions):
                u_th = eval(function)
                u_sim = np.tile(u_th, n_schemes).reshape(n_schemes,  nnx)

                # Iteration of the schemes
                for i_scheme, scheme in enumerate(schemes):
                    its_fd(nt, res, u_sim[i_scheme, :], cfl, scheme)
                    errors[i_mesh, i_scheme, i_func] = L1error(u_th, u_sim[i_scheme, :], ncx)
        # One plot per cfl
        plot_cvg(nnxs, schemes, functions, errors, f'CFL = {cfl:.2f}', cvg_dir / f'cfl_{index}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--schemes', help='Name of the schemes to study', nargs='+')
    parser.add_argument('-d', '--figdir', help='Name of the figures directory')
    args = parser.parse_args()
    main(args)
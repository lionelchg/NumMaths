import os
import numpy as np
import argparse
from numba import njit

from plot import plot_sim, plot_G, plot_cvg
from test_funcs import gaussian, step, packet_wave
from fd_schemes import its_fd, its_limiter_fd
from errors import L1error, L2error, Linferror
from pathlib import Path

lim_list = ['van_leer', 'min_mod', 'superbee']

def run_sims(a: float, cfls: list, xmin: float, xmax: float, nnx: int, 
    n_periods: float, schemes: list, sim_dir: Path):
    """Run simulations for given cfls in 1D domain (xmin, xmax, nnx) for
    n_periods and for each scheme in schemes with figures saved in sim_dir

    :param a: Convection speed
    :type a: float
    :param cfls: Cfls studied
    :type cfls: list
    :param xmin: Left boundary
    :type xmin: float
    :param xmax: Right boundary
    :type xmax: float
    :param nnx: Number of nodes
    :type nnx: int
    :param n_periods: Number of periods
    :type n_periods: float
    :param schemes: Schemes
    :type schemes: list
    :param sim_dir: Figures directory
    :type sim_dir: Path
    """
    # Mesh properties
    Lx, ncx = xmax - xmin, nnx - 1
    x_th = np.linspace(xmin, xmax, 1001)
    dx = (xmax - xmin) / ncx
    x = np.linspace(xmin, xmax, nnx)

    # Center of the domain
    x0 = (xmax + xmin) / 2

    # Number of schemes to test
    n_schemes = len(schemes)

    print('\n-------------------------------------------------------')
    print(f'Launching simulations at nnx = {nnx:d} and a = {a:.2f}')
    print('-------------------------------------------------------')

    # Launching all the cases
    for index, cfl in enumerate(cfls):
        dt = dx * cfl / a
        
        # initialization (number of timesteps required to do a full round)
        u_gauss = np.tile(gaussian(x, x0, 0.3), n_schemes).reshape(n_schemes,  nnx)
        u_step = np.tile(step(x, x0), n_schemes).reshape(n_schemes,  nnx)
        u_2pw = np.tile(packet_wave(x, x0, 0.5), n_schemes).reshape(n_schemes,  nnx)
        u_4pw = np.tile(packet_wave(x, x0, 0.25), n_schemes).reshape(n_schemes,  nnx)
        res = np.zeros_like(x)
        nt = int(n_periods * Lx / a / dt)

        print(f'CFL = {cfl:.2f} - nt = {nt:d}')

        # Iteration of the schemes
        for i_scheme, scheme in enumerate(schemes):
            if scheme in lim_list:
                its_limiter_fd(nt, res, u_gauss[i_scheme, :], cfl, scheme)
                its_limiter_fd(nt, res, u_step[i_scheme, :], cfl, scheme)
                its_limiter_fd(nt, res, u_2pw[i_scheme, :], cfl, scheme)
                its_limiter_fd(nt, res, u_4pw[i_scheme, :], cfl, scheme)
            else:
                its_fd(nt, res, u_gauss[i_scheme, :], cfl, scheme)
                its_fd(nt, res, u_step[i_scheme, :], cfl, scheme)
                its_fd(nt, res, u_2pw[i_scheme, :], cfl, scheme)
                its_fd(nt, res, u_4pw[i_scheme, :], cfl, scheme)

        # One plot per cfl
        plot_sim(x_th, x, x0, u_gauss, u_step, u_2pw, u_4pw, 
                schemes, f'CFL = {cfl:.2f} - dx = {dx:.2e} m - dt = {dt:.2e} s - nits = {nt:d}', 
                sim_dir / f'sim_cfl_{index}')

def run_cvg(a: float, cfls: list, xmin: float, xmax: float, nnxs: np.ndarray, 
    n_periods_cvg: float, schemes: list, functions: list, cvg_dir: Path):
    """ Run convergence for given resolutions and specified schemes in 1D domain
    (xmin, xmax) with profiles specified in functions

    :param a: Convection speed
    :type a: float
    :param cfls: CFLs studied
    :type cfls: list
    :param xmin: Left boundary
    :type xmin: float
    :param xmax: Right boundary
    :type xmax: float
    :param nnxs: Resolutions studied
    :type nnxs: np.ndarray
    :param n_periods_cvg: Number of periods
    :type n_periods_cvg: float
    :param schemes: Schemes studied
    :type schemes: list
    :param functions: 1D profiles
    :type functions: list
    :param cvg_dir: Directory for figures
    :type cvg_dir: Path
    """
    print('\n-------------------------------------------------------')
    print(f'Studying mesh convergence')
    print('-------------------------------------------------------')

    Lx = xmax - xmin
    x0 = (xmin + xmax) / 2
    n_schemes = len(schemes)
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

    # Mesh properties
    xmin, xmax = -1, 1

    # Convection speed
    a = 1.0

    # Launch simulations
    cfls = [0.1, 0.3, 0.5, 0.7, 0.9]
    run_sims(a, cfls, xmin, xmax, 401, 2.0, schemes, sim_dir)

    # Plot the diffusion and disperson errors 
    # from the amplification factors of the schemes
    print(f'\n--> Plotting amplifications factors...')
    for scheme in schemes:
        plot_G(scheme, cfls, sp_dir)
    
    # Mesh convergence of the schemes
    functions = ['gaussian(x, x0, 0.3)', 'step(x, x0)', 'packet_wave(x, x0, 0.5)']
    nnxs = np.array([51, 101, 201, 501])
    run_cvg(a, cfls, xmin, xmax, nnxs, 2.0, schemes, functions, cvg_dir)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--schemes', help='Name of the schemes to study', nargs='+')
    parser.add_argument('-d', '--figdir', help='Name of the figures directory')
    args = parser.parse_args()
    main(args)
import numpy as np

def gaussian(x, x0, sigma_x):
    """ Gaussian test function """
    return np.exp(- (x - x0)**2 / sigma_x**2)

def step(x, x0):
    """ Step test function """
    return np.where(abs(x - x0) < 0.5, 1.0, 0.0)

def packet_wave(x, x0, lam):
    """ Packet wave of spatial period lam """
    return np.where(abs(x - x0) < 0.5, np.sin(2 * np.pi / lam * (x - x0)), 0)


if __name__ == '__main__':
    xmin, xmax, nnx = -0.1, 0.1, 201
    Lx = (xmax - xmin)
    dx = Lx / (nnx - 1)

    # Print the residuals around the central cell
    #  -----------|
    # |     | dx  |
    # |-----i-----|
    # |     |     |
    # |-----------|

    # Geometric parameters
    x0 = 0.0
    sigma_x = 2e-2
    xi, xj = 0.0, dx
    sij = dx / 2
    Vi = dx**2
    
    # Values at the nodes
    ui = gaussian(xi, x0, sigma_x)
    uj = gaussian(xj, x0, sigma_x)

    # Advection diffusion parameters
    V = 100.0
    D = 1e-2

    # Upwind flux - central diffusion
    upwind_cflux = V * ui
    central_dflux = - D * (ui - uj) / dx
    total_fou_cd = upwind_cflux + central_dflux
    print(f'upwind = {upwind_cflux:.4e} - central = {central_dflux:.4e}')
    print(f'total_fou_cd = {total_fou_cd:.4e} - total * sij = {total_fou_cd * sij:.4e}')

    # SG scheme flux
    alpha = dx * V / D
    sg_flux = V * (uj - np.exp(alpha) * ui) / (1 - np.exp(alpha))
    print(f'sg_flux = {sg_flux:.4e} - sg_flux * sij = {sg_flux * sij:.4e}')
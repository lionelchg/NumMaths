import numpy as np
from numba import njit

@njit(cache=True)
def bjs(sigma, scheme):
    if scheme == 'FOU':
        ju = 1
        jd = 0
        coeffs = [sigma, 1 - sigma]
    elif scheme == 'C1' or scheme == 'LW':
        ju = 1
        jd = 1
        coeffs = [0.5 * sigma * (sigma + 1), 1 - sigma**2, 0.5 * sigma * (sigma - 1)]
    elif scheme == 'SOU':
        ju = 2
        jd = 0
        coeffs = [sigma * (sigma - 1) / 2, sigma * (2 - sigma), (1 - sigma) * (2 - sigma) / 2]
    elif scheme == 'FR':
        ju = 2
        jd = 1
        coeffs = [sigma * (sigma - 1) / 4, sigma * (5 - sigma) / 4, (1 - sigma) * (4 + sigma) / 4,
                    sigma * (sigma - 1) / 4]
    elif scheme == 'TOS':
        ju = 2
        jd = 1
        coeffs = [sigma * (sigma**2 - 1) / 6, sigma * (2 - sigma) * (sigma + 1) / 2, 
            (2 - sigma) * (1 - sigma**2) / 2, sigma * (sigma - 1) * (2 - sigma) / 6]
    elif scheme == 'C2':
        ju = 2
        jd = 2
        coeffs = [(sigma - 1) * sigma * (sigma + 1) * (sigma + 2) / 24,
                - (sigma - 2) * sigma * (sigma + 1) * (sigma + 2) / 6,
                (sigma - 2) * (sigma - 1) * (sigma + 1) * (sigma + 2) / 4,
                - (sigma - 2) * (sigma - 1) * sigma * (sigma + 2) / 6,
                (sigma - 2) * (sigma - 1) * sigma * (sigma + 1) / 24]
    

    return coeffs, ju, jd

@njit(cache=True)
def advance_fd(res, u, sigma, scheme, coeffs, ju, jd):
    """ Calculate the scheme advancement for a scheme with ju upwind points
    and jd downwind points with periodic boundary conditions """
    res[:] = 0.0
    for i in range(len(u)):
        for j in range(-ju, jd + 1):
            if i + j < 0 or i + j > len(u) - 1:
                index = (i + j) % len(u)
            else:
                index = i + j
            res[i] -= coeffs[ju + j] * u[index]
        res[i] += u[i]

@njit(cache=True)
def its_fd(nt, res, u, sigma, scheme):
    """ Function to do iterations in finite difference formulation """
    coeffs, ju, jd = bjs(sigma, scheme)
    for _ in range(nt):
        advance_fd(res, u, sigma, scheme, coeffs, ju, jd)
        u -= res
import numpy as np
from numba import njit
import limiters as lmt

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
def periodic_index(index, n):
    """ Return the generalized periodic index associated to index """
    if index < 0:
        per_index = (index - 1) % n
    elif index > n - 1:
        per_index = (index + 1) % n
    else:
        per_index = index
    return per_index

@njit(cache=True)
def advance_fd(res, u, coeffs, ju, jd):
    """ Calculate the scheme advancement for a scheme with ju upwind points
    and jd downwind points with periodic boundary conditions """
    res[:] = 0.0
    for i in range(len(u)):
        for j in range(-ju, jd + 1):
            index = periodic_index(i + j, len(u))
            res[i] -= coeffs[ju + j] * u[index]
        res[i] += u[i]

@njit(cache=True)
def its_fd(nt, res, u, sigma, scheme):
    """ Function to do iterations in finite difference formulation """
    coeffs, ju, jd = bjs(sigma, scheme)
    for _ in range(nt):
        advance_fd(res, u, coeffs, ju, jd)
        u -= res

@njit(cache=True)
def grad_ratio(u, i, i1, im1):
    """ Ratio of gradients at location i """
    if u[i1] == u[i]:
        r_i = 0.0
    elif u[i] == u[im1]:
        r_i = 10.0 * np.sign((u[i1] - u[i]))
    else:
        r_i = (u[i1] - u[i]) / (u[i] - u[im1])
    return r_i

@njit(cache=True)
def lim_value(u, i, i1, im1, limiter):
    r = grad_ratio(u, i, i1, im1)
    if r <= 0: return 0.0 
    if limiter == 'van_leer':
        psi_r_i = (r + abs(r)) / (1 + r)
    elif limiter == 'min_mod':
        psi_r_i = min(r, 1)
    elif limiter == 'superbee':
        psi_r_i = max(0, min(2 * r, 1), min(r, 2))
    return psi_r_i

@njit(cache=True)
def advance_limiter_fd(res, u, sigma, limiter):
    """ Calculate the scheme advancement for a limited scheme based on SOU
    WB scheme with periodic boundary conditions """

    for i in range(len(u)):
        res[i] += u[i]
        i1 = periodic_index(i + 1, len(u))
        im1 = periodic_index(i - 1, len(u))
        im2 = periodic_index(i - 2, len(u))

        psi_r_i = lim_value(u, i, i1, im1, limiter)
        psi_r_im1 = lim_value(u, im1, i, im2, limiter)
        res[i] = sigma * (1 + 0.5 * (1 - sigma) * 
                    psi_r_i) * (u[i] - u[im1]) \
                    - 0.5 * sigma * (1 - sigma) * psi_r_im1 \
                    * (u[im1] - u[im2])

@njit(cache=True)
def its_limiter_fd(nt, res, u, sigma, limiter):
    """ Function to do iterations in finite difference formulation """
    for it in range(nt):
        advance_limiter_fd(res, u, sigma, limiter)
        u -= res
#!/Users/cheng/code/envs/dl/bin/python
import os
import numpy as np
import cmath
import matplotlib.pyplot as plt
from fd_schemes import bjs

def ampl_factor(phi, sigma, scheme):
    """ Return the amplification factor for a scheme with ju upwind points
    jd downwind points and coeffs coefficients """
    coeffs, ju, jd = bjs(sigma, scheme)
    G = np.zeros_like(phi, dtype=np.complex128)
    for j in range(-ju, jd + 1):
        G += coeffs[ju + j] * np.exp(j * 1j * phi)
    return G

def errors(G_num, phi, sigma):
    """ Computation of diffusion and dispersion error for a constant advection
    speed problem """
    diff_err = abs(G_num)
    disp_err = np.zeros_like(phi)
    disp_err[0] = 1
    disp_err[1:] = np.array([- cmath.phase(G_num[i]) / sigma / phi[i] for i in range(1, len(phi))])
    return diff_err, disp_err

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from utils import vec_str, print_mat
from pathlib import Path
import copy

def system_energy(A, b, y):
    """ The energy of a linear system """
    return 0.5 * np.dot(y, np.matmul(A, y)) - np.dot(y, b)

def gradient_method(A: np.ndarray, x_0: np.ndarray, b: np.ndarray,
    P: np.ndarray, maxits: int, tol: float):
    """ Implementation of Gradient Method """
    it = 0
    r0 = b - np.matmul(A, x_0)
    norm_r0 = LA.norm(r0)
    if norm_r0 == 0: norm_r0 = 1.0
    err = norm_r0

    # List of iterates for plotting
    list_xk = list()

    # Initialization of temporary variable for solution
    x_k = x_0
    list_xk.append(x_k)
    r_k = b - np.matmul(A, x_k)
    while (err > tol and it <= maxits):
        alpha_k = np.dot(r_k, r_k) / np.dot(r_k, np.matmul(A, r_k))
        x_k = x_k + alpha_k * r_k
        list_xk.append(x_k)
        r_k = b - np.matmul(A, x_k)
        err = LA.norm(r_k) / norm_r0
        it += 1
        if it % 1 == 0:
            print(f'Iteration {it:3d} - x = [{vec_str(x_k)}] - err = {err:.2e}')
    print(f'Last iteration {it:3d} - x = [{vec_str(x_k)}] - err = {err:.2e}')

    return list_xk

def conjugate_gradient_method(A: np.ndarray, x_0: np.ndarray, b: np.ndarray,
    P: np.ndarray, maxits: int, tol: float):
    """ Implementation of Conjugate-Gradient Method """
    it = 0
    r0 = b - np.matmul(A, x_0)
    norm_r0 = LA.norm(r0)
    if norm_r0 == 0: norm_r0 = 1.0
    err = norm_r0

    # List of iterates for plotting
    list_xk = list()

    # Initialization of temporary variable for solution
    x_k = x_0
    list_xk.append(x_k)
    r_k = b - np.matmul(A, x_k)
    p_k = copy.deepcopy(r_k)
    while (err > tol and it <= maxits):
        alpha_k = np.dot(p_k, r_k) / np.dot(p_k, np.matmul(A, p_k))
        x_k = x_k + alpha_k * p_k
        list_xk.append(x_k)
        r_k = b - np.matmul(A, x_k)

        # Update direction of gradient descent
        beta_k = np.dot(p_k, np.matmul(A, r_k)) / np.dot(p_k, np.matmul(A, p_k))
        p_k = r_k - beta_k * p_k
        err = LA.norm(r_k) / norm_r0
        it += 1
        if it % 1 == 0:
            print(f'Iteration {it:3d} - x = [{vec_str(x_k)}] - err = {err:.2e}')
    print(f'Last iteration {it:3d} - x = [{vec_str(x_k)}] - err = {err:.2e}')

    return list_xk

def ax_prop(ax):
    """ Axes object properties """
    ax.grid(True)
    ax.set_aspect('equal')

if __name__ == '__main__':
    # Figures directory
    fig_dir = Path('figures/gradients_method')
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Linear system parameters
    lambda_1 = 2.0
    lambda_2 = 1.0
    A = np.array([[lambda_1, 0],
                [0, lambda_2]])
    x_sol = np.ones(2)
    b = np.matmul(A, x_sol)

    # Gradient method
    x_0 = np.array([4, 3.5])
    list_xk_gm = gradient_method(A, x_0, b, np.eye(2), 5, 1e-4)
    iterate_xks_gm = np.array(list_xk_gm)
    list_energiesk_gm = [system_energy(A, b, xk) for xk in list_xk_gm]

    # Conjugate gradient method
    x_0 = np.array([4, 3.5])
    list_xk_cg = conjugate_gradient_method(A, x_0, b, np.eye(2), 5, 1e-4)
    iterate_xks_cg = np.array(list_xk_cg)
    list_energiesk_cg = [system_energy(A, b, xk) for xk in list_xk_cg]

    # Compute the system energies over a grid
    ny1, ny2 = 101, 101
    y1, y2 = np.linspace(-1, 5, ny1), np.linspace(-1, 5, ny2)
    Y1, Y2 = np.meshgrid(y1, y2)
    Phi = np.zeros_like(Y1)
    for i in range(ny1):
        for j in range(ny2):
            Phi[j, i] = system_energy(A, b, np.array([Y1[j, i], Y2[j, i]]))

    # Global figure
    fig, axes = plt.subplots(ncols=2, figsize=(8, 6))

    # GM plotting
    cs = axes[0].contour(Y1, Y2, Phi, levels=list_energiesk_gm[::-1], colors='k')
    axes[0].clabel(cs, fontsize=9, inline=True)
    axes[0].plot(iterate_xks_gm[:, 0], iterate_xks_gm[:, 1], color='firebrick')
    ax_prop(axes[0])

    # CG plotting
    cs = axes[1].contour(Y1, Y2, Phi, levels=list_energiesk_cg[::-1], colors='k')
    axes[1].clabel(cs, fontsize=9, inline=True)
    axes[1].plot(iterate_xks_cg[:, 0], iterate_xks_cg[:, 1], color='firebrick')
    ax_prop(axes[1])
    fig.savefig(fig_dir / 'system_energy', bbox_inches='tight')


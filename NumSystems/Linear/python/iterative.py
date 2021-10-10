import sys
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from utils import vec_str
from pathlib import Path

def decomp_jacobi(mat):
    """ Return the iteration matrix and inverse of P for Jacobi iterations """
    # Extract diagonal part of mat
    P = np.diag(np.diag(mat))
    N = P - mat
    Pinv = LA.inv(P)
    B = np.matmul(Pinv,  N)
    return B, Pinv

def decomp_JOR(mat, omega):
    """ Return the iteration matrix and inverse of P for Over-relaxation method """
    B_J, Pinv_J = decomp_jacobi(mat)
    B = omega * B_J + (1 - omega) * np.eye(*B_J.shape)
    Pinv = omega * Pinv_J
    return B, Pinv

def decomp_gauss_seidel(mat):
    """ Return the iteration matrix and inverse of P for Gauss-Seidel iterations """
    P = np.tril(mat)
    N = P - mat
    Pinv = LA.inv(P)
    B = np.matmul(Pinv, N)
    return B, Pinv

def decomp_SOR(mat, omega):
    """ Return the iteration matrix and inverse of P for
    Successive over-relaxation (SOR) iterations """
    D = np.diag(np.diag(mat))
    Dinv = LA.inv(D)
    E = - np.tril(mat, k=-1)
    F = - np.triu(mat, k=1)
    B = np.matmul(LA.inv(np.eye(*mat.shape) - omega * np.matmul(Dinv, E)),
        ((1 - omega) * np.eye(*mat.shape) + omega * np.matmul(Dinv, F)))
    Pinv = D / omega - E
    return B, Pinv

current_mod = sys.modules[__name__]

def iterations(mat: np.ndarray, rhs: np.ndarray, x0: np.ndarray,
            tol: float, method_dict: dict):
    """ Iterations for a general x^{(k+1)} = B x^{(k)} + f iteration
    as described in Numerical Mathematics chap. 4 """
    it, maxits = 0, 20
    r0 = rhs - np.matmul(mat, x0)
    norm_r0 = LA.norm(r0)
    err = norm_r0
    B, Pinv = getattr(current_mod, method_dict['name'])(mat, **method_dict['args'])
    x = x0
    while (err > tol and it <= maxits):
        x = np.matmul(B, x) + np.matmul(Pinv, rhs)
        r = rhs - np.matmul(mat, x)
        err = LA.norm(r) / norm_r0
        it += 1
        if it % 5 == 0:
            print(f'Iteration {it:3d} - x = [{vec_str(x)}] - err = {err:.2e}')
    print(f'Last iteration {it:3d} - x = [{vec_str(x)}] - err = {err:.2e}')

def spectal_radius(mat):
    """ Return the spectral radius of matrix mat """
    return np.max(np.abs(LA.eigvals(mat)))

def rho_method(mat, method_dict):
    """ Return the spectral radius of the iteration matrix of specific method """
    B, _ = getattr(current_mod, method_dict['name'])(mat, **method_dict['args'])
    return spectal_radius(B)

def ax_prop(ax):
    ax.grid(True)

if __name__ == '__main__':
    # Figures directory
    fig_dir = Path('figures')
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Different matrices
    A1 = np.array([[3, 0, 4],
                [7, 4, 2],
                [-1, 1, 2]])

    A2 = np.array([[-3, 3, -6],
                [-4, 7, -8],
                [5, 7, 9]])

    A3 = np.array([[4, 1, 1],
                [2, -9, 0],
                [0, -8, -6]])

    A4 = np.array([[7, 6, 9],
                [4, 5, -4],
                [-7, -3, 8]])

    # Theoretical solution of (1, 1, 1)
    x = np.ones(3)

    mats = {'A1': A1, 'A2': A2, 'A3': A3, 'A4': A4}
    # mats = {'A1': A1}

    jacobi_method = {'name': 'decomp_jacobi', 'args': {}}
    gs_method = {'name': 'decomp_gauss_seidel', 'args': {}}
    methods = {'Jacobi': jacobi_method, 'Gauss-Seidel': gs_method}

    for name_mat, mat in mats.items():
        print(f'Matrix {name_mat}:')

        # right hand side for chosen x
        b = np.matmul(mat, x)

        # Use different iterative methods
        for method_name, method in methods.items():
            x0 = np.array([1, 0, 0])
            rho = rho_method(mat, method)
            print(f'rho_{method_name} = {rho:.2e}')
            iterations(mat, b, x0, 1e-3, method)
            print('')

        # Compute spectral radiuses of JOR
        omega_JOR = np.linspace(0.1, 1.0, 20)
        rhos_JOR = np.zeros_like(omega_JOR)
        for i_JOR, omega in enumerate(omega_JOR):
            JOR_method = {'name': 'decomp_JOR', 'args': {'omega': omega}}
            rhos_JOR[i_JOR] = rho_method(mat, JOR_method)

        # Compute spectral radiuses of SOR
        omega_SOR = np.linspace(0.1, 1.9, 36)
        rhos_SOR = np.zeros_like(omega_SOR)
        for i_SOR, omega in enumerate(omega_SOR):
            SOR_method = {'name': 'decomp_SOR', 'args': {'omega': omega}}
            rhos_SOR[i_SOR] = rho_method(mat, SOR_method)

        # Create associated figures
        fig, axes = plt.subplots(ncols=2, figsize=(8, 5))
        axes[0].plot(omega_JOR, rhos_JOR)
        axes[1].plot(omega_SOR, rhos_SOR)
        ax_prop(axes[0])
        ax_prop(axes[1])
        fig.savefig(fig_dir / f'{name_mat}_rhos', bbox_inches='tight')
        plt.close(fig)
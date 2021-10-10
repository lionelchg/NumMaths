import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from utils import vec_str, print_mat

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

    # Initialization of temporary variable for solution
    x_k = x_0
    r_k = b - np.matmul(A, x_k)
    while (err > tol and it <= maxits):
        alpha_k = np.dot(r_k, r_k) / np.dot(r_k, np.matmul(A, r_k))
        x_k = x_k + alpha_k * r_k
        r_k = b - np.matmul(A, x_k)
        err = LA.norm(r_k) / norm_r0
        it += 1
        if it % 1 == 0:
            print(f'Iteration {it:3d} - x = [{vec_str(x_k)}] - err = {err:.2e}')
    print(f'Last iteration {it:3d} - x = [{vec_str(x_k)}] - err = {err:.2e}')

if __name__ == '__main__':
    lambda_1 = 2.0
    lambda_2 = 1.0
    A = np.array([[lambda_1, 0],
                [0, lambda_2]])
    x_sol = np.ones(2)
    b = np.matmul(A, x_sol)
    print_mat(A, 'A')
    print(f'b = {vec_str(b)}')

    # initial solution for testing iterative method
    x_0 = np.array([4, 5])
    gradient_method(A, x_0, b, np.eye(2), 5, 1e-4)

    # Plot the gradient descent
    ny1, ny2 = 101, 101
    y1, y2 = np.linspace(0, 5, ny1), np.linspace(0, 5, ny2)
    Y1, Y2 = np.meshgrid(y1, y2)
    Phi = np.zeros_like(Y1)
    for i in range(ny1):
        for j in range(ny2):
            Phi[j, i] = system_energy(A, b, np.array([Y1[j, i], Y2[j, i]]))
    fig, ax = plt.subplots()
    cs = ax.contour(Y1, Y2, Phi)
    ax.clabel(cs, fontsize=9, inline=True)
    ax.set_aspect('equal')
    fig.savefig('figures/system_energy', bbox_inches='tight')


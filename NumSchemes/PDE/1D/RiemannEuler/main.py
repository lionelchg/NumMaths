import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from exact_riemann import ExactRiemann, pressure_function

# Set default matplotlib params
mpl.rc('lines', linewidth=1.8)

if __name__ == '__main__':
    # Figures directory
    fig_dir = Path('figures')
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Create f_K figure, create left and right state from Test 1 of Toro chapter 4
    gamma = 1.4
    rho_L, u_L, p_L = 1.0, 0.0, 1.0
    rho_R, u_R, p_R = 0.125, 0.0, 0.2

    # Create figure
    p = np.linspace(0.01, 1.2, 201)
    fig, ax = plt.subplots()
    ax.plot(p, pressure_function(p, rho_L, u_L, p_L, rho_R, -2.5, p_R, gamma), label='$\Delta u$ = -2.5')
    ax.plot(p, pressure_function(p, rho_L, u_L, p_L, rho_R, 0.0, p_R, gamma), label='$\Delta u$ = 0')
    ax.plot(p, pressure_function(p, rho_L, u_L, p_L, rho_R, 2.0, p_R, gamma), label='$\Delta u$ = 2')
    ax.plot([p_L, p_L], [-1.0, 1.0], color='k')
    ax.plot([p_R, p_R], [-1.0, 1.0], color='k')
    ax.plot(p, np.zeros_like(p), color='k')
    ax.grid(True)
    ax.set_xlabel('$p$')
    ax.set_ylabel('$f(p)$')
    ax.legend()
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    fig.savefig(fig_dir / 'pressure_function', bbox_inches='tight')

    # Define classes for Test 1 to 5 from Toro chapter 4
    test1 = ExactRiemann(1.0, 0.0, 1.0, 0.125, 0.0, 0.1, 1.4, 'Test 1')
    test2 = ExactRiemann(1.0, -2.0, 0.4, 1.0, 2.0, 0.4, 1.4, 'Test 2')
    test3 = ExactRiemann(1.0, 0.0, 1e3, 1.0, 0.0, 1e-2, 1.4, 'Test 3')
    test4 = ExactRiemann(1.0, 0.0, 1e-2, 1.0, 0.0, 1e2, 1.4, 'Test 4')
    test5 = ExactRiemann(5.99924, 19.5975, 460.894, 5.99242, -6.19633, 46.095, 1.4, 'Test 5')

    # Solve the values
    test1.solve_riemann()
    test2.solve_riemann(initial_guess='TR')
    test3.solve_riemann()
    test4.solve_riemann()
    test5.solve_riemann()

    # Print solutions
    print(test1)
    print(test2)
    print(test3)
    print(test4)
    print(test5)
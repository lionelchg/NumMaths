import numpy as np
import scipy.constants as co
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path
from exact_riemann import ExactRiemann, pressure_function

# Set default matplotlib params
mpl.rc('lines', linewidth=1.8)

if __name__ == '__main__':
    # Figures directory
    fig_dir = Path('figures/avip')
    fig_dir.mkdir(parents=True, exist_ok=True)
    data_dir = Path('data/avip')
    data_dir.mkdir(parents=True, exist_ok=True)

    # Compute pressure
    rho_L = 1.0
    rho_R = 0.01
    Wair = 0.02885066880e0
    Tair = 300.0
    p_L = rho_L * co.R / Wair * Tair
    p_R = rho_R * co.R / Wair * Tair

    # Define classes for avip cases
    case1 = ExactRiemann(rho_L, 0.0, p_L, rho_R, 0.0, p_R, 1.4, 'Test 1')
    case1.plot_pfunction(fig_dir / 'pressure_function_case1')

    # Solve the values
    case1.solve_riemann(initial_guess=1e4)

    # Print solutions
    print(case1)

    # Plot solutions
    x = np.linspace(0.0, 1.0, 401)
    x0 = 0.5
    case1.plot_solution(x, x0, 3e-4, fig_dir / 'case1')
    case1.print_solution(x, x0, 3e-4, data_dir / 'case1.dat')

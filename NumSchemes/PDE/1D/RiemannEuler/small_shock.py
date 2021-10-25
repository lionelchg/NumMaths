# Test cases in chapter 10 of Toro
import numpy as np
import matplotlib as mpl
from pathlib import Path
from exact_riemann import ExactRiemann, pressure_function

# Set default matplotlib params
mpl.rc('lines', linewidth=1.8)

if __name__ == '__main__':
    # Figures directory
    fig_dir = Path('figures/small_shock')
    fig_dir.mkdir(parents=True, exist_ok=True)
    data_dir = Path('data/small_shock')
    data_dir.mkdir(parents=True, exist_ok=True)

    # Define classes for Test 1 to 5 from Toro chapter 4
    test1 = ExactRiemann(1.0, 0.3, 1.0, 0.5, 0.0, 0.3, 1.4, 'Small shock')

    # Solve the values
    test1.solve_riemann()

    # Print solutions
    print(test1)

    # Write solution
    test1.print_info(data_dir / 'summary.md')

    # Plot solutions
    x = np.linspace(0.0, 1.0, 201)
    x0 = 0.3
    test1.plot_solution(x, x0, 0.25, fig_dir / 'test1')
    test1.print_solution(x, x0, 0.25, data_dir / 'test1.dat')

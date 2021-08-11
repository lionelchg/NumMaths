import numpy as np
import matplotlib as mpl
from pathlib import Path
from exact_riemann import ExactRiemann, pressure_function

# Set default matplotlib params
mpl.rc('lines', linewidth=1.8)

if __name__ == '__main__':
    # Figures directory
    fig_dir = Path('figures/toro_tests')
    fig_dir.mkdir(parents=True, exist_ok=True)
    data_dir = Path('data/toro_tests')
    data_dir.mkdir(parents=True, exist_ok=True)

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

    # Plot solutions
    x = np.linspace(-0.5, 0.5, 401)
    test1.plot_solution(x, 0.25, fig_dir / 'test1')
    test1.print_solution(x, 0.25, data_dir / 'test1.dat')
    test2.plot_solution(x, 0.15, fig_dir / 'test2')
    test2.print_solution(x, 0.15, data_dir / 'test2.dat')
    test3.plot_solution(x, 0.012, fig_dir / 'test3')
    test3.print_solution(x, 0.012, data_dir / 'test3.dat')
    test4.plot_solution(x, 0.035, fig_dir / 'test4')
    test4.print_solution(x, 0.035, data_dir / 'test4.dat')
    test5.plot_solution(x, 0.035, fig_dir / 'test5')
    test5.print_solution(x, 0.035, data_dir / 'test5.dat')
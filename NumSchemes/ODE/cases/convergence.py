# Standard imports
import numpy as np
import matplotlib.pyplot as plt

# Library import
from odesolver.utils import create_dir
from odesolver.solver import ODESim

# Local
from rhssquare import RHSSquare

class ConvergenceODE:
    def __init__(self, tmin, tend, list_nts, schemes, model, init_value, fig_dir):
        self.tmin = tmin
        self.tend = tend
        self.list_nts = list_nts
        self.times = [np.linspace(tmin, tend, ntimes) for ntimes in list_nts]
        self.sims = [ODESim(times, schemes, model, init_value) for times in self.times]
        self.model = model
        self.schemes = schemes
        self.linestyles = ['k-.', 'k--', 'k:']
        self.fig_dir = f'figures/{fig_dir}/'
        create_dir(self.fig_dir)

    def run_convergence(self):
        for sim in self.sims:
            sim.run_schemes()
    
    def plot_errors(self):
        for i_scheme, scheme in enumerate(self.schemes):
            fig, axes = plt.subplots(nrows=2, sharex=True, figsize=(8, 8))
            for i_sim, sim in enumerate(self.sims):
                u_exact = self.model.u_exact(sim.times)
                self.model.ax_plot(axes, sim.times, sim.v[i_scheme, :], u_exact, self.linestyles[i_sim])
            axes[0].plot(sim.times, u_exact, 'k')
            axes[0].legend([rf'$\Delta x$ = {(self.tend - self.tmin) / (nts - 1):.1e}' for nts in self.list_nts] + ['Exact'])
            fig.suptitle(f'Convergence for {scheme}')
            fig.tight_layout(rect=[0, 0.03, 1, 0.97])
            fig.savefig(self.fig_dir + scheme, bbox_inches='tight')

if __name__ == '__main__':
    # Convergence test for u(t) = 1 / (1 + t)
    tmin, tend = 0, 10
    list_nts = [101, 201, 401]
    cvg_sim = ConvergenceODE(tmin, tend, list_nts, ['forwardEuler', 'midpoint', 'multi_step2'], RHSSquare(), 1.0, 'cvg/')
    cvg_sim.run_convergence()
    cvg_sim.plot_errors()
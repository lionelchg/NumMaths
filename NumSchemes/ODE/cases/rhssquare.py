import numpy as np
import matplotlib.pyplot as plt
from odesolver.solver import ODESim

class RHSSquare:
    def __init__(self):
        # Dimension of the system
        self.nd = 1

    def f(self, u, t):
        return - u**2
    
    def u_exact(self, time):
        return 1 / (1 + time)
    
    def ax_plot(self, axes, time, u, u_exact, linestyle):
        axes[0].plot(time, u, linestyle)
        axes[0].set_xlabel('$t$ [s]')
        axes[1].plot(time, np.abs(u[:, 0] - u_exact), linestyle)
    
    def plot(self, time, u, figtitle, figname):
        fig, axes = plt.subplots(nrows=2, sharex=True)
        u_exact = self.u_exact(time)
        axes[0].plot(time, u_exact, 'k--')
        self.ax_plot(axes, time, u, u_exact, 'k')
        axes[0].legend(['Exact', 'Simulation'])
        fig.suptitle(figtitle)
        fig.tight_layout(rect=[0, 0.03, 1, 0.97])
        fig.savefig(figname, bbox_inches='tight')

if __name__ == '__main__':
    # Second test with other model
    tmin, tend, ntimes = 0, 10, 101
    times = np.linspace(tmin, tend, ntimes)
    sim = ODESim(times, ['forwardEuler', 'midpoint', 'multi_step2'], RHSSquare(), 1.0, fig_dir='rhs_square/')
    sim.run_schemes()
    sim.plot()
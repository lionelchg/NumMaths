import numpy as np
import matplotlib.pyplot as plt
from odesolver.solver import ODESim

class Pendulum:
    def __init__(self, L, g):
        # Parameters of the model
        self.L = L
        self.g = g

        # Dimension of the system
        self.nd = 2
    
    def f(self, u, t):
        return np.array([-self.g / self.L * np.sin(u[1]), u[0]])
    
    def plot(self, time, u, figtitle, figname):
        fig, axes = plt.subplots(nrows=2, sharex=True)
        axes[0].plot(time, u[:, 1], 'k')
        axes[0].set_ylabel(r'$\theta$')
        axes[1].plot(time, u[:, 0], 'k')
        axes[1].set_ylabel(r'$\theta_t$')
        axes[1].set_xlabel('$t$ [s]')
        fig.suptitle(figtitle)
        fig.tight_layout(rect=[0, 0.03, 1, 0.97])
        fig.savefig(figname, bbox_inches='tight')

if __name__ == '__main__':
    # Pendulum system model
    tmin, tend, ntimes = 0, 10, 501
    times = np.linspace(tmin, tend, ntimes)
    sim = ODESim(times, ['forwardEuler', 'midpoint'], Pendulum(1, 9.81), 
        np.array([0.0, 45 * np.pi / 180]), fig_dir='pendulum/')
    sim.run_schemes()
    sim.plot()
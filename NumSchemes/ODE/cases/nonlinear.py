import numpy as np
import matplotlib.pyplot as plt
from odesolver.solver import ODESim
from odesolver.utils import make_times

class NonLinear:
    def __init__(self):
        self.nd = 2

    def f(self, u, t):
        return np.array([np.sin(u[1]), np.cos(u[0])])

    def jac_f(self, u, t):
        """ Jacoabian of the function """
        return np.array([
            [0, np.cos(u[1])],
            [- np.sin(u[0]), 0]
        ])

    def R_trapez(self, w, v, t, dt):
        """ Residual for trapezoidal method """
        return w - v - 0.5 * dt * (self.f(w, t + dt) + self.f(v, t))
    
    def jac_R(self, u, t, dt):
        """ Jacobian of the algebraic function verified by the
        trapezoidal method """
        return np.eye(self.nd) - dt / 2 * self.jac_f(u, t)
    
    @staticmethod
    def norm2(u):
        return np.sqrt(np.sum(u**2))
    
    def ftrapez(self, v, t, dt):
        eps = 1.0e-5
        tmp_w = v
        while self.norm2(self.R_trapez(tmp_w, v, t, dt)) >= eps:
            dw = - np.matmul(np.linalg.inv(self.jac_R(tmp_w, t, dt)), self.R_trapez(tmp_w, v, t, dt))
            tmp_w += dw
        print(tmp_w)
        return tmp_w
    
    def plot(self, time, u, figtitle, figname):
        fig, axes = plt.subplots(nrows=2, sharex=True)
        axes[0].plot(time, u[:, 0], 'k')
        axes[0].set_ylabel(r'$u_1$')
        axes[1].plot(time, u[:, 1], 'k')
        axes[1].set_ylabel(r'$u_2$')
        axes[1].set_xlabel('$t$ [s]')
        fig.suptitle(figtitle)
        fig.tight_layout(rect=[0, 0.03, 1, 0.97])
        fig.savefig(figname, bbox_inches='tight')
    
if __name__ == '__main__':
    # Pendulum system model
    tmin, tend, ntimes = 0, 0.1, 2
    times = np.linspace(tmin, tend, ntimes)
    model = NonLinear()
    sim = ODESim(times, ['trapezoidal'], model, 
        np.array([0.0, np.pi / 2]), fig_dir='nonlinear/')
    sim.run_schemes()
    print(sim.v0)
    print(sim.v[0, 1, :])
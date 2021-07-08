import numpy as np
import matplotlib.pyplot as plt
from .utils import create_dir

class ODESim:
    def __init__(self, times, schemes, model, init_value, fig_dir=None):
        # Time related variables
        self.times = times
        self.ntimes = len(times)
        self.dt = times[1] - times[0]

        # Model and schemes
        self.model = model
        self.schemes = schemes
        self.nschemes = len(schemes)

        # variable of interest
        self.v = np.zeros((self.nschemes, self.ntimes, self.model.nd))
        self.v0 = init_value

        # Figures directory
        if not fig_dir is None:
            self.fig_dir = f'figures/{fig_dir}/'
            create_dir(self.fig_dir)
    
    def forwardEuler(self, v):
        """ Use model to apply forward Euler scheme """
        v[0] = self.v0
        for i in range(1, self.ntimes):
            v[i] = v[i - 1] + self.dt * self.model.f(v[i - 1], self.times[i - 1])
    
    def midpoint(self, v):
        """ Use model to apply midpoint formula """
        v[0] = self.v0
        v[1] = v[0] + self.dt * self.model.f(v[0], self.times[0])
        for i in range(2, self.ntimes):
            v[i] = v[i - 2] + 2 * self.dt * self.model.f(v[i - 1], self.times[i - 1])
    
    def multi_step2(self, v):
        """ Most accurate explicit 2multistep method """
        v[0] = self.v0
        v[1] = v[0] + self.dt * self.model.f(v[0], self.times[0])
        for i in range(2, self.ntimes):
            v[i] = - 4 * v[i - 1] + 5 * v[i - 1] + self.dt * \
                        (4 * self.model.f(v[i - 1], self.times[i - 1]) + 2 * self.model.f(v[i - 2], self.times[i - 2]))
    
    def backwardEuler(self, v):
        """ Only works with stiff problem for now """
        v[0] = self.v0
        for i in range(1, self.ntimes):
            v[i] = self.model.fbackwardEuler(v[i - 1], self.times[i - 1], self.dt)
        
    def trapezoidal(self, v):
        """ Only works with nonlinear problem for now """
        v[0] = self.v0
        for i in range(1, self.ntimes):
            v[i] = self.model.ftrapez(v[i - 1], self.times[i - 1], self.dt)
        print(v.shape)
    
    def run_schemes(self):
        """ Apply scheme and plot the results """
        for i_scheme, name_scheme in enumerate(self.schemes):
            scheme = getattr(self, name_scheme)
            scheme(self.v[i_scheme, :])
    
    def plot(self, figname=None):
        for i_scheme, name_scheme in enumerate(self.schemes):
            if figname is None:
                self.model.plot(self.times, self.v[i_scheme, :], f'{name_scheme} - dt = {self.dt:.2e}', self.fig_dir + name_scheme)
            else:
                self.model.plot(self.times, self.v[i_scheme, :], f'{name_scheme} - dt = {self.dt:.2e}', self.fig_dir + figname)

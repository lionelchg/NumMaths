from os import lseek
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns

def smooth_shock(x, t_s):
    """ Smooth solution turning into a shock at t = t_s """
    return np.sin(2 * np.pi * x) / (2 * np.pi * t_s)

def smooth_rarefaction(x, A, k):
    """ Smooth rarefaction profile """
    return A * np.tanh(k * (x - 0.5))

def compression_shock(x):
    """ Compression turning into a shock """
    return np.where(x <= 0, 1.5, np.where(x >= 1, - 1/2, 1.5 - 2 * x))

def ax_prop(ax, xlabel, ylabel, legend=True):
    ax.grid(True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if legend: ax.legend()

def plot_profiles(x:np.ndarray, profiles: dict, t_s: float, figname: Path):
    # Plotting
    fig, axes = plt.subplots(nrows=2, figsize=(6, 8), sharex=True)
    for time, profile in profiles.items():
        axes[0].plot(x, profile, label=f'{time:.2f}')
    ax_prop(axes[0], '$x$', '$u$')

    # Creation of characteristic from initial profile
    init_profile = profiles[0.0]
    nnx = len(x)
    charac_indices = [int(0.1 * (nnx - 1) * i) for i in range(11)]
    for i_charac in charac_indices:
        x_charac = x[i_charac]
        u_charac = init_profile[i_charac]
        t_tmp = np.linspace(0, 1.3 * t_s, 121)
        x_tmp = x_charac + u_charac * t_tmp
        axes[1].plot(x_tmp, t_tmp, 'k')
    ax_prop(axes[1], '$x$', '$t$')
    
    # Formatting of figure and saving
    fig.tight_layout()
    fig.savefig(figname, bbox_inches='tight')

def Linf_norm(y):
    return np.max(np.abs(y))

def burgers_sol(x, fprofile, args_prof, t):
    """ Give the Burgers solution at instant t of fprofile 
    function """
    tol = 1e-6
    prof = fprofile(x, **args_prof)
    tmp_prof = np.zeros_like(prof)
    err = Linf_norm(prof)
    it = 0
    maxits = 200
    while err > tol and it < maxits:
        tmp_prof[:] = prof[:]
        prof = fprofile(x - prof * t, **args_prof)
        err = Linf_norm(prof - tmp_prof)
        it += 1
    return prof

def run_burgers(x:np.ndarray, fprofile, args_profile: dict, 
            time_shock: float, figname: Path):
    """ Run Burgers solutions for fprofile on 1D mesh x """
    profiles = dict()
    profiles[0.0] = fprofile(x, **args_profile)
    instants = np.array([0.25, 0.5, 0.75, 0.9]) * time_shock
    for time in instants:
        profiles[time] = burgers_sol(x, fprofile, args_profile, time)
    plot_profiles(x, profiles, time_shock, figname)

if __name__ == '__main__':
    # Figures directory
    fig_dir = Path('figures')
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Geometry
    xmin, xmax, nnx = 0.0, 1.0, 201
    x = np.linspace(xmin, xmax, nnx)

    # Profile creation
    rarefaction = smooth_rarefaction(x, 1.0, 5.0)
    comp_shock = compression_shock(x)

    # Advance shock
    time_shock = 2.0
    args_shock = {'t_s': time_shock}
    run_burgers(x, smooth_shock, args_shock, 
            time_shock, fig_dir / 'smooth_shock')

    # Rarefaction
    time_shock = 0.3
    args_rare = {'A': 0.2, 'k': 20.0}
    run_burgers(x, smooth_rarefaction, args_rare,
            time_shock, fig_dir / 'smooth_rarefaction')
    
    # Compression-shock
    xmin, xmax, nnx = -0.5, 1.5, 401
    x = np.linspace(xmin, xmax, nnx)
    time_shock = 0.5
    args_comp = {}
    run_burgers(x, compression_shock, args_comp,
            time_shock, fig_dir / 'compression_shock')
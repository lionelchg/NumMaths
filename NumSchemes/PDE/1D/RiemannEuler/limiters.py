import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def van_leer_waf(r, abs_c):
    """ Van Leer limiter function for WAF formulation """
    return np.where(r <= 0, 1.0, 1 - (1 - abs_c) * 2 * r / (1 + r))

def van_albada_waf(r, abs_c):
    """ Van Albada limiter function for WAF formulation """
    return np.where(r <= 0, 1.0, 1 - (1 - abs_c) * r * (1 + r) / (1 + r**2))

def minbee_waf(r, abs_c):
    return np.where(r <= 0, 1.0, np.where(r <= 1.0, 1 - (1 - abs_c)*r, abs_c))

if __name__ == '__main__':
    fig_dir = Path('figures')
    fig_dir.mkdir(parents=True, exist_ok=True)
    cs = [1.0e-4, 0.2, 0.5, 0.8]
    r = np.linspace(0, 2, 101)
    fig, axes = plt.subplots(ncols=len(cs), figsize=(12, 5))
    for i, c in enumerate(cs):
        ax = axes[i]
        ax.plot(r, van_leer_waf(r, abs(c)), label='Van Leer')
        ax.plot(r, van_albada_waf(r, abs(c)), label='Van Albada')
        ax.plot(r, minbee_waf(r, abs(c)), label='Minbee')
        ax.grid(True)
        ax.legend()
        ax.set_title(f'c = {c:.1e}')
    fig.savefig(fig_dir / 'limiters_waf')
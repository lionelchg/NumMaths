# Plotting 2D colormaps of amplification factors
import os
import numpy as np
import argparse

from plot import plot_G, plot_G_2D
from pathlib import Path

if __name__ == '__main__':
    # figures directory
    fig_dir = Path(f'figures/ampl_factors/')

    # Create the directories
    fig_dir.mkdir(parents=True, exist_ok=True)
    
    # Schemes selected
    schemes = ['FOU', 'LW', 'SOU', 'FR', 'TOS']
    print(f"Schemes selected: {' '.join(schemes)}")

    # Plot the diffusion and disperson errors 
    # from the amplification factors of the schemes
    print(f'\n--> Plotting amplifications factors...')
    for scheme in schemes:
        plot_G(scheme, [0.1, 0.3, 0.5, 0.7, 0.9], fig_dir)
        plot_G_2D(scheme, fig_dir)

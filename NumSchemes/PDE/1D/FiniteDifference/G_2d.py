# Plotting 2D colormaps of amplification factors
import os
import numpy as np
import argparse

from plot import plot_G_2D
from pathlib import Path

def main(args):
    # figures directory
    fig_dir = Path(f'figures/{args.figdir}/')

    # Create the directories
    fig_dir.mkdir(parents=True, exist_ok=True)
    
    # Schemes selected
    schemes = args.schemes
    print(f"Schemes selected: {' '.join(schemes)}")

    # Plot the diffusion and disperson errors 
    # from the amplification factors of the schemes
    print(f'\n--> Plotting amplifications factors...')
    for scheme in schemes:
        plot_G_2D(scheme, fig_dir)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--schemes', help='Name of the schemes to study', nargs='+')
    parser.add_argument('-d', '--figdir', help='Name of the figures directory')
    args = parser.parse_args()
    main(args)
import numpy as np
import matplotlib.pyplot as plt
import argparse
from parser import read_data
from pathlib import Path
import re

def plot(data_fns, figname):
    header, data = read_data(data_fns[0], head=True)
    ndata = len(header)

    # Determine number of columns with fixed rows at 2
    nprofiles = ndata - 1
    nrows = 2
    ncols = nprofiles // nrows + (nprofiles % nrows != 0)

    # Figure creation
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4 * nrows, 4 * ncols))
    axes = axes.reshape(-1)
    
    # Fill the figure
    for data_fn in data_fns:
        data = read_data(data_fn)
        scheme = re.search("/(\w*)\.dat", data_fn).group(1)
        for i in range(1, ndata):
            axes[i - 1].plot(data[:, 0], data[:, i], label=scheme)    

    # Ax properties
    for i in range(1, ndata):
        axes[i - 1].grid(True)
        axes[i - 1].legend()

    fig.savefig(figname, bbox_inches='tight')

if __name__ == '__main__':
    fig_dir = Path('figures/')
    fig_dir.mkdir(parents=True, exist_ok=True)
    args = argparse.ArgumentParser(description="Post-process 1D FV data")
    args.add_argument("-d", "--data_fn", type=str, nargs='+',
            required=True, help="Data filename")
    args.add_argument("-f", "--figname", type=str, 
            required=True, help="Data filename")
    args = args.parse_args()

    plot(args.data_fn, fig_dir / args.figname)


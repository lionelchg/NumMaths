import numpy as np
import matplotlib.pyplot as plt
import argparse
from parser import read_data
from pathlib import Path

def plot(data_fn, figname):
    header, data = read_data(data_fn, head=True)
    ndata = len(header)
    fig, ax = plt.subplots()
    for i in range(1, ndata):
        ax.plot(data[:, 0], data[:, i], label=header[i])
    ax.grid(True)
    ax.legend()
    fig.savefig(figname, bbox_inches='tight')

if __name__ == '__main__':
    fig_dir = Path('figures/')
    fig_dir.mkdir(parents=True, exist_ok=True)
    args = argparse.ArgumentParser(description="Post-process 1D FV data")
    args.add_argument("-d", "--data_fn", type=str, 
            required=True, help="Data filename")
    args.add_argument("-f", "--figname", type=str, 
            required=True, help="Data filename")
    args = args.parse_args()

    plot(args.data_fn, fig_dir / args.figname)


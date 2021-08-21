import h5py
import numpy as np
import matplotlib.pyplot as plt
import argparse
import re
from pathlib import Path

def ax_prop_sim(ax):
    ax.legend()
    ax.grid(True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Postprocessing of HDF5 files from C FiniteVolume code')
    parser.add_argument('-d', '--data_fn', type=str, required=True,
                help='Name of the HDF5 filename')
    args = parser.parse_args()

    h5file = h5py.File(args.data_fn, 'r')
    casename = re.search("/(\w*)\.h5", args.data_fn).group(1)
    fig_dir = Path('figures') / casename
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Common x vector
    x = np.array(h5file['x'])
    cfls = [cfl_str for cfl_str in h5file if cfl_str != 'x']

    for cfl_str in cfls:
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(9, 9), sharex=True)
        axes = axes.reshape(-1)
        for scheme_str in h5file[cfl_str]:
            dset = h5file[cfl_str][scheme_str]
            data = np.array(dset)
            scheme_name = dset.attrs['scheme'].decode('ascii')
            for i, ax in enumerate(axes):
                ax.plot(x, data[i, :], label=scheme_name)

        for ax in axes:
            ax_prop_sim(ax)
        
        fig.savefig(fig_dir / cfl_str, bbox_inches='tight')
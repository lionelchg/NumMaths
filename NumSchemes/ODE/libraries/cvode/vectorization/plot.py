import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from itertools import product
import os

from run import nloops, use_jac

data = dict()
data["true"] = { "serial": list(), "vect": list() }
data["false"] = { "serial": list(), "vect": list() }

for nloop, jac in tqdm(product(nloops, use_jac), total=len(nloops)*len(use_jac)):
    with open(f"out/out_{nloop}_{jac}.txt", "r") as file:
        lines = file.readlines()
    first = True
    for line in lines:
        if first and line.startswith("time_total"):
            data[jac]["serial"].append(float(line.split(" ")[1]))
            first = False
            # print("Found serial in ", nloop, jac)
        elif not first and line.startswith("time_total"):
            data[jac]["vect"].append(float(line.split(" ")[1]))
            # print("Found vect in ", nloop, jac)

# Convert to Numpy
for jac in data:
    for method in data[jac]:
        data[jac][method] = np.asfarray(data[jac][method])

# Plot
fig, axarr = plt.subplots(1, 2, figsize=(10,4))
ax, ax2 = axarr.ravel()

for jac in data:
    for method in data[jac]:
        ax.loglog(nloops, data[jac][method], label=f"{method} jac {jac}")
    ax2.semilogx(nloops, data[jac]["serial"] / data[jac]["vect"], label=f"jac {jac}")
ax.legend()
ax.set_xlabel("nloops")
ax.set_ylabel("t [s]")
ax.set_title("Serial and vectorized solver performance")
ax2.legend()
ax2.set_xlabel("nloops")
ax2.set_ylabel("speedup")
ax2.set_title("Vectorized vs serial speedup")

fig.tight_layout()
fig.savefig("comp.png", dpi=150)
os.system("imgcat comp.png")

import os
from itertools import product
from subprocess import run
from tqdm import tqdm

nloops = [
    "1", "2", "5", "10", "20", "40", "60", "80", "100", "150", "200", "300", "400", "500",
    "650", "800", "1000", "1500", "2000"
]
# nloops = [
#     "1", "2", "5", "10", "20", "40"
# ]
use_jac = ["true", "false"]

if __name__ == "__main__":
    for nloop, jac in tqdm(product(nloops, use_jac), total=len(nloops)*len(use_jac)):
        cmd = [
            "./brusselator_vect",
            nloop,
            jac,
            ">",
            f"out/out_{nloop}_{jac}.txt",
        ]
        run(" ".join(cmd), shell=True)
